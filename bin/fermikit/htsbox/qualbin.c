#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "sam.h"
#include "kstring.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	char *name, *seq;
	uint8_t *qual;
	bam1_t *b;
} qz_record_t;

typedef struct {
	kseq_t *fq;
	BGZF *bam;
	bam_hdr_t *hdr;
	bam1_t *b;
} qz_infile_t;

typedef struct {
	kstring_t str;
	BGZF *fp;
} qz_outfile_t;

static qz_infile_t *qz_open(const char *fn, int is_bam)
{
	qz_infile_t *f;
	f = (qz_infile_t*)calloc(1, sizeof(qz_infile_t));
	if (!is_bam) {
		gzFile fp;
		fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
		f->fq = kseq_init(fp);
	} else {
		f->bam = fn && strcmp(fn, "-")? bgzf_open(fn, "r") : bgzf_dopen(fileno(stdin), "r");
		f->hdr = bam_hdr_read(f->bam);
		f->b = bam_init1();
	}
	return f;
}

static void qz_close(qz_infile_t *f)
{
	if (f->bam) {
		bam_destroy1(f->b);
		bam_hdr_destroy(f->hdr);
		bgzf_close(f->bam);
	}
	if (f->fq) {
		gzclose(f->fq->f->f);
		kseq_destroy(f->fq);
	}
	free(f);
}

static int qz_read(qz_infile_t *f, qz_record_t *r)
{
	int ret = -1, i;
	if (f->fq) {
		ret = kseq_read(f->fq);
		if (ret >= 0) {
			r->name = strdup(f->fq->name.s);
			r->seq = strdup(f->fq->seq.s);
			r->qual = 0;
			if (f->fq->qual.l) {
				r->qual = (uint8_t*)strdup(f->fq->qual.s);
				for (i = 0; i < f->fq->qual.l; ++i)
					r->qual[i] -= 33;
			}
		}
	} else if (f->bam) {
		ret = bam_read1(f->bam, f->b);
		if (ret >= 0) {
			r->b = (bam1_t*)calloc(1, sizeof(bam1_t));
			bam_copy1(r->b, f->b);
			r->seq = 0;
			r->qual = bam_get_qual(r->b);
		}
	}
	return ret;
}

static void qz_write_free(qz_outfile_t *f, qz_record_t *r)
{
	if (r->b) {
		bam_write1(f->fp, r->b);
		bam_destroy1(r->b);
	} else {
		f->str.l = 0;
		kputc(r->qual? '@' : '>', &f->str);
		kputs(r->name, &f->str); kputc('\n', &f->str);
		kputs(r->seq,  &f->str);
		if (r->qual) {
			int i, l;
			l = r->b? r->b->core.l_qseq : strlen(r->seq);
			kputsn("\n+\n", 3, &f->str);
			for (i = 0; i < l; ++i) r->qual[i] += 33;
			kputs((char*)r->qual, &f->str);
		}
		puts(f->str.s);
		free(r->name); free(r->seq); free(r->qual);
	}
}

static void qz_alter(qz_record_t *r, int mode)
{
	int i, l;
	uint8_t *qual = r->qual;
	l = r->b? r->b->core.l_qseq : strlen(r->seq);
	if (mode == 3) {
		for (i = 0; i < l; ++i)
			qual[i] = qual[i] >= 40? 40 : qual[i] >= 20? 30 : 10;
	} else if (mode == 2) {
		for (i = 0; i < l; ++i)
			qual[i] = qual[i] >= 20? 30 : 10;
	} else if (mode == 1) {
		for (i = 0; i < l; ++i) qual[i] = 25;
	} else if (mode == 0) {
		if (!r->b) {
			free(r->qual); r->qual = 0;
		} else for (i = 0; i < l; ++i) qual[i] = 255;
	} else if (mode == 7) { // see Table 1 in "Reducing Whole-Genome Data Storage Footprint"
		for (i = 0; i < l; ++i) {
			if (qual[i] < 2) {} // do nothing
			else if (qual[i] < 10) qual[i] = 6;
			else if (qual[i] < 20) qual[i] = 15;
			else if (qual[i] < 25) qual[i] = 22;
			else if (qual[i] < 30) qual[i] = 27;
			else if (qual[i] < 40) qual[i] = 37;
			else qual[i] = 40;
		}
	}
}

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	qz_infile_t *in;
	qz_outfile_t *out;
	int is_bam, chunk_size, n_threads, mode;
} pipeline_t;

typedef struct {
	long n_rec;
	qz_record_t *rec;
	pipeline_t *p;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
	step_t *step = (step_t*)_data;
	qz_record_t *r = &step->rec[i];
	qz_alter(r, step->p->mode);
}

static void *worker_pipeline(void *shared, int step, void *in) // kt_pipeline() callback
{
	pipeline_t *p = (pipeline_t*)shared;
	if (step == 0) { // step 0: read lines into the buffer
		step_t *s;
		s = (step_t*)calloc(1, sizeof(step_t));
		s->rec = (qz_record_t*)calloc(p->chunk_size, sizeof(qz_record_t));
		s->p = p;
		for (s->n_rec = 0; s->n_rec < p->chunk_size; ++s->n_rec)
			if (qz_read(p->in, &s->rec[s->n_rec]) < 0) break;
		if (s->n_rec) return s;
		free(s->rec); free(s);
	} else if (step == 1) { // step 1: reverse lines
		kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_rec);
		return in;
	} else if (step == 2) { // step 2: write the buffer to output
		step_t *s = (step_t*)in;
		int i;
		for (i = 0; i < (int)s->n_rec; ++i)
			qz_write_free(p->out, &s->rec[i]);
		free(s->rec); free(s);
	}
	return 0;
}

static int usage(FILE *out)
{
	fprintf(stderr, "Usage: qualbin [options] <in.bam>|<in.fq.gz>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -t INT    number of threads [1]\n");
	fprintf(stderr, "  -n INT    number of records in buffer [1000000]\n");
	fprintf(stderr, "  -m INT    number of bins (0, 1, 2, 3 or 7) [2]\n");
	fprintf(stderr, "  -b        the input is a BAM file\n");
	fprintf(stderr, "  -u        output uncompressed BAM (force -b)\n");
	return out == stdout? 0 : 1;
}

int main_qualbin(int argc, char *argv[])
{
	int c, is_un_bam = 0;
	pipeline_t p;
	memset(&p, 0, sizeof(pipeline_t));
	p.mode = 2; p.n_threads = 1; p.chunk_size = 1000000;
	while ((c = getopt(argc, argv, "udt:bn:m:")) >= 0) {
		if (c == 'm') p.mode = atoi(optarg);
		else if (c == 'b') p.is_bam = 1;
		else if (c == 'u') p.is_bam = is_un_bam = 1;
		else if (c == 't') p.n_threads = atoi(optarg);
		else if (c == 'n') p.chunk_size = atoi(optarg);
	}
	if (optind + 1 > argc) return usage(stderr);

	p.in = qz_open(argv[optind], p.is_bam);
	p.out = (qz_outfile_t*)calloc(1, sizeof(qz_outfile_t));
	if (p.is_bam) {
		char mode[3] = "w";
		if (is_un_bam) strcat(mode, "u");
		p.out->fp = bgzf_dopen(fileno(stdout), mode);
		bgzf_mt(p.out->fp, 2, 256);
		bam_hdr_write(p.out->fp, p.in->hdr);
	}
	kt_pipeline(2, worker_pipeline, &p, 3);
	if (p.is_bam) bgzf_close(p.out->fp);
	else free(p.out->str.s);
	free(p.out);
	qz_close(p.in);
	return 0;
}
