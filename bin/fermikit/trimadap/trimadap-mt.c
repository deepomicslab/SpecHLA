#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <assert.h>
#include <zlib.h>
#include "ksw.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define VERSION "r9"

/***************
 * CMD options *
 ***************/

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct {
	int type, len;
	uint8_t *seq;
	uint64_t cnt;
} ta_adap_t;

typedef struct {
	int sa, sb, go, ge; // scoring; not exposed to the command line, for now
	int min_sc, min_len, min_trim_len;
	int n_threads;
	int chunk_size;
	double max_diff;
	int n_adaps, m_adaps;
	ta_adap_t *adaps;
	int8_t mat[25];
	kseq_t *ks;
} ta_opt_t;

void ta_opt_set_mat(int sa, int sb, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? sa : -sb;
		mat[k++] = 0; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = 0;
}

void ta_opt_init(ta_opt_t *opt)
{
	memset(opt, 0, sizeof(ta_opt_t));
	opt->sa = 1, opt->sb = 2, opt->go = 1, opt->ge = 3;
	opt->min_sc = 15, opt->min_len = 8, opt->min_trim_len = 1000000, opt->max_diff = .15f;
	opt->n_threads = 1, opt->chunk_size = 10000000;
	ta_opt_set_mat(opt->sa, opt->sb, opt->mat);
}

void ta_opt_add_adap(ta_opt_t *opt, int type, const char *adap)
{
	ta_adap_t *p;
	uint8_t *q;
	if (opt->m_adaps == opt->n_adaps) {
		opt->m_adaps = opt->m_adaps? opt->m_adaps<<1 : 4;
		opt->adaps = realloc(opt->adaps, opt->m_adaps * sizeof(ta_adap_t));
	}
	p = &opt->adaps[opt->n_adaps++];
	p->len = strlen(adap);
	assert(p->len * opt->sa < 256);
	p->seq = (uint8_t*)strdup(adap);
	p->type = type;
	p->cnt = 0;
	for (q = p->seq; *q; ++q)
		*q = seq_nt4_table[*q];
}

void ta_opt_default_adaps(ta_opt_t *opt)
{
	if (opt->n_adaps) return;
	ta_opt_add_adap(opt, 5, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT");
	ta_opt_add_adap(opt, 3, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
	ta_opt_add_adap(opt, 3, "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT");
	ta_opt_add_adap(opt, 3, "ATCTCGTATGCCGTCTTCTGCTTG");
}

void ta_opt_open(ta_opt_t *opt, const char *fn)
{
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
	opt->ks = kseq_init(fp);
}

void ta_opt_free(ta_opt_t *opt)
{
	int i;
	gzFile fp = opt->ks->f->f;
	for (i = 0; i < opt->n_adaps; ++i)
		free(opt->adaps[i].seq);
	kseq_destroy(opt->ks);
	gzclose(fp);
	free(opt->adaps);
}

/*************
 * FASTQ I/O *
 *************/

typedef struct {
	int l_seq;
	char *name, *seq, *qual;
} bseq1_t;

bseq1_t *bseq_read(kseq_t *ks, int chunk_size, int *n_)
{
	int size = 0, m, n;
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->name = strdup(ks->name.s);
		s->seq = strdup(ks->seq.s);
		s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
		s->l_seq = ks->seq.l;
		size += seqs[n++].l_seq;
		if (size >= chunk_size) break;
	}
	*n_ = n;
	return seqs;
}

/****************
 * Trim one seq *
 ****************/

void ta_trim1(ta_opt_t *opt, char *seq)
{
	int i, j, k;
	kstring_t _str = {0,0,0}, *str = &_str;
	str->l = strlen(seq);
	str->s = malloc(str->l);
	for (i = 0; i < str->l; ++i)
		str->s[i] = seq_nt4_table[(uint8_t)seq[i]];
	for (j = 0; j < opt->n_adaps; ++j) {
		kswr_t r;
		double diff;
		int type;
		ta_adap_t *p = &opt->adaps[j];
		r = ksw_align(p->len, p->seq, str->l, (uint8_t*)str->s, 5, opt->mat, opt->go, opt->ge, KSW_XBYTE|KSW_XSTART|(opt->min_len * opt->sa), 0);
		++r.te; ++r.qe; // change to 0-based
		k = r.qe - r.qb < r.te - r.tb? r.qe - r.qb : r.te - r.tb;
		diff = (double)(k * opt->sa - r.score) / opt->sb / k;
		//printf("%d:%.3f [%d,%d):%d <=> [%d,%d):%d\n", r.score, diff, r.qb, r.qe, p->len, r.tb, r.te, (int)str.l);
		if (r.qb <= r.tb && p->len - r.qe <= str->l - r.te) { // contained
			if (r.qb * opt->sa > opt->sa + opt->sb) continue;
			if ((p->len - r.qe) * opt->sa > opt->sa + opt->sb) continue;
			type = 1;
		} else if (r.qb <= r.tb) { // 3'-end overlap
			if (r.qb * opt->sa > opt->sa + opt->sb) continue;
			if ((str->l - r.te) * opt->sa > opt->sa + opt->sb) continue;
			type = 2;
		} else { // 5'-end overlap
			if ((p->len - r.qe) * opt->sa > opt->sa + opt->sb) continue;
			if (r.tb * opt->sa > opt->sa + opt->sb) continue;
			type = 3;
		}
		if (p->type == 5) {
			if (r.tb == 0 && r.qe == p->len && (r.te - r.tb) * opt->sa == r.score)
				type = 4;
		} else if (p->type == 3) {
			if (r.qb == 0 && r.te == str->l && (r.te - r.tb) * opt->sa == r.score)
				type = 4;
		}
		if (type == 4) { // exact match
			if (r.te - r.tb < opt->min_len) continue;
		} else { // inexact match
			if (r.score < opt->min_sc || diff > opt->max_diff) continue;
		}
		__sync_fetch_and_add(&p->cnt, 1);
		if (p->type == 5) {
			k = r.te + (p->len - r.qe);
			k = k < str->l? k : str->l;
			for (i = 0; i < k; ++i) seq[i] = 'X';
		} else if (p->type == 3) {
			k = r.tb > r.qb? r.tb - r.qb : 0;
			for (i = k; i < str->l; ++i) seq[i] = 'X';
		}
	}
	free(str->s);
}

static void apply_trim(int min_trim, int l_seq, char *seq, char *qual)
{
	int i;
	while (l_seq > min_trim && seq[l_seq-1] == 'X') --l_seq;
	for (i = 0; l_seq - i > min_trim && seq[i] == 'X';) ++i;
	if (i > 0) {
		memmove(seq, seq + i, l_seq - i + 1);
		if (qual) memmove(qual, qual + i, l_seq - i + 1);
		l_seq -= i;
	}
	seq[l_seq] = 0;
	if (qual) qual[l_seq] = 0;
}

/**********************
 * Callback functions *
 **********************/

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	ta_opt_t *opt;
} data_for_t;

static void worker_for(void *_data, long i, int tid)
{
	data_for_t *data = (data_for_t*)_data;
	ta_trim1(data->opt, data->seqs[i].seq);
}

static void *worker_pipeline(void *shared, int step, void *_data)
{
	int i;
	ta_opt_t *opt = (ta_opt_t*)shared;
	if (step == 0) {
		data_for_t *ret;
		ret = calloc(1, sizeof(data_for_t));
		ret->seqs = bseq_read(opt->ks, opt->chunk_size, &ret->n_seqs);
		ret->opt = opt;
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		data_for_t *data = (data_for_t*)_data;
		kt_for(opt->n_threads, worker_for, data, data->n_seqs);
		return data;
	} else if (step == 2) {
		data_for_t *data = (data_for_t*)_data;
		for (i = 0; i < data->n_seqs; ++i) {
			bseq1_t *s = &data->seqs[i];
			if (opt->min_trim_len < s->l_seq)
				apply_trim(opt->min_trim_len, s->l_seq, s->seq, s->qual);
			putchar(s->qual? '@' : '>'); puts(s->name);
			puts(s->seq);
			if (s->qual) {
				puts("+"); puts(s->qual);
			}
			free(s->seq); free(s->qual); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

/*****************
 * Main function *
 *****************/

int main(int argc, char *argv[])
{
	int c, i, j;
	ta_opt_t opt;

	ta_opt_init(&opt);
	while ((c = getopt(argc, argv, "5:3:s:p:l:t:v")) >= 0) {
		if (c == '5' || c == '3') ta_opt_add_adap(&opt, c - '0', optarg);
		else if (c == 's') opt.min_sc = atoi(optarg);
		else if (c == 'd') opt.max_diff = atof(optarg);
		else if (c == 'l') opt.min_len = atoi(optarg);
		else if (c == 'p') opt.n_threads = atoi(optarg);
		else if (c == 't') opt.min_trim_len = atoi(optarg);
		else if (c == 'v') {
			puts(VERSION);
			return 0;
		}
	}

	if (opt.n_adaps == 0) ta_opt_default_adaps(&opt);

	if (optind == argc && isatty(fileno(stdin))) {
		fprintf(stderr, "Usage: trimadap-mt [options] <in.fq>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -5 STR     5'-end adapter\n");
		fprintf(stderr, "  -3 STR     3'-end adapter\n");
		fprintf(stderr, "  -l INT     min length [%d]\n", opt.min_len);
		fprintf(stderr, "  -s INT     min score [%d]\n", opt.min_sc);
		fprintf(stderr, "  -t INT     trim down [don't trim]\n");
		fprintf(stderr, "  -d FLOAT   max difference [%.3f]\n", opt.max_diff);
		fprintf(stderr, "  -p INT     number of trimmer threads [%d]\n", opt.n_threads);
		fprintf(stderr, "  -v         print version number\n");
		return 1; // FIXME: memory leak
	}

	ta_opt_open(&opt, optind < argc? argv[optind] : 0);
	kt_pipeline(2, worker_pipeline, &opt, 3);
	for (j = 0; j < opt.n_adaps; ++j) {
		ta_adap_t *p = &opt.adaps[j];
		fprintf(stderr, "%-15ld ", (long)p->cnt);
		for (i = 0; i < p->len; ++i) fputc("ACGTN"[(int)p->seq[i]], stderr);
		fputc('\n', stderr);
	}
	ta_opt_free(&opt);
	return 0;
}
