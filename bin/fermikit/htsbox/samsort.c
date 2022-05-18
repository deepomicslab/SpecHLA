#include "sam.h"
#include "ksort.h"

typedef struct {
	uint64_t pos, idx;
	bam1_t *b;
} aln_t;

#define aln_lt(a, b) ((a).pos < (b).pos || ((a).pos == (b).pos && (a).idx < (b).idx))
KSORT_INIT(aln, aln_t, aln_lt)

int main_samsort(int argc, char *argv[])
{
	samFile *in;
	int c, clevel = -1, sam_in = 0, n_threads = 2, ignore_sam_err = 0;
	char moder[8], modew[8];
	bam_hdr_t *h;
	uint64_t k, n_aln = 0, m_aln = 0;
	bam1_t *b;
	aln_t *aln = 0;
	BGZF *fp = 0;

	while ((c = getopt(argc, argv, "SIl:t:")) >= 0) {
		if (c == 'S') sam_in = 1;
		else if (c == 'I') ignore_sam_err = 1;
		else if (c == 'l') clevel = atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: samsort [-S] [-l clevel] [-t nThreads] <in.bam>|<in.sam>\n");
		return 1;
	}

	strcpy(moder, "r");
	if (!sam_in) strcat(moder, "b");
	in = sam_open(argv[optind], moder, 0);

	h = sam_hdr_read(in);
	h->ignore_sam_err = ignore_sam_err;
	b = bam_init1();
	while (sam_read1(in, h, b) >= 0) {
		aln_t *a;
		if (n_aln == m_aln) {
			m_aln = m_aln? m_aln<<1 : 64;
			aln = (aln_t*)realloc(aln, m_aln * sizeof(aln_t));
		}
		a = &aln[n_aln++];
		a->idx = n_aln - 1;
		a->pos = (uint64_t)b->core.tid << 32 | b->core.pos;
		a->b = bam_init1();
		bam_copy1(a->b, b);
	}
	bam_destroy1(b);
	sam_close(in);

	ks_introsort(aln, n_aln, aln);

	strcpy(modew, "w");
	if (clevel >= 0 && clevel <= 9) sprintf(modew+1, "%d", clevel);
	fp = bgzf_dopen(fileno(stdout), modew);
	bgzf_mt(fp, n_threads, 256);

	bam_hdr_write(fp, h);
	for (k = 0; k < n_aln; ++k) {
		bam_write1(fp, aln[k].b);
		bam_destroy1(aln[k].b);
	}
	bgzf_close(fp);
	free(aln);

	bam_hdr_destroy(h);
	return 0;
}
