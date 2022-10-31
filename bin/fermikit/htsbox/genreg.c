#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include "sam.h"

int main_genreg(int argc, char *argv[])
{
	BGZF *fp;
	bam_hdr_t *h;
	bam1_t *b;
	double size = 1e9, gap = 10, x;
	char *p;
	int64_t s = 0;
	int c, start = -1, end = -1, tid = -1;

	while ((c = getopt(argc, argv, "n:g:")) >= 0) {
		if (c == 'n' || c == 'g') {
			x = strtod(optarg, &p);
			if (*p == 'G' || *p == 'g') x *= 1e9;
			else if (*p == 'M' || *p == 'm') x *= 1e6;
			else if (*p == 'K' || *p == 'k') x *= 1e3;
			if (c == 'n') size = x;
			else gap = x;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: genreg [-n nBases=1g] [-g gapSize=%.0lf] <sorted.bam>\n", gap);
		return 1;
	}

	fp = strcmp(argv[optind], "-")? bgzf_open(argv[optind], "r") : bgzf_dopen(fileno(stdin), "r");
	h = bam_hdr_read(fp);
	b = bam_init1();

	while (bam_read1(fp, b) >= 0) {
		int st, en;
		if (b->core.tid < 0) break;
		if ((b->core.flag&BAM_FSECONDARY)) continue;
		st = b->core.pos;
		en = st + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
		if (st == en) en = st + 1;
		if (en > h->target_len[b->core.tid])
			en = h->target_len[b->core.tid];
		if (b->core.tid != tid || (s + (en - st) > size && st - end > gap)) {
			if (tid >= 0) printf("%s\t%d\t%d\n", h->target_name[tid], start, end);
			tid = b->core.tid, start = st, end = en, s = en - st;
		} else {
			end = end > en? end : en;
			s += en - st;
		}
	}
	if (tid >= 0) printf("%s\t%d\t%d\n", h->target_name[tid], start, end);

	bam_destroy1(b);
	bam_hdr_destroy(h);
	bgzf_close(fp);
	return 0;
}
