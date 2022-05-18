#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "boxver.h"
#include "sam.h"
#include "kstring.h"
#include "ksort.h"
#include "faidx.h"

typedef struct {
	int is_bam, print_bp, min_len, min_sc, min_q, max_gap, is_vcf, min_tip_q;
	float mask_level;
} cmdopt_t;

typedef struct {
	int64_t L, n_un, l_un, n_dropped;
	int n_b[5], n_bg[5];
	int n, m;
	int *len;
} stat_t;

typedef struct {
	int tid, pos, len, qlen, rlen, flag, mapq, qbeg, clip[2], ins, del, nm, score, tip_q[2];
	double diff;
} aln_t;

#define aln_lt_b(a, b) ((a).tid < (b).tid || ((a).tid == (b).tid && (a).pos < (b).pos))
KSORT_INIT(s2b, aln_t, aln_lt_b)

#define aln_lt_c(a, b) ((a).qbeg < (b).qbeg)
KSORT_INIT(s2c, aln_t, aln_lt_c)

#define intr_lt(a, b) ((a) > (b))
KSORT_INIT(intr, int, intr_lt)

static void count_break(int c[5], int n_aa, aln_t *aa, const cmdopt_t *o)
{
	int i, b[5] = { n_aa, 0, 0, 0, 0 };
	for (i = 0; i < n_aa; ++i) {
		aln_t *p = &aa[i];
		if (p->mapq < o->min_q) continue;
		++b[1];
		if (p->qlen >= 100) {
			++b[2];
			if (p->qlen >= 200) {
				++b[3];
				if (p->qlen >= 500) ++b[4];
			}
		}
	}
	for (i = 0; i < 5; ++i)
		if (b[i]) c[i] += b[i] - 1;
}

static void print_break_points(int n_aa, aln_t *aa, const cmdopt_t *o, const bam_hdr_t *h, const char *name, faidx_t *fai)
{
	int i, n_high = 0;
	if (n_aa < 2) return;
	for (i = n_high = 0; i < n_aa; ++i)
		if (aa[i].mapq > o->min_q) aa[n_high++] = aa[i];
	n_aa = n_high;
	if (n_aa < 2) return;
	ks_introsort(s2c, n_aa, aa);
	if (!o->is_vcf) {
		printf(">%s\n", name);
		for (i = 0; i < n_aa; ++i) { // print evidence
			aln_t *p = &aa[i];
			printf("#\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%.4f\t%d\t%d\t%d\n", p->qbeg, p->qbeg + p->qlen, "+-"[p->flag>>4&1],
					h->target_name[p->tid], p->pos, p->pos + p->rlen, p->mapq, p->diff, p->score, p->tip_q[0], p->tip_q[1]);
		}
	}
	for (i = 1; i < n_aa; ++i) {
		aln_t *q = &aa[i-1], *p = &aa[i];
		int min_mapq = q->mapq < p->mapq? q->mapq : p->mapq;
		int min_sc = q->score < p->score? q->score : p->score;
		int min_tip_q = q->tip_q[1] < p->tip_q[0]? q->tip_q[1] : p->tip_q[0];
		int qgap = p->qbeg - (q->qbeg + q->qlen);
		int qt_end, pt_start, rlen;
		char type;
		qt_end = (q->flag&16)? q->pos : q->pos + q->rlen;
		pt_start = (p->flag&16)? p->pos + p->rlen : p->pos;
		if (q->tid == p->tid) { // same chr
			if ((q->flag&16) == (p->flag&16)) { // same strand
				aln_t *r;
				int tmp;
				if (p->flag&16) r = p, p = q, q = r, tmp = qt_end, qt_end = pt_start, pt_start = tmp;
				type = q->pos >= p->pos? 'C' : pt_start - qt_end > qgap? 'D' : 'I';
			} else type = 'V';
		} else type = 'X';
		if (!o->is_vcf) {
			if (type == 'D' || type == 'I') {
				printf("%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", type, h->target_name[q->tid], qt_end, pt_start,
						qgap - (pt_start - qt_end), qgap, min_mapq, min_sc, min_tip_q);
			} else {
				printf("%c\t%s\t%d\t%c\t%s\t%d\t%c\t%d\t%d\t%d\n", type, 
						h->target_name[q->tid], qt_end, "+-"[!!(q->flag&16)],
						h->target_name[p->tid], pt_start, "+-"[!!(p->flag&16)],
						qgap, min_mapq, min_sc);
			}
		} else if (type == 'D' || type == 'I') {
			int max_start, max_end, len = qgap - (pt_start - qt_end);
			const char *type_str = 0;
			if (type == 'D') {
				len = -len;
				type_str = "DEL";
				max_start = qt_end < pt_start - len? qt_end : pt_start - len;
				max_end = pt_start > qt_end + len? pt_start : qt_end + len;
			} else {
				type_str = "INS";
				max_start = qt_end < pt_start? qt_end : pt_start;
				max_end = qt_end > pt_start? qt_end : pt_start;
			}
			char *ref = fai ? faidx_fetch_seq(fai, h->target_name[q->tid], max_start, max_start, &rlen) : "N";
			printf("%s\t%d\t.\t%c\t<%s>\t30\t%s\tSVTYPE=%s;END=%d;SVLEN=%d;QGAP=%d;MINMAPQ=%d;MINSC=%d;MINTIPQ=%d\n", h->target_name[q->tid],
					max_start+1, toupper(ref[0]), type_str, min_tip_q < o->min_tip_q? "LowSupp" : ".", type_str, max_end, len, qgap, min_mapq, min_sc, min_tip_q);
		} else {
			int dir = (p->flag&16)? '[' : ']';
			char *ref = fai ? faidx_fetch_seq(fai, h->target_name[q->tid], qt_end, qt_end, &rlen) : "N";
			char *join = fai ? faidx_fetch_seq(fai, h->target_name[p->tid], pt_start, pt_start, &rlen) : "N";
			printf("%s\t%d\t.\t%c\t", h->target_name[q->tid], qt_end + 1, toupper(ref[0]));
			if (q->flag&16) printf("%c%s:%d%c%c", dir, h->target_name[p->tid], pt_start + 1, dir, toupper(join[0]));
			else printf("%c%c%s:%d%c", toupper(join[0]), dir, h->target_name[p->tid], pt_start + 1, dir);
			printf("\t30\t%s\tSVTYPE=COMPLEX;QGAP=%d;MINMAPQ=%d;MINSC=%d;MINTIPQ=%d\n",
					min_tip_q < o->min_tip_q? "LowSupp" : ".", qgap, min_mapq, min_sc, min_tip_q);
		}
	}
}

static void analyze_aln(int n_aa, aln_t *aa, stat_t *s, const cmdopt_t *o, const bam_hdr_t *h, const char *name, faidx_t *fai)
{
	int n_tmp, i;
	aln_t *p, *tmp = 0;
	// if asked to print break points only
	if (o->print_bp) {
		if (n_aa > 1) print_break_points(n_aa, aa, o, h, name, fai);
		return;
	}
	// special treatment of unmapped
	if (n_aa == 1 && (aa[0].flag&4)) {
		++s->n_un; s->l_un += aa[0].len;
		return;
	}
	tmp = (aln_t*)alloca(n_aa * sizeof(aln_t));
	// apply mask_level
	if (n_aa > 1) { // multi-part alignment; check if this is a BWA-SW problem
		aln_t *p, *q;
		for (n_tmp = 0, p = aa; p < aa + n_aa; ++p) {
			int dropped = 0;
			for (q = tmp; q < tmp + n_tmp; ++q) {
				int beg = p->qbeg > q->qbeg? p->qbeg : q->qbeg;
				int end = p->qbeg + p->qlen < q->qbeg + q->qlen? p->qbeg + p->qlen : q->qbeg + q->qlen;
				if (beg < end && (double)(end - beg) > p->qlen * o->mask_level) {
					dropped = 1;
					break;
				}
			}
			if (!dropped) tmp[n_tmp++] = *p;
			else ++s->n_dropped;
		}
		memcpy(aa, tmp, n_tmp * sizeof(aln_t));
		n_aa = n_tmp;
		count_break(s->n_b, n_aa, aa, o);
	}
	if (s->n + n_aa > s->m) {
		s->m = s->n + n_aa;
		kroundup32(s->m);
		s->len = (int*)realloc(s->len, sizeof(int) * s->m);
	}
	for (p = aa; p < aa + n_aa; ++p) s->len[s->n++] = p->qlen;
	if (n_aa) { // still multi-part
		ks_introsort(s2b, n_aa, aa);
		for (i = 1; i < n_aa; ++i) {
			aln_t *p = &aa[i], *q = &aa[i-1];
			if (p->tid == q->tid && (p->flag&16) == (q->flag&16)) {
				int gapr = p->pos - (q->pos + q->rlen);
				int gapq = p->clip[0] - (q->clip[0] + q->qlen);
				if (gapr < 0) gapr = -gapr;
				if (gapq < 0) gapq = -gapq;
				if (gapr < o->max_gap && gapq < o->max_gap) {
					p->qlen = p->clip[0] + p->qlen - q->clip[0]; p->clip[0] = q->clip[0];
					p->rlen = p->pos + p->rlen - q->pos; p->pos = q->pos;
					q->flag |= 4;
				}
			}
		}
		for (n_tmp = i = 0; i < n_aa; ++i)
			if ((aa[i].flag&4) == 0) tmp[n_tmp++] = aa[i];
		memcpy(aa, tmp, n_tmp * sizeof(aln_t));
		n_aa = n_tmp;
		count_break(s->n_bg, n_aa, aa, o);
	}
}

int main_abreak(int argc, char *argv[])
{
	cmdopt_t o;
	samFile *in;
	bam_hdr_t *h;
	int c, k, n_aa = 0, m_aa = 0;
	aln_t a, *aa = 0;
	kstring_t last, out;
	stat_t s;
	bam1_t *b;
	faidx_t *fai = NULL;
	char *fname = 0;
	
	memset(&a, 0, sizeof(aln_t));
	memset(&s, 0, sizeof(stat_t));
	memset(&last, 0, sizeof(kstring_t));
	memset(&out, 0, sizeof(kstring_t));
	memset(&o, 0, sizeof(cmdopt_t));
	o.min_len = 150; o.min_q = 10; o.mask_level = 0.5; o.max_gap = 500; o.min_tip_q = 10;
	while ((c = getopt(argc, argv, "ul:bq:m:g:pcs:d:f:")) >= 0)
		if (c == 'b') o.is_bam = 1;
		else if (c == 'u') o.print_bp = 1, o.min_sc = 80, o.min_q = 40;
		else if (c == 'p') o.print_bp = 1;
		else if (c == 's') o.min_sc = atoi(optarg);
		else if (c == 'l') o.min_len = atoi(optarg);
		else if (c == 'q') o.min_q = atoi(optarg);
		else if (c == 'm') o.mask_level = atof(optarg);
		else if (c == 'g') o.max_gap = atoi(optarg);
		else if (c == 'c') o.is_vcf = o.print_bp = 1;
		else if (c == 'd') o.min_tip_q = atoi(optarg);
		else if (c == 'f') { fname = optarg; fai = fai_load(fname); }
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   htscmd abreak [options] <unsrt.sam>|<unsrt.bam>\n\n");
		fprintf(stderr, "Options: -b        assume the input is BAM (default is SAM)\n");
		fprintf(stderr, "         -l INT    exclude contigs shorter than INT [%d]\n", o.min_len);
		fprintf(stderr, "         -s INT    exclude alignemnts with score less than INT [%d]\n", o.min_sc);
		fprintf(stderr, "         -q INT    exclude alignments with mapQ below INT [%d]\n", o.min_q);
		fprintf(stderr, "         -d INT    filter calls with min flanking depth below INT in VCF [%d]\n", o.min_tip_q);
		fprintf(stderr, "         -p        print break points\n");
		fprintf(stderr, "         -c        VCF output (force -p)\n");
		fprintf(stderr, "         -u        unitig SV calling mode (-pq40 -s80)\n");
		fprintf(stderr, "         -f FILE   faidx indexed reference fasta (for -u)\n\n");
		fprintf(stderr, "         -m FLOAT  exclude aln overlapping another long aln by FLOAT fraction (effective w/o -p) [%g]\n", o.mask_level);
		fprintf(stderr, "         -g INT    join alignments separated by a gap shorter than INT bp (effective w/o -p) [%d]\n\n", o.max_gap);
		fprintf(stderr, "Note: recommended BWA-MEM setting is '-x intractg'. In the default output:\n\n");
		fprintf(stderr, "        >qName\n");
		fprintf(stderr, "        #      qStart  qEnd   strand   tName     tStart   tEnd     mapQ     perBaseDiv  alnScore\n");
		fprintf(stderr, "        [DI]   tName   tEnd1  tEnd2    inDelLen  qGapLen  minMapQ  minSc\n");
		fprintf(stderr, "        [CXV]  tName1  tEnd1  strand1  tName2    tEnd2    strand2  qGapLen  minMapQ     minSc\n\n");
		return 1;
	}

	in = sam_open(argv[optind], o.is_bam? "rb" : "r", 0);
	h = sam_hdr_read(in);
	b = bam_init1();

	if (o.is_vcf) {
		printf("##fileformat=VCFv4.1\n");
		printf("##source=htsbox-abreak-%s\n", HTSBOX_VERSION);
		if (fai) {
			printf("##reference=%s\n", fname);
			int i, n = faidx_fetch_nseq(fai);
			for (i=0; i<n; i++) {
				const char *seq = faidx_iseq(fai,i);
				int len = faidx_seq_len(fai, seq);
				printf("##contig=<ID=%s,length=%d>\n", seq, len);
			}
		}
		printf("##ALT=<ID=DEL,Description=\"Deletion\">\n");
		printf("##ALT=<ID=INS,Description=\"Insertion\">\n");
		printf("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">\n");
		printf("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
		printf("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">\n");
		printf("##INFO=<ID=QGAP,Number=1,Type=Integer,Description=\"Length of gap on the query sequence\">\n");
		printf("##INFO=<ID=MINMAPQ,Number=1,Type=Integer,Description=\"Min flanking mapping quality\">\n");
		printf("##INFO=<ID=MINSC,Number=1,Type=Integer,Description=\"Min flanking alignment score\">\n");
		printf("##INFO=<ID=MINTIPQ,Number=1,Type=Integer,Description=\"Min quality/depth flanking the break point\">\n");
		printf("##FILTER=<ID=LowSupp,Description=\"MINTIPQ < %d\">\n", o.min_tip_q);
		printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
	}

	while (sam_read1(in, h, b) >= 0) {
		const uint32_t *cigar = bam_get_cigar(b);
		const uint8_t *qual = bam_get_qual(b);
		uint8_t *tmp = 0;
		int clip_q[2], is_rev = !!(b->core.flag&BAM_FREVERSE);
		if (last.s == 0 || strcmp(last.s, bam_get_qname(b))) {
			if (last.s) analyze_aln(n_aa, aa, &s, &o, h, last.s, fai);
			last.l = 0;
			kputs(bam_get_qname(b), &last);
			n_aa = 0;
		}
		a.qlen = a.rlen = a.ins = a.del = a.clip[0] = a.clip[1] = 0;
		clip_q[0] = qual[0], clip_q[1] = qual[b->core.l_qseq - 1];
		for (k = 0; k < b->core.n_cigar; ++k) {
			int op = bam_cigar_op(cigar[k]);
			int oplen = bam_cigar_oplen(cigar[k]);
			if (oplen == 0) continue;
			if ((bam_cigar_type(op)&1) && op != BAM_CSOFT_CLIP) a.qlen += oplen;
			if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
				a.clip[!!k] = oplen;
				if (op == BAM_CSOFT_CLIP)
					clip_q[!!k] = !k? qual[oplen] : qual[a.qlen + a.clip[0] - 1];
			}
			if (bam_cigar_type(op)&2) a.rlen += oplen;
			if (op == BAM_CINS) ++a.ins;
			else if (op == BAM_CDEL) ++a.del;
		}
		a.len = a.qlen + a.clip[0] + a.clip[1];
		if (a.len == 0) a.len = b->core.l_qseq;
		a.tid = b->core.tid; a.pos = b->core.pos; a.flag = b->core.flag; a.mapq = b->core.qual;
		a.qbeg = a.clip[is_rev];
		a.tip_q[0] = clip_q[is_rev], a.tip_q[1] = clip_q[!is_rev];
		if ((tmp = bam_aux_get(b, "NM")) != 0) {
			a.nm = bam_aux2i(tmp);
			a.diff = (double)a.nm / (a.qlen + a.del);
		} else a.nm = -1, a.diff = -1.;
		a.score = (tmp = bam_aux_get(b, "AS")) != 0? bam_aux2i(tmp) : a.qlen - a.ins; // if "AS" absent, use the total "M" length as score
		if (a.len >= o.min_len && a.score >= o.min_sc) {
			if (n_aa == m_aa) {
				m_aa = m_aa? m_aa<<1 : 8;
				aa = (aln_t*)realloc(aa, sizeof(aln_t) * m_aa);
			}
			aa[n_aa++] = a;
		}
	}
	analyze_aln(n_aa, aa, &s, &o, h, last.s, fai);
	bam_hdr_destroy(h);
	sam_close(in);
	if (fai) fai_destroy(fai);
	if (!o.print_bp) {
		uint64_t L = 0, tmp;
		int N50;
		ks_introsort(intr, s.n, s.len);
		for (k = 0; k < s.n; ++k) L += s.len[k];
		for (k = 0, tmp = 0; k < s.n; ++k)
			if ((tmp += s.len[k]) >= L/2) break;
		N50 = s.len[k];
		printf("Number of unmapped contigs: %ld\n", (long)s.n_un);
		printf("Total length of unmapped contigs: %ld\n", (long)s.l_un);
		//printf("Number of alignments dropped due to excessive overlaps: %ld\n", (long)s.n_dropped); // this is caused by a bwa-sw only bug
		printf("Mapped contig bases: %ld\n", (long)L);
		printf("Mapped N50: %d\n", N50);
		printf("Number of break points: %d\n", s.n_b[0]);
		printf("Number of Q%d break points longer than (0,100,200,500)bp: (%d,%d,%d,%d)\n", o.min_q, s.n_b[1], s.n_b[2], s.n_b[3], s.n_b[4]);
		printf("Number of break points after patching gaps short than %dbp: %d\n", o.max_gap, s.n_bg[0]);
		printf("Number of Q%d break points longer than (0,100,200,500)bp after gap patching: (%d,%d,%d,%d)\n", o.min_q, s.n_bg[1], s.n_bg[2], s.n_bg[3], s.n_bg[4]);
	}
	return 0;
}
