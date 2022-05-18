#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

/****************************
 *** Hard coded constants ***
 ****************************/

#define FMC_NOHIT_PEN 63
#define FMC_MAX_PATHS 8
#define FMC_Q_MAX_OUT 33
#define FMC_Q_MAX     41
#define FMC_Q_1       25 // IMPORTANT: FMC_Q_1*2 > FMC_Q_MAX
#define FMC_Q_NULL    31
#define FMC_SI_GAP    3

/******************
 *** Parameters ***
 ******************/

int fmc_verbose = 3;

#define FMC_KMER_MAGIC "FCK\2" // IMPORTANT: change this magic whenever fmc_collect_opt_t is changed!!!

typedef struct {
	int k:16, suf_len:16;
	int min_occ:16, max_ec_depth:16;
	int q1_depth, dummy;
	double a1, a2, err, prior;
} fmc_collect_opt_t;

typedef struct {
	fmc_collect_opt_t c;
	int n_threads, ecQ, defQ;
	int max_heap_size;
	int max_penalty_diff;
	int show_ori_name;
	int max_dist4;
	int drop_reads;
	int64_t batch_size;
} fmc_opt_t;

void fmc_opt_init(fmc_opt_t *opt)
{
	memset(opt, 0, sizeof(fmc_opt_t));
	opt->c.k = 17;
	opt->c.suf_len = 1;
	opt->c.min_occ = 3;
	opt->c.q1_depth = 17; // if there q1_depth bases but only one 2nd-best, correct regardless of the quality
	opt->c.max_ec_depth = opt->c.min_occ - 1; // if there are more than max_ec_depth 2nd-best bases, don't correct
	opt->c.a1 = 0.05;
	opt->c.a2 = 10;
	opt->c.err = 0.005;
	opt->c.prior = 0.99;

	opt->n_threads = 1;
	opt->defQ = 17;
	opt->ecQ = 30;
	opt->max_heap_size = 256;
	opt->max_penalty_diff = 60;
	opt->batch_size = (1ULL<<28) - (1ULL<<20);
	opt->max_dist4 = 8;
}

void kt_for(int n_threads, void (*func)(void*,long,int), void *shared, long n_items);
double cputime(void);
double realtime(void);
void liftrlimit(void);

/****************************
 *** Consensus generation ***
 ****************************/

#include <math.h>

double kf_lgamma(double z) // log(gamma(z))
{
	double x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

double fmc_beta_binomial(int n, int k, double a, double b) // beta-binomial density
{
	double x, y, z;
	x = lgamma(n + 1) - (lgamma(k + 1) + lgamma(n - k + 1));
	y = lgamma(k + a) + lgamma(n - k + b) - lgamma(n + a + b);
	z = lgamma(a + b) - (lgamma(a) + lgamma(b));
	return exp(x + y + z);
}

uint8_t *fmc_precal_qtab(int max, double e1, double e2, double a1, double a2, double prior1)
{ // precalculate a lookup table q(n,k), the consensus quality given n bases with k differences
	int n, k;
	uint8_t *qtab;
	double b1 = a1 * (1 - e1) / e1, b2 = a2 * (1 - e2) / e2;

	assert(FMC_Q_MAX>>1 < FMC_Q_1);
	qtab = calloc(max * max, 1);
	for (n = 1; n < max; ++n) {
		uint8_t *qn = &qtab[n*max];
		//fprintf(stderr, "=> %d <=", n);
		for (k = 0; k < n; ++k) {
			double p1, p2;
			int q;
			p1 = fmc_beta_binomial(n, k, a1, b1);
			p2 = fmc_beta_binomial(n, k, a2, b2);
			q = (int)(-4.343 * log(1. - p1 * prior1 / (p1 * prior1 + p2 * (1-prior1))) + .499);
			qn[k] = (q < FMC_Q_MAX? q : FMC_Q_MAX) >> 1;
			//if (q) fprintf(stderr, "\t%d:%d", k, q);
		}
		//fprintf(stderr, "\n");
	}
	return qtab;
}

/******************
 *** Hash table ***
 ******************/

#include "kvec.h"
#include "khash.h"

#define fmc_cell_get_key(x) ((x)>>28)
#define fmc_cell_get_val(x, is_right) ((x)>>((is_right)?14:0)&0x3fff)
#define fmc_cell_set_keyval(key, val0, val1) ((x)<<28|(val1)<<14|(val0))
#define fmc_cell_set_val(b1, b2, q1, q2) ((b1)<<12|(b2)<<10|(q1)<<5|(q2))

#define fmc_cell_get_b1(v) ((v)>>12&0x3)
#define fmc_cell_get_b2(v) ((v)>>10&0x3)
#define fmc_cell_get_q1(v) (((v)>>5&0x1f)<<1)
#define fmc_cell_get_q2(v) (((v)&0x1f)<<1)

#define fmc_cell_has_b1(v) (((v)>>5&0x1f) != FMC_Q_NULL)
#define fmc_cell_has_b2(v) (((v)&0x1f) != FMC_Q_NULL)

static inline uint64_t hash_64(uint64_t key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}

#define fmc_hash_func(a) hash_64(fmc_cell_get_key(a))
#define fmc_eq_func(a, b) (fmc_cell_get_key(a) == fmc_cell_get_key(b))
KHASH_INIT(fmc, uint64_t, char, 0, fmc_hash_func, fmc_eq_func)

typedef khash_t(fmc) fmc_hash_t;
typedef kvec_t(uint64_t) fmc64_v;

typedef struct {
	uint64_t suf:63, missing:1;
	uint64_t x;
} longcell_t;

#define fmc_longcell_get_key(a) ((a).suf<<36 | fmc_cell_get_key((a).x))
#define fmc_hash_func_long(a) hash_64(fmc_longcell_get_key(a))
#define fmc_eq_func_long(a, b) (fmc_longcell_get_key(a) == fmc_longcell_get_key(b))
KHASH_INIT(kache, longcell_t, char, 0, fmc_hash_func_long, fmc_eq_func_long)

typedef khash_t(kache) kmercache_t;

/*********************************
 *** Collect k-mer information ***
 *********************************/

#include "rld0.h"

typedef kvec_t(rldintv_t) rldintv_v;

rldintv_t *fmc_traverse(const rld_t *e, int depth) // traverse FM-index up to $depth
{
	rldintv_v stack = {0,0,0};
	rldintv_t *p, *ret;
	uint64_t x;

	ret = calloc(1<<depth*2, sizeof(rldintv_t));
	kv_pushp(rldintv_t, stack, &p);
	p->x[0] = p->x[1] = 0, p->x[2] = e->mcnt[0], p->info = 0;
	x = 0;
	while (stack.n) {
		rldintv_t top = kv_pop(stack);
		if (top.info>>2 > 0) {
			int shift = ((top.info>>2) - 1) << 1;
			x = (x & ~(3ULL<<shift)) | (uint64_t)(top.info&3)<<shift;
		}
		if (top.info>>2 != depth) {
			int c;
			rldintv_t t[6];
			rld_extend(e, &top, t, 1);
			for (c = 1; c < 5; ++c) {
				if (t[c].x[2] == 0) continue;
				t[c].info = ((top.info>>2) + 1) << 2 | (c - 1);
				kv_push(rldintv_t, stack, t[c]);
			}
		} else ret[x] = top, ret[x].info = x;
	}
	free(stack.a);
	return ret;
}

int fmc_intv2tip(uint8_t *qtab[2], const rldintv_t t[6], int max_ec_depth, int q1_depth, int min_occ) // given a "tip" compute the consensus and the quality
{
	int c, max_c, max_c2, q1, q2;
	uint64_t max, max2, rest, rest2, sum;
	for (c = 1, max = max2 = 0, max_c = max_c2 = 1, sum = 0; c <= 4; ++c) { // get the top 2 bases
		if (t[c].x[2] > max) max2 = max, max_c2 = max_c, max = t[c].x[2], max_c = c;
		else if (t[c].x[2] > max2) max2 = t[c].x[2], max_c2 = c;
		sum += t[c].x[2];
	}
	rest = sum - max; rest2 = sum - max - max2; // the following always stands: rest2 <= rest <= sum and 0 <= max2 <= max <= sum
	if (sum == 0) {
		q1 = q2 = FMC_Q_NULL;
	} else if (rest <= 1 && sum >= q1_depth) { // if rest=1, max=sum-1 and max2=1, then rest2=0
		q1 = FMC_Q_1;
		q2 = rest? FMC_Q_1 : FMC_Q_NULL;
	} else if (max < min_occ || rest > max_ec_depth) {
		q1 = q2 = 0;
	} else {
		if (sum > 255) {
			rest  = (int)(255. * rest  / sum + .499);
			rest2 = (int)(255. * rest2 / sum + .499);
			sum = 255;
		}
		q1 = qtab[0][sum<<8|rest];
		q2 = rest? qtab[1][sum<<8|rest2] : FMC_Q_NULL;
	}
	return fmc_cell_set_val(4-max_c, 4-max_c2, q1, q2);
}

void fmc_collect1(const rld_t *e, uint8_t *qtab[2], int suf_len, int depth, int min_occ, int max_ec_depth, int q1_depth, const rldintv_t *start, fmc64_v *a)
{
	rldintv_v stack = {0,0,0};
	uint64_t x = 0, *p;

	kv_push(rldintv_t, stack, *start);
	stack.a[0].info = 0;
	a->n = 0;
	while (stack.n) {
		rldintv_t top = kv_pop(stack);
		if (top.info>>2 > 0) {
			int shift = ((top.info>>2) - 1) << 1;
			x = (x & ~(3ULL<<shift)) | (uint64_t)(top.info&3)<<shift;
		}
		if (top.info>>2 == depth) { // reach the length; collect info at the two tips
			rldintv_t t[6];
			int val[2];
			kv_pushp(uint64_t, *a, &p);
			rld_extend(e, &top, t, 1); // backward tip
			val[0] = fmc_intv2tip(qtab, t, max_ec_depth, q1_depth, min_occ);
			rld_extend(e, &top, t, 0); // forward tip
			val[1] = fmc_intv2tip(qtab, t, max_ec_depth, q1_depth, min_occ);
			*p = fmc_cell_set_keyval(x, val[0], val[1]);
		} else {
			int c, end = (suf_len + (top.info>>2)) == (suf_len + depth) / 2? 2 : 4;
			rldintv_t t[6];
			rld_extend(e, &top, t, 1);
			for (c = 1; c <= end; ++c) {
				if (t[c].x[2] < min_occ) continue;
				t[c].info = ((top.info>>2) + 1) << 2 | (c - 1);
				kv_push(rldintv_t, stack, t[c]);
			}
		}
	}
	free(stack.a);
}

typedef struct {
	const fmc_opt_t *opt;
	const rld_t *e;
	uint8_t *qtab[2];
	rldintv_t *suf;
	fmc64_v *kmer;
	int depth;
} for_collect_t;

static void collect_func(void *shared, long i, int tid)
{
	for_collect_t *s = (for_collect_t*)shared;
	fmc_collect1(s->e, s->qtab, s->opt->c.suf_len, s->depth, s->opt->c.min_occ, s->opt->c.max_ec_depth, s->opt->c.q1_depth, &s->suf[i], &s->kmer[i]);
}

void fmc_kmer_stat(int suf_len, const fmc64_v *a)
{
	int i, n_suf = 1<<suf_len*2;
	int64_t tot = 0;
	for (i = 0; i < n_suf; ++i)
		tot += a[i].n<<1;
	fprintf(stderr, "[M::%s] %ld k-mers\n", __func__, (long)tot);
}

fmc64_v *fmc_collect(fmc_opt_t *opt, const char *fn_fmi)
{
	rld_t *e;
	double tc, tr;
	int depth = opt->c.k - opt->c.suf_len, n_suf = 1 << opt->c.suf_len*2;
	for_collect_t f;

	assert(0 < depth && depth <= 18);

	fprintf(stderr, "[M::%s] reading the FMD-index... ", __func__);
	tc = cputime(); tr = realtime();
	e = rld_restore(fn_fmi);
	fprintf(stderr, " in %.3f sec (%.3f CPU sec)\n", realtime() - tr, cputime() - tc);

	fprintf(stderr, "[M::%s] collecting high occurrence k-mers... ", __func__);
	tc = cputime(); tr = realtime();
	f.e = e; f.opt = opt; f.depth = depth;
	f.suf = fmc_traverse(e, opt->c.suf_len);
	f.qtab[0] = fmc_precal_qtab(1<<8, opt->c.err, 0.5,      opt->c.a1, opt->c.a2, opt->c.prior);
	f.qtab[1] = fmc_precal_qtab(1<<8, opt->c.err, 0.333333, opt->c.a1, opt->c.a2, opt->c.prior);
	f.kmer = calloc(n_suf, sizeof(fmc64_v));
	kt_for(opt->n_threads, collect_func, &f, n_suf);
	rld_destroy(e);
	free(f.qtab[0]); free(f.qtab[1]); free(f.suf);
	fprintf(stderr, "in %.3f sec (%.3f CPU sec)\n", realtime() - tr, cputime() - tc);

	fmc_kmer_stat(opt->c.suf_len, f.kmer);
	return f.kmer;
}

/****************************
 *** Write/read kmer list ***
 ****************************/

void fmc_kmer_write(FILE *fp, const fmc_opt_t *opt, const fmc64_v *a)
{
	int i, n = 1<<opt->c.suf_len*2;
	fwrite(FMC_KMER_MAGIC, 1, 4, fp);
	fwrite(&opt->c, sizeof(fmc_collect_opt_t), 1, fp);
	for (i = 0; i < n; ++i) {
		fwrite(&a[i].n, sizeof(size_t), 1, fp);
		fwrite(a[i].a, 8, a[i].n, fp);
	}
}

fmc64_v *fmc_kmer_read(FILE *fp, fmc_opt_t *opt)
{
	int i, n;
	char magic[4];
	fmc64_v *a;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, FMC_KMER_MAGIC, 4) != 0) {
		fprintf(stderr, "[E::%s] invalid file magic\n", __func__);
		return 0;
	}
	fread(&opt->c, sizeof(fmc_collect_opt_t), 1, fp);
	n = 1<<opt->c.suf_len*2;
	a = malloc(sizeof(fmc64_v) * n);
	for (i = 0; i < n; ++i) {
		fread(&a[i].n, sizeof(size_t), 1, fp);
		a[i].m = a[i].n;
		a[i].a = malloc(8 * a[i].n);
		fread(a[i].a, 8, a[i].n, fp);
	}
	return a;
}

fmc_hash_t **fmc_kmer2hash(const fmc_opt_t *opt, fmc64_v *a)
{
	int i, n = 1 << opt->c.suf_len*2;
	fmc_hash_t **h;
	double t;
	t = cputime();
	fprintf(stderr, "[M::%s] constructing the hash table... ", __func__);
	h = calloc(n, sizeof(void*));
	for (i = 0; i < n; ++i) {
		size_t j;
		int absent;
		fmc64_v *ai = &a[i];
		h[i] = kh_init(fmc);
		kh_resize(fmc, h[i], (int)(ai->n / .7 + 1.));
		for (j = 0; j < ai->n; ++j)
			kh_put(fmc, h[i], ai->a[j], &absent);
		assert(kh_size(h[i]) == ai->n);
		free(ai->a);
	}
	fprintf(stderr, "in %.3f CPU sec\n", cputime() - t);
	free(a);
	return h;
}

/************************
 *** Sequence reading ***
 ************************/

#include <zlib.h>
#include "kseq.h"
KSEQ_DECLARE(gzFile)

extern unsigned char seq_nt6_table[128];

typedef struct {
	int n, m;
	char **s, **q, **name;
} fmc_batch_t;

fmc_batch_t *fmc_batch_read(kseq_t *ks, int64_t batch_size)
{
	int64_t l = 0;
	fmc_batch_t *b;
	b = calloc(1, sizeof(fmc_batch_t));
	while (l < batch_size && kseq_read(ks) >= 0) {
		if (b->n >= b->m) {
			b->m = b->n + 1;
			kroundup32(b->m);
			b->s = realloc(b->s, b->m * sizeof(void*));
			b->q = realloc(b->q, b->m * sizeof(void*));
			b->name = realloc(b->name, b->m * sizeof(void*));
		}
		b->s[b->n] = strdup(ks->seq.s);
		b->q[b->n] = ks->qual.l? strdup(ks->qual.s) : 0;
		b->name[b->n++] = strdup(ks->name.s);
		l += ks->seq.l;
	}
	if (l == 0) {
		free(b);
		return 0;
	} else return b;
}

void fmc_batch_destroy(fmc_batch_t *b)
{
	int i;
	for (i = 0; i < b->n; ++i) {
		free(b->name[i]); free(b->s[i]); free(b->q[i]);
	}
	free(b->name); free(b->s); free(b->q); free(b);
}

#define STATE_N 0
#define STATE_M 1
#define STATE_I 2
#define STATE_D 3

typedef struct { // NOTE: unaligned memory
	uint8_t b:4, state:4;
	uint8_t ob:4, f:4; // ob/oq: original base/quality
	uint8_t q, oq, min_diff, cov;
	int i;
} ecbase_t;

typedef kvec_t(ecbase_t) ecseq_t;

int fmc_seq_conv(const char *s, const char *q, int defQ, ecseq_t *seq)
{
	int i, l;
	l = strlen(s);
	kv_resize(ecbase_t, *seq, l);
	seq->n = l;
	for (i = 0; i < l; ++i) {
		ecbase_t *c = &seq->a[i];
		c->b = c->ob = seq_nt6_table[(int)s[i]] - 1;
		c->q = q? q[i] - 33 : defQ;
		c->q = c->oq = c->q < FMC_Q_MAX? c->q : FMC_Q_MAX;
		c->q = c->q > 5? c->q : 5;
		c->state = STATE_N;
		c->i = i;
		c->f = 0;
		c->min_diff = 0xff;
	}
	return l;
}

void fmc_seq_cpy_no_del(ecseq_t *dst, const ecseq_t *src)
{
	int i;
	kv_resize(ecbase_t, *dst, src->n);
	dst->n = 0;
	memcpy(dst->a, src->a, src->n * sizeof(ecbase_t));
	for (i = 0; i < src->n; ++i)
		if (src->a[i].state != STATE_D)
			dst->a[dst->n++] = src->a[i];
}

static inline ecbase_t ecbase_comp(const ecbase_t *b)
{
	ecbase_t r = *b;
	r.b = b->b < 4? 3 - b->b : 4;
	r.ob = b->ob < 4? 3 - b->ob : 4;
	return r;
}

void fmc_seq_revcomp(ecseq_t *seq)
{
	int i;
	for (i = 0; i < seq->n>>1; ++i) {
		ecbase_t tmp;
		tmp = ecbase_comp(&seq->a[i]);
		seq->a[i] = ecbase_comp(&seq->a[seq->n - 1 - i]);
		seq->a[seq->n - 1 - i] = tmp;
	}
	if (seq->n&1) seq->a[i] = ecbase_comp(&seq->a[i]);
}

/************************
 *** Error correction ***
 ************************/

#include "ksort.h"

typedef struct {
	uint64_t kmer[2];
	int penalty;
	int k; // pos on the stack
	int i; // pos of the next-to-add base on the sequence 
	uint32_t state:3, last_solid:29;
	uint64_t ec_pos4;
} echeap1_t;

#define echeap1_lt(a, b) ((a).penalty > (b).penalty)
KSORT_INIT(ec, echeap1_t, echeap1_lt)

typedef struct {
	int parent;
	int i, penalty;
	uint8_t state, base, dummy1, dummy2; // dummy for memory alignment
} ecstack1_t;

typedef kvec_t(echeap1_t)  echeap_t;
typedef kvec_t(ecstack1_t) ecstack_t;

typedef struct {
	ecseq_t ori, tmp[2], seq, ec_for;
	echeap_t heap;
	ecstack_t stack;
	kmercache_t *cache;
} fmc_aux_t;

fmc_aux_t *fmc_aux_init()
{
	fmc_aux_t *a;
	a = calloc(1, sizeof(fmc_aux_t));
	a->cache = kh_init(kache);
	return a;
}

void fmc_aux_destroy(fmc_aux_t *a)
{
	free(a->seq.a); free(a->ori.a); free(a->tmp[0].a); free(a->tmp[1].a);
	free(a->heap.a); free(a->stack.a);
	kh_destroy(kache, a->cache);
	free(a);
}

static inline void append_to_kmer(int k, uint64_t kmer[2], int a)
{
	kmer[0] = (kmer[0]<<2 & ((1ULL<<(k<<1)) - 1)) | a;
	kmer[1] = kmer[1]>>2 | (uint64_t)(3 - a) << ((k-1)<<1);
}

static inline int kmer_lookup(int k, int suf_len, uint64_t kmer[2], fmc_hash_t **h, kmercache_t *cache)
{
	int absent, i = (kmer[0]>>(k>>1<<1)&3) < 2? 0 : 1;
	khint_t kh, kc;
	longcell_t x, *p;

	x.suf = kmer[i] & ((1<<(suf_len<<1)) - 1);
	x.x = kmer[i] >> (suf_len<<1) << 28;
	x.missing = 0; // This line suppresses a gcc warning. It has no effect.
	kc = kh_put(kache, cache, x, &absent);
	p = &kh_key(cache, kc);
	if (absent) {
		fmc_hash_t *g = h[x.suf];
		kh = kh_get(fmc, g, x.x);
		if (kh == kh_end(g)) p->missing = 1;
		else p->missing = 0, p->x = kh_key(g, kh);
	}
	if (fmc_verbose >= 6) {
		int i, which = (kmer[0]>>(k>>1<<1)&3) < 2? 0 : 1;
		int val = p->missing? -1 : fmc_cell_get_val(p->x, !which);
		fprintf(stderr, "?? ");
		for (i = k-1; i >= 0; --i) fputc("ACGT"[kmer[0]>>2*i&3], stderr); fprintf(stderr, " - ");
		for (i = k-1; i >= 0; --i) fputc("ACGT"[kmer[1]>>2*i&3], stderr);
		fprintf(stderr, " - [%c] %lx", "+-"[which], (long)kmer[which]);
		if (p->missing) fprintf(stderr, " - NOHIT\n");
		else fprintf(stderr, " - %c%d\n", "ACGTN"[fmc_cell_get_b1(val)], fmc_cell_get_q1(val));
	}
	return p->missing? -1 : fmc_cell_get_val(p->x, !i);
}

static inline void update_aux(int k, fmc_aux_t *a, const echeap1_t *p, int b, int state, int penalty, int is_diff, int is_solid)
{
	ecstack1_t *q;
	echeap1_t *r;
	// update the stack
	kv_pushp(ecstack1_t, a->stack, &q);
	q->parent = p->k;
	q->i = p->i;
	q->base = b;
	q->state = state;
	q->penalty = p->penalty + penalty;
	// update the heap
	kv_pushp(echeap1_t, a->heap, &r);
	r->penalty = q->penalty;
	r->k = a->stack.n - 1;
	r->kmer[0] = p->kmer[0], r->kmer[1] = p->kmer[1];
	r->state = state;
	r->i = state == STATE_I? p->i : p->i + 1;
	r->last_solid = is_solid && is_diff? r->i : p->last_solid;
	r->ec_pos4 = is_diff? p->ec_pos4<<16 | (p->i + 1) : p->ec_pos4;
	if (fmc_verbose >= 6)fprintf(stderr, "+> [%d] (%d,%c), ipen=%d, state=%c, w=%d\n", r->k, r->i, "ACGTN"[b], penalty, "NMID"[state], q->penalty);
	if (state != STATE_D) append_to_kmer(k, r->kmer, b);
	ks_heapup_ec(a->heap.n, a->heap.a);
}

static void path_backtrack(const ecstack_t *a, int start, const ecseq_t *o, ecseq_t *s)
{
	int i = start, last = -1;
	s->n = 0;
	while (i >= 0) {
		ecbase_t *c;
		ecstack1_t *p = &a->a[i];
		int is_match = (p->state == STATE_M || p->state == STATE_N);
		kv_pushp(ecbase_t, *s, &c);
		c->state = p->state;
		c->i = o->a[p->i].i;
		c->b = p->state == STATE_D? 4 : p->base;
		c->q = is_match && o->a[p->i].b == p->base? o->a[p->i].q : 0; // if we make a correction, reduce the base quality to 0.
		c->ob = p->state == STATE_I? 4 : o->a[p->i].ob; // no ob/oq for inserted bases
		c->oq = p->state == STATE_I? 0 : o->a[p->i].oq;
		c->f = (p->state != STATE_N);
		c->min_diff = is_match? o->a[p->i].min_diff : 0xff;
		last = p->i;
		i = p->parent;
	}
	assert(last > 0);
	for (i = last - 1; i >= 0; --i)
		kv_push(ecbase_t, *s, o->a[i]);
	for (i = 0; i < s->n>>1; ++i) { // reverse
		ecbase_t tmp = s->a[i];
		s->a[i] = s->a[s->n - 1 - i];
		s->a[s->n - 1 - i] = tmp;
	}
}

static void path_adjustq(int diff, ecseq_t *s1, const ecseq_t *s2)
{
	int i1 = 0, i2 = 0;
	while (i1 < s1->n && i2 < s2->n) {
		ecbase_t *b1;
		const ecbase_t *b2;
		b1 = &s1->a[i1];
		b2 = &s2->a[i2];
		if (b1->b != b2->b)
			b1->min_diff = b1->min_diff < diff? b1->min_diff : diff;
		if (b1->state == STATE_I && b2->state != STATE_I) ++i1;
		else if (b2->state == STATE_I && b1->state != STATE_I) ++i2;
		else ++i1, ++i2;
	}
	for (; i1 < s1->n; ++i1) {
		ecbase_t *b = &s1->a[i1];
		b->min_diff = b->min_diff < diff? b->min_diff : diff;
	}
}

static int kmer_cov(const fmc_opt_t *opt, ecseq_t *seq, fmc_hash_t **h, kmercache_t *cache)
{
	int i, j, l, v, n_si, in_si;
	uint64_t kmer[2], kfirst[2];
	// compute the k-mer coverage for k-mers contained in the read
	kmer[0] = kmer[1] = kfirst[0] = kfirst[1] = 0;
	for (i = 0; i < seq->n; ++i) seq->a[i].cov = 0;
	for (i = l = 0; i < seq->n; ++i) {
		ecbase_t *p = &seq->a[i];
		if (p->b > 3) l = 0, kmer[0] = kmer[1] = 0;
		else ++l, append_to_kmer(opt->c.k, kmer, p->b);
		if (l >= opt->c.k && (v = kmer_lookup(opt->c.k, opt->c.suf_len, kmer, h, cache)) >= 0)
			for (j = 0; j < opt->c.k; ++j)
				++seq->a[i - j].cov;
		if (i == opt->c.k - 1)
			kfirst[0] = kmer[0], kfirst[1] = kmer[1];
	}
	/* // The following gives better k-mer coverage towards the ends of reads, but it is a bit slow and not that useful
	// extend to the 3'-end
	if (l >= opt->c.k && v >= 0) {
		for (i = 0; i < opt->c.k - 1; ++i) {
			if (v < 0 || !fmc_cell_has_b1(v) || fmc_cell_has_b2(v)) break;
			append_to_kmer(opt->c.k, kmer, fmc_cell_get_b1(v));
			v = kmer_lookup(opt->c.k, opt->c.suf_len, kmer, h, cache);
			for (j = 0; j < opt->c.k - i - 1; ++j)
				++seq->a[seq->n - 1 - j].cov;
		}
	}
	// extend to the 5'-end
	if (seq->a[0].cov > 0) {
		kmer[0] = kfirst[1], kmer[1] = kfirst[0]; // we are extending on the reverse strand
		v = kmer_lookup(opt->c.k, opt->c.suf_len, kmer, h, cache);
		for (i = 0; i < opt->c.k - 1; ++i) {
			if (v < 0 || !fmc_cell_has_b1(v) || fmc_cell_has_b2(v)) break;
			append_to_kmer(opt->c.k, kmer, fmc_cell_get_b1(v));
			v = kmer_lookup(opt->c.k, opt->c.suf_len, kmer, h, cache);
			for (j = 0; j < opt->c.k - i - 1; ++j)
				++seq->a[j].cov;
		}
	}
	*/
	// compute n_si
	for (i = 0, n_si = 0, in_si = 0; i < seq->n; ++i) {
		if (seq->a[i].cov >= opt->c.k - FMC_SI_GAP) {
			if (!in_si) ++n_si, in_si = 1;
		} else if (in_si) in_si = 0;
	}
	if (fmc_verbose >= 6) {
		fprintf(stderr, "SI ");
		for (i = 0; i < seq->n; ++i) {
			if (i) fputc(',', stderr);
			fprintf(stderr, "%d", seq->a[i].cov);
		}
		fputc('\n', stderr);
	}
	return n_si;
}

typedef struct {
	int penalty, n_paths, n_failures;
} correct1_stat_t;

static correct1_stat_t fmc_correct1_aux(const fmc_opt_t *opt, fmc_hash_t **h, fmc_aux_t *a)
{
	echeap1_t z;
	int l, path_end[FMC_MAX_PATHS], n_paths = 0, max_i = 0, n_failures = 0;
	correct1_stat_t s;

	assert(a->ori.n < 0x10000);
	a->heap.n = a->stack.n = 0;
	s.penalty = s.n_paths = s.n_failures = 0;
	// find the first k-mer
	memset(&z, 0, sizeof(echeap1_t));
	for (z.i = 0, l = 0; z.i < a->seq.n;) {
		if (a->seq.a[z.i].b > 3) l = 0, z.kmer[0] = z.kmer[1] = 0;
		else ++l, append_to_kmer(opt->c.k, z.kmer, a->seq.a[z.i].b);
		if (++z.i == a->seq.n) break;
		if (l >= opt->c.k && kmer_lookup(opt->c.k, opt->c.suf_len, z.kmer, h, a->cache) >= 0) break;
	}
	if (z.i == a->seq.n) return s;
	z.last_solid = 0; z.ec_pos4 = 0; z.k = -1; // the first k-mer is not on the stack
	kv_push(echeap1_t, a->heap, z);
	// search for the best path
	while (a->heap.n) {
		ecbase_t *c;
		int val, is_excessive;
		z = a->heap.a[0];
		a->heap.a[0] = kv_pop(a->heap);
		ks_heapdown_ec(0, a->heap.n, a->heap.a);
		if (n_paths && z.penalty > a->stack.a[path_end[0]].penalty + opt->max_penalty_diff) break;
		if (z.i == a->seq.n) { // end of sequence
			if (fmc_verbose >= 6) fprintf(stderr, "** penalty=%d\n", z.penalty);
			path_end[n_paths++] = z.k;
			if (n_paths == FMC_MAX_PATHS) break;
			continue;
		}
		if (fmc_verbose >= 6) {
			int i;
			fprintf(stderr, "<- [%d] (%d,%c%d), size=%ld, penalty=%d, state=%d, ec_pos4=[", z.k, z.i, "ACGTN"[a->seq.a[z.i].b], a->seq.a[z.i].q,
					a->heap.n, z.penalty, z.k>=0? a->stack.a[z.k].state : -1);
			for (i = 3; i >= 0; --i)
				if (z.ec_pos4>>(i*16)&0xffff) {
					fprintf(stderr, "%d", (int)(z.ec_pos4>>(i*16)&0xffff) - 1);
					if (i) fputc(',', stderr);
				}
			fprintf(stderr, "]\n");
		}
		c = &a->seq.a[z.i];
		max_i = max_i > z.i? max_i : z.i;
		is_excessive = (a->heap.n >= max_i * 3);
		val = kmer_lookup(opt->c.k, opt->c.suf_len, z.kmer, h, a->cache);
		if (val >= 0 && fmc_cell_has_b1(val)) { // present in the hash table
			int b1 = fmc_cell_get_b1(val);
			int b2 = fmc_cell_has_b2(val)? fmc_cell_get_b2(val) : 4;
			int q1 = fmc_cell_get_q1(val);
			int q2 = fmc_cell_get_q2(val);
			int is_solid = (opt->ecQ > 0 && c->q >= opt->ecQ);
			if (b1 == c->b) { // read base matching the consensus
				update_aux(opt->c.k, a, &z, c->b, STATE_M, 0, 0, is_solid);
				if (b2 != 4 && q1 < 20 && !is_solid && !is_excessive) update_aux(opt->c.k, a, &z, b2, STATE_M, q1 > c->oq? q1 : c->oq, 1, 0);
			} else if (c->b > 3) { // read base is "N"
				update_aux(opt->c.k, a, &z, b1, STATE_M, 3, 0, 0); // "N" is not counted as a diff
				if (b2 < 4 && !is_excessive) update_aux(opt->c.k, a, &z, b2, STATE_M, q1, 0, 0);
			} else if (is_solid && (z.i < z.last_solid + opt->c.k || q1>>1 < FMC_Q_1)) { // base is solid and the evidence is weak; reject
				update_aux(opt->c.k, a, &z, c->b, STATE_N, q1, 0, q1 >= opt->ecQ? 1 : 0); // TODO: should we set is_solid???
			} else if (z.ec_pos4>>48 && z.i < (z.ec_pos4>>48) + opt->max_dist4) {
				update_aux(opt->c.k, a, &z, c->b, STATE_N, q1, 0, q1 >= opt->ecQ? 1 : 0); // TODO: should we set is_solid???
			} else if (b2 >= 4 || b2 == c->b) { // no second base or the second base is the read base; two branches
				if (!is_excessive || q1 <= c->q)
					update_aux(opt->c.k, a, &z, c->b, STATE_M, q1,   0, 0);
				if (!is_excessive || q1 >= c->q)
					update_aux(opt->c.k, a, &z, b1,   STATE_M, c->q, 1, is_solid);
			} else { // we are looking at three different bases
				if (!is_excessive || q1 + q2 <= c->q)
					update_aux(opt->c.k, a, &z, c->b, STATE_M, q1 + q2, 0, 0);
				if (!is_excessive || q1 + q2 >= c->q)
					update_aux(opt->c.k, a, &z, b1,   STATE_M, c->q,    1, is_solid);
				if (!is_excessive)
					update_aux(opt->c.k, a, &z, b2,   STATE_M, c->q > q1? c->q : q1, 1, is_solid);
			}
		} else {
			update_aux(opt->c.k, a, &z, c->b < 4? c->b : lrand48()&4, STATE_N, FMC_NOHIT_PEN, 0, 0); // not present in the hash table
			if (n_paths == 0) ++s.n_failures;
			if (++n_failures > a->seq.n && a->heap.n > 1) a->heap.n = 1, n_failures = 0;
		}
		if (fmc_verbose >= 6) fprintf(stderr, "//\n");
	}
	// backtrack
	s.n_paths = n_paths;
	if (n_paths) {
		int j;
		if (fmc_verbose >= 6) { // duplicated computation, but this is for debugging only
			for (j = 0; j < n_paths; ++j) {
				int i;
				fprintf(stderr, "%.2d ", j);
				path_backtrack(&a->stack, path_end[j], &a->seq, &a->tmp[0]);
				for (i = 0; i < a->tmp[0].n; ++i) {
					int s = a->tmp[0].a[i].state;
					fputc(s == STATE_D? '-' : s == STATE_I? "acgtn"[a->tmp[0].a[i].b] : "ACGTN"[a->tmp[0].a[i].b], stderr);
				}
				fprintf(stderr, "\t%d\n", a->stack.a[path_end[j]].penalty - a->stack.a[path_end[0]].penalty);
			}
			fprintf(stderr, "//\n");
		}
		s.penalty = a->stack.a[path_end[0]].penalty;
		path_backtrack(&a->stack, path_end[0], &a->seq, &a->tmp[0]);
		for (j = 1; j < n_paths; ++j) {
			path_backtrack(&a->stack, path_end[j], &a->seq, &a->tmp[1]);
			path_adjustq(a->stack.a[path_end[j]].penalty - a->stack.a[path_end[0]].penalty, &a->tmp[0], &a->tmp[1]);
		}
		fmc_seq_cpy_no_del(&a->seq, &a->tmp[0]);
	}
	return s;
}

int fmc_cns_ungap(ecseq_t *s1, const ecseq_t *s2)
{
	int i, n = 0;
	assert(s1->n == s2->n);
	for (i = 0; i < s1->n; ++i) {
		ecbase_t *b1 = &s1->a[i];
		const ecbase_t *b2 = &s2->a[i];
		if (b2->state == STATE_N) continue;
		if (b1->state == STATE_N) {
			*b1 = *b2;
		} else if (b1->b == b2->b) { // same correction in both directions
			b1->f = b1->f > b2->f? b1->f : b2->f;
			b1->min_diff = b1->min_diff > b2->min_diff? b1->min_diff : b2->min_diff;
		} else { // conflicting correction; don't make a correction
			b1->f = 0; b1->min_diff = 0;
			b1->b = b1->ob; b1->q = b1->oq;
			++n;
		}
	}
	return n;
}

typedef struct {
	int n_diff, q_diff, n_paths[2], n_failures[2];
	int penalty, n_conflict, n_si, to_drop;
} fmc_ecstat_t;

void fmc_correct1(const fmc_opt_t *opt, fmc_hash_t **h, char **s, char **q, fmc_aux_t *a, fmc_ecstat_t *ecs)
{
	fmc_aux_t *_a = 0;
	int i;
	correct1_stat_t st[2];

	memset(ecs, 0, sizeof(fmc_ecstat_t));
	if (a == 0) a = _a = fmc_aux_init();
	kh_clear(kache, a->cache);
	fmc_seq_conv(*s, *q, opt->defQ, &a->ori);
	// forward strand
	fmc_seq_cpy_no_del(&a->seq, &a->ori);
	st[0] = fmc_correct1_aux(opt, h, a);
	fmc_seq_cpy_no_del(&a->ec_for, &a->seq);
	// reverse strand
	fmc_seq_cpy_no_del(&a->seq, &a->ori);
	fmc_seq_revcomp(&a->seq);
	st[1] = fmc_correct1_aux(opt, h, a);
	fmc_seq_revcomp(&a->seq);
	ecs->n_conflict = fmc_cns_ungap(&a->ec_for, &a->seq);
	fmc_seq_cpy_no_del(&a->seq, &a->ec_for);
	// generate final stats
	ecs->n_paths[0] = st[0].n_paths; ecs->n_failures[0] = st[0].n_failures;
	ecs->n_paths[1] = st[1].n_paths; ecs->n_failures[1] = st[1].n_failures;
	ecs->penalty = st[0].penalty + st[1].penalty;
	if (a->seq.n > a->ori.n) {
		*s = realloc(*s, a->seq.n + 1);
		*q = realloc(*q, a->seq.n + 1);
	} else if (!*q) *q = calloc(a->seq.n + 1, 1);
	ecs->n_si = kmer_cov(opt, &a->seq, h, a->cache);
	ecs->to_drop = (ecs->n_si == 0 || ecs->n_failures[0] > a->seq.n || ecs->n_failures[1] > a->seq.n);
	// write the sequence
	for (i = 0; i < a->seq.n; ++i) {
		ecbase_t *b = &a->seq.a[i];
		int qual;
		if (b->f || a->ori.a[b->i].b < 4) { // FIXME: this block is not quite right in the presence of gaps.
			if (b->b != a->ori.a[b->i].b)
				++ecs->n_diff, ecs->q_diff += a->ori.a[b->i].q;
			(*s)[i] = b->b == a->ori.a[b->i].b? "ACGTN"[b->b] : "acgtn"[b->b];
			qual = b->f? b->min_diff : b->q;
			(*q)[i] = (qual < FMC_Q_MAX_OUT? qual : FMC_Q_MAX_OUT) + 33;
		} else (*s)[i] = 'N', (*q)[i] = 33;
	}
	(*s)[i] = (*q)[i] = 0;
	if (_a) fmc_aux_destroy(_a);
}

typedef struct {
	const fmc_opt_t *opt;
	fmc_hash_t **h;
	char **name, **s, **q;
	fmc_ecstat_t *ecs;
	fmc_aux_t **a;
	int64_t start;
} for_correct_t;

static void correct_func(void *data, long i, int tid)
{
	for_correct_t *f = (for_correct_t*)data;
	if (fmc_verbose >= 5)
		fprintf(stderr, ">%s tid:%d heap:%ld stack:%ld kmercache:%d\n",
				f->name[i], tid, f->a[tid]->heap.m, f->a[tid]->stack.m, kh_n_buckets(f->a[tid]->cache));
	fmc_correct1(f->opt, f->h, &f->s[i], &f->q[i], f->a[tid], &f->ecs[i]);
}

void fmc_correct(const fmc_opt_t *opt, fmc_hash_t **h, int n, char **s, char **q, char **name, char **last_name, int64_t *last_id)
{
	for_correct_t f;
	int i, j;
	double tr, tc;
	char *la = *last_name;
	int64_t li = *last_id;

	if (n <= 0) return;
	tr = realtime(), tc = cputime();
	f.a = calloc(opt->n_threads, sizeof(void*));
	f.opt = opt, f.h = h, f.name = name, f.s = s, f.q = q;
	f.ecs = calloc(n, sizeof(fmc_ecstat_t));
	for (i = 0; i < opt->n_threads; ++i)
		f.a[i] = fmc_aux_init();
	if (opt->n_threads == 1) {
		for (i = 0; i < n; ++i)
			correct_func(&f, i, 0);
	} else kt_for(opt->n_threads, correct_func, &f, n);
	for (i = 0; i < n; ++i) {
		fmc_ecstat_t *s = &f.ecs[i];
		int is_same = 0; // whether the current read in the same template as the last one
		int64_t id;
		char *ni = f.name[i];
		if (la) {
			for (j = 0; la[j] && ni[j] && la[j] == ni[j]; ++j);
			if ((la[j] == 0 && ni[j] == 0) || (j > 0 && isdigit(la[j]) && isdigit(ni[j]) && la[j-1] == '/' && la[j+1] == 0 && ni[j+1] == 0))
				is_same = 1;
		}
		id = is_same? li : li + 1;
		la = ni, li = id;
		if (opt->drop_reads && s->to_drop) continue;
		if (opt->show_ori_name) printf("@%s", ni);
		else printf("@%ld", (long)id);
		printf(" ec:Z:%d_%d_%d_%d_%d:%d_%d:%d\n", s->n_si, s->n_diff, s->q_diff, s->n_conflict, s->n_paths[0], s->n_paths[1], s->n_failures[0], s->n_failures[1]);
		puts(f.s[i]); putchar('+'); putchar('\n');
		puts(f.q[i]);
	}
	free(*last_name);
	*last_name = strdup(la); *last_id = li;
	for (i = 0; i < opt->n_threads; ++i)
		fmc_aux_destroy(f.a[i]);
	free(f.a); free(f.ecs);
	fprintf(stderr, "[M::%s] corrected %d reads in %.3f sec (%.3f CPU sec)\n", __func__, n, realtime() - tr, cputime() - tc);
}

/*********************
 *** Main function ***
 *********************/

int main_correct(int argc, char *argv[])
{
	int c;
	fmc_opt_t opt;
	fmc64_v *kmer;
	char *fn_kmer = 0, *last_name = 0;
	int64_t last_id = 0;

	liftrlimit();

	fmc_opt_init(&opt);
	while ((c = getopt(argc, argv, "DOk:o:t:h:v:p:e:q:w:")) >= 0) {
		if (c == 'k') opt.c.k = atoi(optarg);
		else if (c == 'd') opt.c.q1_depth = atoi(optarg);
		else if (c == 'o') opt.c.min_occ = atoi(optarg), opt.c.max_ec_depth = opt.c.min_occ - 1;
		else if (c == 'p') opt.c.prior = atof(optarg);
		else if (c == 'e') opt.c.err = atof(optarg);
		else if (c == 't') opt.n_threads = atoi(optarg);
		else if (c == 'h') fn_kmer = optarg;
		else if (c == 'v') fmc_verbose = atoi(optarg);
		else if (c == 'q') opt.ecQ = atoi(optarg);
		else if (c == 'O') opt.show_ori_name = 1;
		else if (c == 'D') opt.drop_reads = 1;
		else if (c == 'w') opt.max_dist4 = atoi(optarg);
	}
	if (!(opt.c.k&1)) {
		++opt.c.k;
		fprintf(stderr, "[W::%s] -k must be an odd number; change -k to %d\n", __func__, opt.c.k);
	}
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi2 correct [options] index.fmd [reads.fq]\n\n");
		fprintf(stderr, "Options: -t INT     number of threads [1]\n");
		fprintf(stderr, "         -k INT     k-mer length [%d]\n", opt.c.k);
		fprintf(stderr, "         -o INT     min occurrence for a solid k-mer [%d]\n", opt.c.min_occ);
		fprintf(stderr, "         -d INT     correct singletons out of INT bases [%d]\n\n", opt.c.q1_depth);
		fprintf(stderr, "         -h FILE    get solid k-mer list from FILE [null]\n");
		fprintf(stderr, "         -q INT     protect Q>INT bases unless they occur once [%d]\n", opt.ecQ);
		fprintf(stderr, "         -w INT     no more than 4 corrections per INT-bp window [%d]\n", opt.max_dist4);
		fprintf(stderr, "         -D         drop error-prone reads\n");
		fprintf(stderr, "         -O         print the original read name\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Notes: If reads.fq is absent, this command dumps the list of solid k-mers.\n");
		fprintf(stderr, "       The dump can be loaded later with option -h.\n\n");
		return 1;
	}
	opt.c.suf_len = opt.c.k > 18? opt.c.k - 18 : 1;

	if (fn_kmer) {
		FILE *fp;
		fp = strcmp(fn_kmer, "-")? fopen(fn_kmer, "rb") : stdin;
		assert(fp);
		kmer = fmc_kmer_read(fp, &opt);
		fclose(fp);
	} else kmer = fmc_collect(&opt, argv[optind]);

	if (optind + 2 > argc) {
		int i;
		fmc_kmer_write(stdout, &opt, kmer);
		for (i = 0; i < 1<<opt.c.suf_len*2; ++i)
			free(kmer[i].a);
		free(kmer);
		return 0;
	} else {
		fmc_hash_t **h;
		int i;
		kseq_t *ks;
		gzFile fp;
		fmc_batch_t *b;

		h = fmc_kmer2hash(&opt, kmer); // kmer is deallocated here
		fp = gzopen(argv[optind+1], "r");
		ks = kseq_init(fp);
		while ((b = fmc_batch_read(ks, opt.batch_size)) != 0) {
			fmc_correct(&opt, h, b->n, b->s, b->q, b->name, &last_name, &last_id);
			fmc_batch_destroy(b);
		}
		free(last_name);
		kseq_destroy(ks);
		gzclose(fp);
		for (i = 0; i < 1<<opt.c.suf_len*2; ++i)
			kh_destroy(fmc, h[i]);
		free(h);
	}
	return 0;
}

int main_inspectk(int argc, char *argv[])
{
	rld_t *e;
	int j;
	if (argc < 3) {
		fprintf(stderr, "Usage: fermi2 inspectk <index.fmd> <kmer1> [...]\n");
		return 1;
	}
	e = rld_restore_mmap(argv[1]);
	for (j = 2; j < argc; ++j) {
		int i, len;
		rldintv_t s, t[6];
		char *aj = argv[j];
		len = strlen(aj);
		s.x[0] = s.x[1] = 0; s.x[2] = e->mcnt[0];
		printf("%d\t", len);
		for (i = len - 1; i >= 0; --i) {
			int c = seq_nt6_table[(int)aj[i]];
			rld_extend(e, &s, t, 1);
			if (t[c].x[2] == 0) break;
			s = t[c];
			rld_extend(e, &s, t, 1);
		}
		printf("%d\t", len - 1 - i);
		rld_extend(e, &s, t, 1);
		for (i = 1; i <= 4; ++i) {
			if (i != 1) putchar(':');
			printf("%ld", (long)t[i].x[2]);
		}
		printf("\t%s\t", aj);
		rld_extend(e, &s, t, 0);
		for (i = 4; i >= 1; --i) {
			if (i != 4) putchar(':');
			printf("%ld", (long)t[i].x[2]);
		}
		putchar('\n');
	}
	rld_destroy(e);
	return 0;
}
