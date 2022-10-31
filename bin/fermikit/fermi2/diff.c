#include <assert.h>
#include <string.h>
#include "rld0.h"
#include "kvec.h"

#define SUF_LEN 5 

typedef struct { size_t n, m; rldintv_t *a; } rldintv_v;

static rldintv_t descend(const rld_t *e, int suf_len, int suf)
{
	int i;
	rldintv_t ok[6], ik;
	ik.x[0] = ik.x[1] = 0; ik.x[2] = e->mcnt[0];
	for (i = 0; i < suf_len; ++i) {
		rld_extend(e, &ik, ok, 1);
		ik = ok[(suf>>i*2&3) + 1];
	}
	return ik;
}

static void collect_tips(const rld_t *e, int min_k, int min_occ, uint64_t *sub, const rldintv_t *_ik, rldintv_v *stack)
{
	stack->n = 0;
	kv_push(rldintv_t, *stack, *_ik);
	while (stack->n) {
		uint64_t k, *p, x;
		rldintv_t ik, ok[6];
		int c;
		ik = kv_pop(*stack);
		if (ik.info <= min_k && ik.x[2] < min_occ) continue;
		rld_extend(e, &ik, ok, 1);
		if (ok[0].x[2] && ik.info >= min_k) {
			for (k = 0; k < ok[0].x[2]; ++k) {
				x = k + ok[0].x[0];
				p = sub + (x>>6);
				x = 1ULL<<(x&0x3f);
				__sync_or_and_fetch(p, x);
			}
		}
		for (c = 1; c <= 4; ++c)
			if (ok[c].x[2]) {
				ok[c].info = ik.info + 1;
				kv_push(rldintv_t, *stack, ok[c]);
			}
	}
}

static void contrast_core(const rld_t *eqry, const rld_t *eref, uint64_t *sub, int min_k, int max_k, int min_occ, int suf_len, int suf)
{
	rldintv_v stack[2], tstack;
	rldintv_t ik[2], ok[2][6];
	const rld_t *e[2];
	int i, c;

	e[0] = eref; e[1] = eqry;
	kv_init(tstack); // temporary stack
	for (i = 0; i < 2; ++i) {
		kv_init(stack[i]);
		ik[i] = descend(e[i], suf_len, suf);
		ik[i].info = suf_len;
	}
	if (ik[0].x[2] == 0 || ik[1].x[2] < min_occ) return;
	kv_push(rldintv_t, stack[0], ik[0]);
	kv_push(rldintv_t, stack[1], ik[1]);
	while (stack[0].n) { // stack[0] and stack[1] are always of the same size
		ik[0] = kv_pop(stack[0]);
		ik[1] = kv_pop(stack[1]); // it is always true that ik[1].x[2] >= min_occ
		if (ik[0].x[2] == 0) collect_tips(e[1], min_k, min_occ, sub, &ik[1], &tstack);
		else if (ik[1].info >= max_k) continue;
		else {
			rld_extend(e[0], &ik[0], ok[0], 1);
			rld_extend(e[1], &ik[1], ok[1], 1);
			for (c = 1; c <= 4; ++c) {
				if (ok[1][c].x[2] < min_occ) continue;
				ok[1][c].info = ik[1].info + 1;
				kv_push(rldintv_t, stack[0], ok[0][c]);
				kv_push(rldintv_t, stack[1], ok[1][c]);
			}
		}
	}
	free(stack[0].a); free(stack[1].a); free(tstack.a);
}

static void occflt_core(const rld_t *e, uint64_t *sub, int min_k, int max_k, int min_occ, int suf_len, int suf)
{
	rldintv_v stack, tstack;
	rldintv_t ik, ok[6];

	kv_init(tstack); kv_init(stack);
	ik = descend(e, suf_len, suf);
	ik.info = suf_len;
	kv_push(rldintv_t, stack, ik);
	while (stack.n) {
		ik = kv_pop(stack);
		if (ik.info < max_k) {
			int c;
			rld_extend(e, &ik, ok, 1);
			for (c = 1; c <= 4; ++c) {
				if (ok[c].x[2] == 0) continue;
				if (ok[c].x[2] < min_occ) {
					collect_tips(e, min_k, min_occ, sub, &ok[c], &tstack);
					continue;
				}
				ok[c].info = ik.info + 1;
				kv_push(rldintv_t, stack, ok[c]);
			}
		}
	}
	free(stack.a); free(tstack.a);
}

typedef struct {
	int min_k, max_k, min_occ;
	const rld_t *eref, *eqry;
	uint64_t *sub;
} shared_t;

static void worker(void *data, long i, int tid)
{
	shared_t *d = (shared_t*)data;
	if (!d->eref) occflt_core(d->eqry, d->sub, d->min_k, d->max_k, d->min_occ, SUF_LEN, i);
	else contrast_core(d->eqry, d->eref, d->sub, d->min_k, d->max_k, d->min_occ, SUF_LEN, i);
}

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

uint64_t *fm_diff(const rld_t *eqry, const rld_t *eref, int min_k, int max_k, int min_occ, int n_threads)
{
	shared_t d;
	assert(min_k > SUF_LEN && min_k < max_k);
	d.sub = calloc((eqry->mcnt[1] + 63) / 64, 8);
	d.min_k = min_k, d.max_k = max_k, d.min_occ = min_occ, d.eref = eref, d.eqry = eqry;
	kt_for(n_threads, worker, &d, 1<<SUF_LEN*2);
	return d.sub;
}

#include <unistd.h>

int main_diff(int argc, char *argv[])
{
	int c, min_k = 25, max_k = 51, min_occ = 2, n_threads = 1;
	uint64_t n_seqs, *bits;
	rld_t *eqry = 0, *eref = 0;
	while ((c = getopt(argc, argv, "k:K:o:t:")) >= 0) {
		if (c == 'k') min_k = atoi(optarg);
		else if (c == 'K') max_k = atoi(optarg);
		else if (c == 'o') min_occ = atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
	}
	if (optind == argc) {
		if (strcmp(argv[0], "diff") == 0)
			fprintf(stderr, "Usage: fermi2 diff [-k minK=%d] [-K maxK=%d] [-o minOcc=%d] [-t nThreads=1] <query.rld> <ref.rld>\n", min_k, max_k, min_occ);
		else fprintf(stderr, "Usage: fermi2 %s [-k minK=%d] [-K maxK=%d] [-o minOcc=%d] [-t nThreads=1] <query.rld>\n", argv[0], min_k, max_k, min_occ);
		return 1;
	}
	eqry = rld_restore(argv[optind]);
	if (optind + 1 < argc)
		eref = rld_restore(argv[optind+1]);
	n_seqs = eqry->mcnt[1];
	bits = fm_diff(eqry, eref, min_k, max_k, min_occ, n_threads);
	rld_destroy(eref);
	rld_destroy(eqry);
	fwrite(&n_seqs, 8, 1, stdout);
	fwrite(bits, 8, (n_seqs + 63) / 64, stdout);
	free(bits);
	return 0;
}
