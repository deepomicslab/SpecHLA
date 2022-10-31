#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include "kstring.h"
#include "kvec.h"
#include "rld0.h"

static int dfs_verbose = 3;

/******************
 *** DFS engine ***
 ******************/

#define DFS_SUF_LEN 5

typedef struct {
	int64_t k, l;
	int d, c;
} elem_t;

typedef struct {
	uint64_t c[6];
} fmint6_t;

typedef void (*fmdfs_f)(void *data, int tid, int k, char *path, const fmint6_t *size, int *cont);
typedef void (*fmdfs2_f)(void *data, int tid, int k, char *path, const rldintv_t *ik, const rldintv_t *ok, int *cont);

void fm_dfs_core(int n, rld_t *const*e, int is_half, int max_k, int suf_len, int suf, fmdfs_f func, void *data, int tid)
{ // this routine is similar to fmc_collect1()
	int i, j, c;
	fmint6_t *size, *tk, *tl;
	elem_t *t;
	char *_path, *path;
	uint64_t ok[6], ol[6];
	kvec_t(elem_t) stack = {0,0,0};

	assert((max_k&1) || !is_half);
	t = alloca(sizeof(elem_t) * n);
	size = alloca(sizeof(fmint6_t) * n);
	tk = alloca(sizeof(fmint6_t) * n);
	tl = alloca(sizeof(fmint6_t) * n);
	_path = alloca(max_k + 2);
	path = _path + 1;
	kv_resize(elem_t, stack, n * max_k * 6);
	// descend
	for (i = 0; i < n; ++i) {
		elem_t *p;
		kv_pushp(elem_t, stack, &p);
		p->k = 0, p->l = e[i]->mcnt[0];
		for (j = 0; j < suf_len; ++j) {
			c = (suf>>j*2&3) + 1;
			rld_rank2a(e[i], p->k, p->l, ok, ol);
			p->k = e[i]->cnt[c] + ok[c];
			p->l = e[i]->cnt[c] + ol[c];
		}
		p->d = suf_len, p->c = (suf>>(suf_len-1)*2&3) + 1;
	}
	for (j = 0; j < suf_len; ++j)
		path[max_k - j - 1] = "ACGT"[suf>>j*2&3];
	path[max_k] = 0;
	// traverse
	while (stack.n) {
		int end, cont = 0x1E;
		for (i = n - 1; i >= 0; --i) t[i] = kv_pop(stack);
		if (t->d > max_k) continue;
		path[max_k - t->d] = "\0ACGTN"[t->c];
		for (i = 0; i < n; ++i) {
			assert(t[i].k < e[i]->mcnt[0]);
			rld_rank2a(e[i], t[i].k, t[i].l, tk[i].c, tl[i].c);
			for (c = 0; c < 6; ++c) {
				size[i].c[c] = tl[i].c[c] - tk[i].c[c];
				if (size[i].c[c] == 0) cont &= ~(1<<c);
			}
		}
		func(data, tid, t->d, path + (max_k - t->d), size, &cont);
		end = t->d == max_k>>1 && is_half? 2 : 4;
		for (c = 1; c <= end; ++c) {
			if ((cont>>c&1) == 0) continue;
			for (i = 0; i < n; ++i) {
				elem_t *p;
				kv_pushp(elem_t, stack, &p);
				p->k = e[i]->cnt[c] + tk[i].c[c];
				p->l = e[i]->cnt[c] + tl[i].c[c];
				p->d = t->d + 1;
				p->c = c;
			}
		}
	}
	free(stack.a);
}

void fm_dfs2_core(int n, rld_t *const*e, int is_half, int max_k, int suf_len, int suf, fmdfs2_f func, void *data, int tid)
{ // this is the bidirectional version of fm_dfs_core(), requiring a bidirectional FM-index
	int i, j, c;
	rldintv_t *t, *o;
	char *_path, *path;
	kvec_t(rldintv_t) stack = {0,0,0};

	// check if the input is correct
	assert((max_k&1) || !is_half);
	for (i = 0; i < n; ++i) // check bidirectionality
		assert(e[i]->mcnt[2] == e[i]->mcnt[5] && e[i]->mcnt[3] == e[i]->mcnt[4]);
	// allocation
	t = alloca(sizeof(rldintv_t) * n);
	o = alloca(sizeof(rldintv_t) * n * 6);
	_path = alloca(max_k + 2);
	path = _path + 1;
	kv_resize(rldintv_t, stack, n * max_k * 6);
	// descend
	for (i = 0; i < n; ++i) {
		rldintv_t *p, t[6];
		kv_pushp(rldintv_t, stack, &p);
		p->x[0] = p->x[1] = p->info = 0, p->x[2] = e[i]->mcnt[0];
		for (j = 0; j < suf_len; ++j) {
			rld_extend(e[i], p, t, 1);
			*p = t[(suf>>j*2&3) + 1];
		}
		p->info = suf_len<<8 | ((suf>>(suf_len-1)*2&3) + 1);
	}
	for (j = 0; j < suf_len; ++j)
		path[max_k - j - 1] = "ACGT"[suf>>j*2&3];
	path[max_k] = 0;
	// traverse
	while (stack.n) {
		int end, cont = 0x1E, depth;
		rldintv_t *oi;
		for (i = n - 1; i >= 0; --i) t[i] = kv_pop(stack);
		depth = t->info>>8;
		if (depth > max_k) continue;
		path[max_k - depth] = "\0ACGTN"[t->info&0xff];
		for (i = 0, oi = o; i < n; ++i, oi += 6) {
			rld_extend(e[i], &t[i], oi, 1);
			for (c = 0; c < 6; ++c)
				if (oi[c].x[2] == 0) cont &= ~(1<<c);
		}
		func(data, tid, depth, path + (max_k - depth), t, o, &cont);
		end = depth == max_k>>1 && is_half? 2 : 4;
		for (c = 1; c <= end; ++c) {
			if ((cont>>c&1) == 0) continue;
			for (i = 0, oi = o; i < n; ++i, oi += 6) {
				rldintv_t *p;
				kv_pushp(rldintv_t, stack, &p);
				*p = oi[c];
				p->info = (depth + 1) << 8 | c;
			}
		}
	}
	free(stack.a);
}

typedef struct {
	int n, max_k, suf_len, is_half;
	rld_t *const*e;
	void *data;
	fmdfs_f func;
	fmdfs2_f func2;
} shared_t;

static void dfs_worker(void *data, long suf, int tid)
{
	shared_t *d = (shared_t*)data;
	if (d->func) fm_dfs_core(d->n, d->e, d->is_half, d->max_k, d->suf_len, suf, d->func, d->data, tid);
	if (d->func2) fm_dfs2_core(d->n, d->e, d->is_half, d->max_k, d->suf_len, suf, d->func2, d->data, tid);
	if (dfs_verbose >= 4)
		fprintf(stderr, "[M::%s] processed suffix %ld in thread %d\n", __func__, suf, tid);
}

void fm_dfs(int n, rld_t *const*e, int is_half, int max_k, int n_threads, fmdfs_f func, fmdfs2_f func2, void *data)
{
	extern void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
	shared_t d;
	int n_suf;
	d.n = n, d.e = e, d.data = data, d.func = func, d.func2 = func2, d.max_k = max_k; d.is_half = is_half;
	d.suf_len = max_k>>1 < DFS_SUF_LEN? max_k>>1 : DFS_SUF_LEN;
	n_suf = 1<<d.suf_len*2;
	n_threads = n_threads < n_suf? n_threads : n_suf;
	kt_for(n_threads, dfs_worker, &d, n_suf);
}

/*************
 *** Count ***
 *************/

typedef struct {
	const rld_t *e;
	int len, min_occ, bidir, bifur_only;
	kstring_t *str;
} dfs_count_t;

static void dfs_count(void *data, int tid, int k, char *path, const fmint6_t *size, int *cont)
{
	dfs_count_t *d = (dfs_count_t*)data;
	int c;
	uint64_t sum = 0;
	for (c = 0; c < 6; ++c) {
		if (size->c[c] < d->min_occ) *cont &= ~(1<<c);
		sum += size->c[c];
	}
	if (k < d->len) return;
	printf("%s\t%ld\n", path, (long)sum);
}

static void dfs_count2(void *data, int tid, int k, char *path, const rldintv_t *ik, const rldintv_t *ok, int *cont)
{
	dfs_count_t *d = (dfs_count_t*)data;
	int c;
	rldintv_t rk[6];
	kstring_t *s = &d->str[tid];
	for (c = 0; c < 6; ++c)
		if (ok[c].x[2] < d->min_occ) *cont &= ~(1<<c);
	if (k < d->len) return;
	rld_extend(d->e, ik, rk, 0);
	if (d->bifur_only) { // check bifurcation
		int n[2];
		n[0] = n[1] = 0;
		for (c = 1; c <= 4; ++c)
			if (rk[c].x[2]) ++n[0];
		for (c = 1; c <= 4; ++c)
			if (ok[c].x[2]) ++n[1];
		if (n[0] < 2 && n[1] < 2) return; // no bifurcation; don't print
	}
	s->l = 0;
	for (c = 0; c < 6; ++c) {
		if (c) kputc(':', s);
		kputl(ok[c].x[2], s);
	}
	kputc('\t', s); kputs(path, s); kputc('\t', s);
	kputl(rk[0].x[2], s);
	for (c = 4; c >= 1; --c) {
		kputc(':', s);
		kputl(rk[c].x[2], s);
	}
	kputc(':', s); kputl(rk[5].x[2], s);
	puts(s->s);
}

int main_count(int argc, char *argv[])
{
	int i, c, n_threads = 1;
	dfs_count_t d;
	rld_t *e;
	memset(&d, 0, sizeof(dfs_count_t));
	d.len = 51, d.min_occ = 1;
	while ((c = getopt(argc, argv, "2bk:o:t:")) >= 0) {
		if (c == 'k') d.len = atoi(optarg);
		else if (c == 'o') d.min_occ = atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == '2') d.bidir = 1;
		else if (c == 'b') d.bifur_only = d.bidir = 1;
	}
	if (d.bifur_only && d.min_occ < 2) d.min_occ = 2; // in the -b mode, we need to see at least 2 k-mers
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi2 count [options] <in.fmd>\n\n");
		fprintf(stderr, "Options: -k INT      k-mer length [%d]\n", d.len);
		fprintf(stderr, "         -o INT      min occurence [%d]\n", d.min_occ);
		fprintf(stderr, "         -t INT      number of threads [%d]\n", n_threads);
		fprintf(stderr, "         -b          only print bifurcating k-mers (force -2)\n");
		fprintf(stderr, "         -2          bidirectional counting\n");
		fprintf(stderr, "\n");
		return 1;
	}
	d.str = calloc(n_threads, sizeof(kstring_t));
	d.e = e = rld_restore(argv[optind]);
	if (!(d.len&1)) {
		++d.len;
		if (dfs_verbose >= 2)
			fprintf(stderr, "[W::%s] %d is an even number; change k to %d\n", __func__, d.len-1, d.len);
	}
	if (d.bidir) fm_dfs(1, &e, 1, d.len, n_threads, 0, dfs_count2, &d);
	else fm_dfs(1, &e, 1, d.len, n_threads, dfs_count, 0, &d);
	rld_destroy(e);
	for (i = 0; i < n_threads; ++i) free(d.str[i].s);
	free(d.str);
	return 0;
}
