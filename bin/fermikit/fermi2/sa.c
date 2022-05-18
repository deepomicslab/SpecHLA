#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "fermi2.h"
#include "kvec.h"

static void sa_gen1(const rld_t *e, fmsa_t *sa, int64_t k, uint64_v *buf)
{
	int c, mask = (1<<sa->ss) - 1;
	uint64_t ok[e->asize1], k0 = k, l = 0;
	size_t i;
	buf->n = 0;
	do {
		++l;
		c = rld_rank1a(e, k + 1, ok);
		k = e->cnt[c] + ok[c] - 1;
		if (c) {
			if (((k - e->mcnt[1]) & mask) == 0) {
				int64_t x = (k - e->mcnt[1]) >> sa->ss;
				sa->ssa[x] = l;
				kv_push(uint64_t, *buf, x);
			}
		} else sa->r2i[k] = k0;
	} while (c);
	for (i = 0; i < buf->n; ++i)
		sa->ssa[buf->a[i]] = (l - 1 - sa->ssa[buf->a[i]]) << sa->ms | k0;
}

typedef struct {
	const rld_t *e;
	fmsa_t *sa;
	uint64_v *buf;
} worker_t;

static void worker(void *data, long i, int tid)
{
	worker_t *w = (worker_t*)data;
	sa_gen1(w->e, w->sa, i, &w->buf[tid]);
}

fmsa_t *fm_sa_gen(const rld_t *e, int ssa_shift, int n_threads)
{
	extern void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
	fmsa_t *sa;
	worker_t *w;
	int i;

	sa = calloc(1, sizeof(fmsa_t));
	sa->ss = ssa_shift;
	sa->m = e->mcnt[1];
	for (sa->ms = 1; 1LL<<sa->ms < sa->m; ++sa->ms);
	sa->n_ssa = (e->mcnt[0] - e->mcnt[1] + (1<<sa->ss) - 1) >> sa->ss;
	sa->r2i = calloc(sa->m, 8);
	sa->ssa = calloc(sa->n_ssa, 8);

	w = calloc(1, sizeof(worker_t));
	w->buf = calloc(n_threads, sizeof(uint64_v));
	w->sa = sa;
	w->e = e;

	kt_for(n_threads, worker, w, sa->m);

	for (i = 0; i < n_threads; ++i) free(w->buf[i].a);
	free(w->buf); free(w);
	return sa;
}

void fm_sa_destroy(fmsa_t *sa)
{
	free(sa->r2i); free(sa->ssa); free(sa);
}

int64_t fm_sa(const rld_t *e, const fmsa_t *sa, int64_t k, int64_t *si)
{
	int c, mask = (1<<sa->ss) - 1;
	int64_t x = 0;
	uint64_t ok[e->asize1];
	*si = -1;
	if (k >= e->mcnt[0]) return -1;
	while (k < e->mcnt[1] || ((k - e->mcnt[1]) & mask)) {
		++x;
		c = rld_rank1a(e, k + 1, ok);
		k = e->cnt[c] + ok[c] - 1;
		if (c == 0) {
			*si = sa->r2i[k];
			return x - 1;
		}
	}
	k = (k - e->mcnt[1]) >> sa->ss;
	*si = sa->ssa[k] & ((1ULL<<sa->ms) - 1);
	return x + (sa->ssa[k] >> sa->ms);
}

int fm_sa_dump(const fmsa_t *sa, const char *fn)
{
	uint32_t y;
	FILE *fp;
	fp = fn && strcmp(fn, "-")? fopen(fn, "wb") : fdopen(fileno(stdout), "wb");
	if (fp == 0) return -1;
	y = sa->ss; fwrite(&y, 4, 1, fp);
	y = sa->ms; fwrite(&y, 4, 1, fp);
	fwrite(&sa->m, 8, 1, fp);
	fwrite(&sa->n_ssa, 8, 1, fp);
	fwrite(sa->r2i, 8, sa->m, fp);
	fwrite(sa->ssa, 8, sa->n_ssa, fp);
	fclose(fp);
	return 0;
}

fmsa_t *fm_sa_restore(const char *fn)
{
	FILE *fp;
	uint32_t y;
	fmsa_t *sa;
	fp = fn && strcmp(fn, "-")? fopen(fn, "rb") : fdopen(fileno(stdin), "rb");
	if (fp == 0) return 0;
	sa = calloc(1, sizeof(fmsa_t));
	fread(&y, 4, 1, fp); sa->ss = y;
	fread(&y, 4, 1, fp); sa->ms = y;
	fread(&sa->m, 8, 1, fp);
	fread(&sa->n_ssa, 8, 1, fp);
	sa->r2i = calloc(sa->m, 8);
	sa->ssa = calloc(sa->n_ssa, 8);
	if (sa->ssa == 0 || sa->r2i == 0) {
		free(sa->r2i); free(sa->ssa); free(sa);
		return 0;
	}
	fread(sa->r2i, 8, sa->m, fp);
	fread(sa->ssa, 8, sa->n_ssa, fp);
	fclose(fp);
	return sa;
}

int main_sa(int argc, char *argv[])
{
	int c, n_threads = 1, ssa_shift = 6;
	fmsa_t *sa;
	rld_t *e;
	char *fn = 0;

	while ((c = getopt(argc, argv, "t:s:o:")) >= 0) {
		if (c == 't') n_threads = atoi(optarg);
		else if (c == 's') ssa_shift = atoi(optarg);
		else if (c == 'o') fn = optarg;
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fermi2 sa [-t nThreads=%d] [-s stepShift=%d] <in.fmd>\n", n_threads, ssa_shift);
		return 1;
	}
	e = rld_restore(argv[optind]);
	if (e == 0) {
		fprintf(stderr, "[E::%s] failed to load the FM-index\n", __func__);
		return 1;
	}
	sa = fm_sa_gen(e, ssa_shift, n_threads);
	//{int64_t i;for (i = 0; i < sa->n_ssa; ++i) fprintf(stderr, "S(%lld)=%lld\n", (i<<sa->ss) + e->mcnt[1], sa->ssa[i]>>sa->ms);}

	fm_sa_dump(sa, fn);
	rld_destroy(e);
	fm_sa_destroy(sa);
	return 0;
}
