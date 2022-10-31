#ifndef FERMI2_H
#define FERMI2_H

#include <stdint.h>
#include "rld0.h"

#define FM_VERSION "r178"

#define fmd_comp(a) ((a) >= 1 && (a) <= 4? 5 - (a) : (a))
#define fmd_set_intv(e, c, ik) ((ik).x[0] = (e)->cnt[(int)(c)], (ik).x[2] = (e)->cnt[(int)(c)+1] - (e)->cnt[(int)(c)], (ik).x[1] = (e)->cnt[fmd_comp(c)], (ik).info = 0)
#define fmd_empty_intv(e, ik) ((ik).x[0] = (ik).x[1] = 0, (ik).x[2] = (e)->mcnt[0], (ik).info = 0)

typedef struct {
	rldintv_t ik, ok[2][6];
} fmdsmem_t;

typedef struct { size_t n, m; rldintv_t *a; } rldintv_v;
typedef struct { size_t n, m; fmdsmem_t *a; } fmdsmem_v;
typedef struct { size_t n, m; uint64_t  *a; } uint64_v;

typedef struct {
	int ms;
	int ss;
	int64_t m, n_ssa;
	uint64_t *r2i; // rank -> index
	uint64_t *ssa; // sampled suffix array
} fmsa_t;

#ifdef __cplusplus
extern "C" {
#endif

fmsa_t *fm_sa_gen(const rld_t *e, int ssa_shift, int n_threads);
int fm_sa_dump(const fmsa_t *sa, const char *fn);
fmsa_t *fm_sa_restore(const char *fn);
void fm_sa_destroy(fmsa_t *sa);

int64_t fm_sa(const rld_t *e, const fmsa_t *sa, int64_t k, int64_t *si);

#ifdef __cplusplus
}
#endif

#endif
