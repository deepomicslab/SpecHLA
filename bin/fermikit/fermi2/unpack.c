#include <stdio.h>
#include <zlib.h>
#include "rld0.h"
#include "kstring.h"
#include "kseq.h"
KSTREAM_DECLARE(gzFile, gzread)

int64_t fm_retrieve(const rld_t *e, uint64_t x, kstring_t *s)
{
	uint64_t k = x, *ok;
	ok = alloca(8 * e->asize);
	s->l = 0;
	while (1) {
		int c = rld_rank1a(e, k + 1, ok);
		k = e->cnt[c] + ok[c] - 1;
		if (c == 0) return k;
		kputc(c, s);
	}
}

char **hts_readlines(const char *fn, int64_t *_n)
{
	int64_t m = 0, n = 0;
	int dret;
	char **s = 0;
	gzFile fp;
	if ((fp = gzopen(fn, "r")) != 0) { // read from file
		kstream_t *ks;
		kstring_t str;
		str.s = 0; str.l = str.m = 0;
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
			if (str.l == 0) continue;
			if (m == n) {
				m = m? m<<1 : 16;
				s = (char**)realloc(s, m * sizeof(void*));
			}
			s[n++] = strdup(str.s);
		}
		ks_destroy(ks);
		gzclose(fp);
		s = (char**)realloc(s, n * sizeof(void*));
		free(str.s);
	} else if (*fn == ':') { // read from string
		const char *q, *p;
		for (q = p = fn + 1;; ++p)
			if (*p == ',' || *p == 0) {
				if (m == n) {
					m = m? m<<1 : 16;
					s = (char**)realloc(s, m * sizeof(void*));
				}
				s[n] = (char*)calloc(p - q + 1, 1);
				strncpy(s[n++], q, p - q);
				q = p + 1;
				if (*p == 0) break;
			}
	} else return 0;
	s = (char**)realloc(s, n * sizeof(void*));
	*_n = n;
	return s;
}

#include <unistd.h>

static void print1(const rld_t *e, int i, kstring_t *s)
{
	int64_t k;
	int j, tmp;
	k = fm_retrieve(e, i, s);
	for (j = 0; j < s->l; ++j)
		s->s[j] = "$ACGTN"[(int)s->s[j]];
	for (j = 0; j < s->l>>1; ++j) // reverse
		tmp = s->s[j], s->s[j] = s->s[s->l-1-j], s->s[s->l-1-j] = tmp;
	fwrite(s->s, 1, s->l, stdout);
	printf("\t%ld\n", (long)k);
}

int main_unpack(int argc, char *argv[])
{
	int64_t i, n;
	int c;
	rld_t *e;
	kstring_t str = {0,0,0};

	while ((c = getopt(argc, argv, "")) >= 0);
	if (optind == argc) {
		fprintf(stderr, "Usage: fermi2 unpack <reads.rld> [list|file]\n");
		return 1;
	}
	e = rld_restore(argv[optind]);
	if (optind + 1 < argc) {
		char *p, **list;
		list = hts_readlines(argv[optind+1], &n);
		if (list == 0) return 1;
		for (i = 0; i < n; ++i) {
			if (isdigit(list[i][0])) {
				int64_t x = strtol(list[i], &p, 10);
				if (x < e->mcnt[1]) print1(e, x, &str);
			}
			free(list[i]);
		}
		free(list);
	} else for (i = 0; i < e->mcnt[1]; ++i) print1(e, i, &str);
	free(str.s);
	rld_destroy(e);
	return 0;
}
