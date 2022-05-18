#include <zlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include "kstring.h"
#include "kseq.h"
KSEQ_INIT2(, gzFile, gzread)

unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

void seq_char2nt6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l; ++i)
		s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5;
}

void seq_reverse(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		s[l-1-i] = s[i]; s[i] = tmp;
	}
}

void seq_comp6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l; ++i)
		s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}

void seq_revcomp6(int l, unsigned char *s)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
		s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		s[i] = tmp;
	}
	if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}

static void write_seq(const kseq_t *seq, kstring_t *out)
{
	kputc(seq->qual.l? '@' : '>', out);
	kputsn(seq->name.s, seq->name.l, out);
	if (seq->comment.l) {
		kputc(' ', out);
		kputsn(seq->comment.s, seq->comment.l, out);
	}
	kputc('\n', out);
	kputsn(seq->seq.s, seq->seq.l, out);
	if (seq->qual.l) {
		kputsn("\n+\n", 3, out);
		kputsn(seq->qual.s, seq->qual.l, out);
	}
	kputc('\n', out);
}

int main_interleave(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *seq[2];
	kstring_t str;

	if (argc < 3) {
		fprintf(stderr, "Usage: fermi interleave <in1.fq> <in2.fq>\n");
		return 1;
	}
	str.l = str.m = 0; str.s = 0;
	fp1 = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	fp2 = strcmp(argv[2], "-")? gzopen(argv[2], "r") : gzdopen(fileno(stdin), "r");
	seq[0] = kseq_init(fp1);
	seq[1] = kseq_init(fp2);
	while (kseq_read(seq[0]) >= 0) {
		if (kseq_read(seq[1]) < 0) break; // one file ends
		str.l = 0;
		if (seq[0]->name.l > 2 && seq[0]->name.s[seq[0]->name.l-2] == '/' && isdigit(seq[0]->name.s[seq[0]->name.l-1]))
			seq[0]->name.s[(seq[0]->name.l -= 2)] = 0; // trim tailing "/[0-9]$"
		seq[1]->name.l = 0;
		kputsn(seq[0]->name.s, seq[0]->name.l, &seq[1]->name); // make sure two ends having the same name
		write_seq(seq[0], &str);
		write_seq(seq[1], &str);
		fputs(str.s, stdout);
	}
	kseq_destroy(seq[0]); gzclose(fp1);
	kseq_destroy(seq[1]); gzclose(fp2);
	free(str.s);
	return 0;
}
