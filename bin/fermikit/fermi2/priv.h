#ifndef FM_PRIV_H
#define FM_PRIV_H

extern int fm_verbose;
extern unsigned char seq_nt6_table[256];

double cputime(void);
double realtime(void);
void liftrlimit(void);

void seq_comp6(int l, unsigned char *s);
void seq_reverse(int l, unsigned char *s);
void seq_revcomp6(int l, unsigned char *s);

#endif
