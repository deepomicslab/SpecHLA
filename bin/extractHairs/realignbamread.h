#ifndef __REALIGNBAMREAD_H__
#define __REALIGNBAMREAD_H__

#include "nw.h"
#include <assert.h>
#include <stdlib.h>
#include "bamread.h"
#include "hapfragments.h"
#include "readvariant.h"
//#include "seqan/align.h"
//
//#include "seqan/score.h"
#include "edlib.h"
//#include "ssw_cpp.h"
 // anchor sequences must have unique kmers of this length
// find the variants that are covered by the read and determine the alleles at each of those variants

extern int VERBOSE;
extern int REALIGN_ALL;
extern int VARIANTS;
extern int PARSEINDELS;
extern int* fcigarlist;
extern int SUM_ALL_ALIGN;
extern int HOMOZYGOUS;
extern std::unordered_map<std::string, std::vector<int>> SUPPORT_READS;


 // max number of SNVs in a short haplotype
// given a=log10(x) and b=log10(y), returns log10(x+y)
#define addlogs(a, b) (((a) > (b)) ? ((a) + log10(1.0 + pow(10.0, (b) - (a)))) : ((b) + log10(1.0 + pow(10.0, (a) - (b)))))
// given a=log10(x) and b=log10(y), returns log10(x-y)
#define subtractlogs(a, b) (((a) > (b)) ? ((a) + log10(1.0 - pow(10, (b) - (a)))) : ((b) + log10(1.0 - pow(10.0, (a) - (b)))))

int compare_strings(const void* a, const void* b);
unsigned int rand_interval(unsigned int min, unsigned int max);
int parse_cigar(struct alignedread* read,REFLIST* reflist,int* fcigarlist);
int test_complexity(char* seq, int k);
int realign_HAPs(struct alignedread* read, REFLIST* reflist, int positions[], VARIANT* varlist, int* snplst, int n_snps, FRAGMENT* fragment);
int compare_read_HAPs(struct alignedread* read,VARIANT* varlist,int* snplst, int n_snps, int hap_pos[], int* fcigarlist,int fcigs,int f1, int f2, REFLIST* reflist, FRAGMENT* fragment);

int realign_and_extract_variants_read(struct alignedread* read,HASHTABLE* ht,CHROMVARS* chromvars,VARIANT* varlist,int paired,FRAGMENT* fragment,int chrom,REFLIST* reflist);

float blast_score(const char* seq1, const char* seq2);

char* reverse_dna(const char* s);
char* complement_dna(const char* s);
#endif // __REALIGNBAMREAD_H__