#ifndef _READVARIANT_H
#define _READVARIANT_H
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include "hashtable.h"
#include "readfasta.h"
//#define _GNU_SOURCE

#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <unordered_map>
#include <map>
#include <string>
#include "vector"


extern FILE* fragment_file; // FILE to which the fragments will be output, if NULL, output to stdout

extern int TRI_ALLELIC;

extern int BSIZE;
extern int PRINT_FRAGMENTS;
extern bool VCF_PHASED;
extern int MATE_AT_SAME;
extern int MINBNDIS;
extern char* SUPPORT_READS_TAG;
extern char* REF_READS_TAG;
//int VARIANTS = 0;

#define BNDTYPE_PAIRED 0
#define BNDTYPE_SINGLE_END 1
#define BNDTYPE_MUTLIPLE_MATE 2
#define BNDTYPE_INTRA_CHROMOSOME 3
#define BND_INS 4
#define BND_INV 5
#define BND_DEL 6
#define BND_DUP 7

#define BNDDIRECT_TRP 0
#define BNDDIRECT_TLP 1 
#define BNDDIRECT_PLT 2
#define BNDDIRECT_PRT 3

extern int BND_RANGE;
extern int BLAST_REGION_LEN;
extern int SAMPLE_IDX;

extern std::unordered_map<std::string, std::vector<int>> SUPPORT_READS;
extern std::unordered_map<std::string, std::vector<int>> REF_READS;

typedef struct {
    char* id;           //variant id
    char* chrom;        //chromomsome 
    int position;       //pos
    int ss; // pos in vcf file
    short altalleles;   //alt allele no.
    char* RA;  // reference alleles
    char* AA; // alternate alleles
    double* GLL; // genotype likelihoods added 11/25/13
    char* genotype; // encoded as integers 0 1 2 3 4 5 6 7
    short type;
    // changed this to char* on April 3 2012
    char* allele1;
    char* allele2; // temporary for SNPs
    char heterozygous; // only heterozygous variants will be used for printing out HAIRS
    int depth;
    int A1, A2;
    int H1, H2;
    int shift; // for indels in low-complexity sequence, rightshift in variant posiiton
    int bnd;
    int bnd_pos; // pos for BND
    int bnd_pair_distance;  // distance between a pair of bnd 
    int bnd_type;
    int bnd_direction;
    int bnd_mate_pos;
    int bnd_sv_len;
    char * bnd_ins_seq;
    char * bnd_mate_chrom;
    char* bnd_seq; // sequence around bnd site
    char* ref_seq;
    std::map<int, int> *snp0_dup_region;
    std::map<int, int> *snp1_dup_region;
    // total reads covering this variant (haploid/diploid, A1-> reads supporting reference allele (single-read)
    //	float L11,L12,L22; // genotype likelihoods for three possible genotypes
    int phase_set; // phase_set for phased vcf
} VARIANT;

// information about the variants on each chromosome

typedef struct {
    int variants;
    int first;
    int last;
    int blocks;
    int* intervalmap;
} CHROMVARS;

/*
typedef struct
{
        char allele; char qv; int varid;  // allele is 0/1    varid is index to varlist[varid] gives all information about the variant
} allele;

typedef struct
{
        char* id; 	int variants; allele* alist;
        int blocks; int paired; int matepos;

} FRAGMENT;

int compare_fragments(const void *a,const void *b);
 */

int count_variants(char* vcffile, char* sampleid, int* samplecol);

int count_variants_hts(char* vcffile, char* sampleid, int* samplecol);

int parse_variant(VARIANT* variant, char* buffer, int samplecol);

int parse_variant_hts(VARIANT* variant, const bcf1_t * record, const bcf_hdr_t *header, int chromosome);

int read_variantfile(char* vcffile, VARIANT* varlist, HASHTABLE* ht, int* hetvariants, int samplecol);

int read_variantfile_hts(char *vcffile, VARIANT *varlist, HASHTABLE *ht, int *hetvariants, std::unordered_map<std::string,std::pair<int, int>>& BNDs);

void build_intervalmap(CHROMVARS* chromvars, int chromosomes, VARIANT* varlist, int variants);

int calculate_rightshift(VARIANT* varlist, int ss, REFLIST* reflist);

int parse_bnd(VARIANT* variant, int chromosome);

void bnd_to_ref_seq(VARIANT *variant, REFLIST* reflist, int chromosome);
#endif
