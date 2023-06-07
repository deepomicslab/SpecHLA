#ifndef __PRINT_CLUSTER_H__
#define __PRINT_CLUSTER_H__

#include "bamread.h"
#include "hashtable.h"
#include "hapfragments.h"

extern int MIN_CLUSTER_DISTANCE;

int get_chrom_name(struct alignedread* read, HASHTABLE* ht, REFLIST* reflist);

void find_matepair(struct alignedread** readlist, int s, int e) ;

int generate_single_fragment(struct alignedread** readlist, int s, int e, int length, double read_density, FRAGMENT* flist, VARIANT* varlist);


int print_read(struct alignedread** readlist, int i, int prevpos, FRAGMENT* flist, VARIANT* varlist);


void print_reads_window(struct alignedread** readlist, int s, int e, FRAGMENT* flist, VARIANT* varlist, int if_variant_read);

int print_clusters(struct alignedread** readlist, int s, int e, FRAGMENT* flist, VARIANT* varlist);

int init_clusters(struct alignedread** readlist, int s, int e);

#endif // __PRINT_CLUSTER_H__