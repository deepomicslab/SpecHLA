#ifndef __FOSMIDBAM_HAIRS_H__
#define __FOSMIDBAM_HAIRS_H__


#include <stdint.h>
#include "print_clusters.h"



#define BS1 20

struct BLOCK_READ {
    int index;
    int cluster;
    int firstread;
    int lastread;
    short reads;
    int start, end;
    float GC; // GC percentage of window // GC content of the full fragment is important |  both GC-rich fragments and AT-rich fragments are underrepresented 
    float mappability; // mappability of short reads in this window 
    short subcounts[BS1];

    double bscore;
    int previous;
    int previous_cluster;
    int reads_window;
    double score_window;
    double mean, variance;
    double bscore1;
    int previous1;
    double score_window1;
    // bed file with 100bp windows -> read count for unique and non-uniquely mapping reads 
};

int cluster_reads(struct alignedread** readlist, int s, int e, FRAGMENT* flist, VARIANT* varlist, REFLIST* reflist);

void process_chunk(struct alignedread** readlist, int s, int e, FRAGMENT* flist, VARIANT* varlist, REFLIST* reflist);

int parse_bamfile_fosmid(char* bamfile, HASHTABLE* ht, CHROMVARS* chromvars, VARIANT* varlist, REFLIST* reflist) ;

#endif // __FOSMIDBAM_HAIRS_H__