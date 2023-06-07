#ifndef _ESTIMATE_HMM_PARAM_H
#define _ESTIMATE_HMM_PARAM_H

#include <unordered_map>
#include "bamread.h"
#include "readfasta.h"
#include "map"
#include "string"
extern int MIN_MQ;
extern int VERBOSE;

const int BTI[] = {
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 1, 0, 2,  0, 0, 0, 3,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  4, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  // A | C | G | T
        0, 1, 0, 2,  0, 0, 0, 3,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  4, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  // a | c | g | t 
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
};

typedef struct // probabilities are in log10space 
{
        int states; // match, insertion,deletion 
        float** TRS; // transition probabilities between states 
        float** MEM; // probability of A->C, A->G (match state)
        float match; float mismatch; float insertion; float deletion; // emission probs
	// ***probabilities in log10 space ***
        float** lTRS; // transition probabilities between states 
        float** lMEM; // probability of A->C, A->G (match state)
        float lmatch; float lmismatch; float linsertion; float ldeletion; // emission probs
} Align_Params;
extern Align_Params *AP;
Align_Params* init_params();
int estimate_counts_read(struct alignedread* read, REFLIST* reflist,int* emission_counts,int* trans_counts,int* indel_lengths);
void print_error_params(int* emission_counts,int* trans_counts,int* indel_lengths,Align_Params* AP);
int realignment_params(char* bamfile,REFLIST* reflist,char* regions,Align_Params* AP);
#endif