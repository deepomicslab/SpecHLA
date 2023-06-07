#ifndef _REALIGN_PAIR_HMM_H
#define _REALIGN_PAIR_HMM_H

#include "estimate_hmm_params.h"
// can do without log-space since sum will not be very small...

// 11/28/2018 
// translated from https://github.com/pjedge/longshot/blob/master/src/realignment.rs

//define states of HMM, emission probability vector and transition probability matrix 
// states = match, insertion, deletion 
// simple dynamic programming over 2-D matrix between v and w (read and haplotype) 

const int ALLOW_END_GAPS = 0; // for anchored alignment...
const int BAND_WIDTH = 20; // default 

//#define logsum(a, b) (((a) > (b)) ? ((a) + log(1.0 + exp(b-a))) : ((b) + log(1.0 + exp( (a) - (b)))))
#define logsum(a, b) (((a) > (b)) ? ((a) + log10(1.0 + pow(10.0, (b) - (a)))) : ((b) + log10(1.0 + pow(10.0, (a) - (b)))))

double logsum1(double a, double b);
double sum_all_alignments_fast(char* v, char* w,Align_Params* params, int min_band_width);
double sum_all_alignments_logspace(char* v, char* w,Align_Params* params, int min_band_width);
#endif