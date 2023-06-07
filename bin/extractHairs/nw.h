#ifndef __NW_H__
#define __NW_H__

extern int VERBOSE;

extern float MATCH;
extern float MISMATCH;
extern float GAP_OPEN;
extern float GAP_EXTEND;
extern float INSERTION_OPEN;
extern float INSERTION_EXTEND;
extern float DELETION_OPEN;
extern float DELETION_EXTEND;

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <stddef.h>
#include<float.h>
#include<ctype.h>

double nw(char* str1, char* str2,int VERBOSE);

#endif // __NW_H__