/* rand.h */

#ifndef _LRF_RAND_HEADER_
#define _LRF_RAND_HEADER_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

double gasdev();
long int clock_seed( void );
float *vrnd( float *x, int n, char *dist );
double *dvrndn(double *x,int n);

#endif
