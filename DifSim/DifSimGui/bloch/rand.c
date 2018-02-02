/**<rand>*************************************************/

#include "rand.h"

/**<dvrndn>*************************************************

    Title: RANDom Normal sequence

  Purpose: Normal random number generator

     Call: double dvrndn(x,n)
	        double x;
		int    n;

      Ref: Numerical Recipes, pg 289

    Notes: Uses unix c-library function drand48()
           for generation of uniformly distributed
	   variates, rather than num rec ran1().

*/
double gasdev()
{
  static int iset = 0;
  static float gset;
  double fac,rsq,v1,v2;

  if (iset == 0) {
    do {
      v1 = 2.0 * drand48() - 1.0;
      v2 = 2.0 * drand48() - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);

    fac = sqrt(-2.0 * log(rsq)/rsq);
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}

/**<clock_seed>*************************************************/

long int clock_seed( void )
{
  long   int iseed;
  static char time_line[80];
  struct tm *times;
  long	 time_sec;

  time (&time_sec);
  times = localtime (&time_sec);

  sprintf (time_line, "%2d%2d%2d\n",
	   times->tm_hour,
	   times->tm_min,
	   times->tm_sec);
  
  iseed = (long) atof(time_line);
  return (iseed);
}

float *vrnd( float *x, int n, char *dist )
{
  int i;
  float *xold = x;

  srand48(clock_seed());

  if(!strcmp(dist,"uniform")) {

    for(i=0; i<n; i++) xold[i] = (float) drand48();

  } else if(!strcmp(dist,"normal")) {

    for(i=0; i<n; i++) xold[i] = (float) gasdev();

  } else {

    printf("vrnd: Unknown distribution requested!\n\a");
    exit(0);

  }
  return(xold);
}

double *dvrndn(double *x,int n)
{
  int i;
  double *xold = x;

  for(i=0; i<n; i++) xold[i] = gasdev();

  return(xold);
}
