/* vector.c
   Vector routines.
 */
#include "vector.h"

/**<zarg.c>****************************************************

     Title: Float Complex ARGument

   Purpose: Calculate phase of float complex number

      Call: dcomplex zarg(dcomplex z) 

    Output: arg(z)

*/
float zarg(complex z)
{
  float angle ;

  angle = (z.r!=0.0 || z.i!=0.0) ? atan2f(z.i,z.r) : 0.0 ;

  return(angle);
}
/**<dzarg.c>****************************************************

     Title: Double Complex ARGument

   Purpose: Calculate phase of double complex number

      Call: dcomplex dzarg(dcomplex z) 

    Output: arg(z)

*/
double dzarg(dcomplex z)
{
  double angle ;

  angle = (z.dr!=0.0 || z.di!=0.0) ? atan2(z.di,z.dr) : 0.0 ;

  return(angle);
}

/**<vmse>**************************************************

    Title: Float Vector Mean SQuare Error

  Purpose: Calculates mean square error between two
           float vectors x[],y[] of length n.
	   Internal calculation is float precision,
	   sum is returned in float precision
 
     Call: float vmse(float *x, float *y, int n)

   mse = sum from i=0 to n-1 { (x[i]-y[i])**2 }

*/
float vmse(float *x, float *y, int n)
{
  float sum = 0.0;
  float temp ;
  float mse ;
  float dn = n;

  if(n>0)
    {
      do {
	temp = *x++ - *y ;
	y++ ;		/* hack to improve VAX code */
	sum += temp*temp ;
      } while(--n > 0) ;
    }
  mse = sum/dn;
  return(mse);
}
/**<dvmse>**************************************************

    Title: Double Vector Mean SQuare Error

  Purpose: Calculates mean square error between two
           double vectors x[],y[] of length n.
	   Internal calculation is double precision,
	   sum is returned in double precision
 
     Call: double dvmse(double *x, double *y, int n)

   mse = sum from i=0 to n-1 { (x[i]-y[i])**2 }

*/
double dvmse(double *x, double *y, int n)
{
  double sum = 0.0;
  double temp ;
  double mse ;
  double dn = n;

  if(n>0)
    {
      do {
	temp = *x++ - *y ;
	y++ ;		/* hack to improve VAX code */
	sum += temp*temp ;
      } while(--n > 0) ;
    }
  mse = sum/dn;
  return(mse);
}

/**<dvmsq>**************************************************

    Title: Double Vector Mean SQuare

  Purpose: Compute msq of double vector.
 
     Call: double dvmsq( double int *x, int n )

*/
double dvmsq( double *x, int n )
{
  double sum = 0.;
  double dn  = n;
  double msq,tmp;

  if(n>0)
    {
      do {
	tmp  = *x++ ;
	sum += tmp*tmp;
      } while(--n > 0) ;
    }
  msq = sum/dn;
  
  return(msq);
}

/**<dvsqrt>**************************************************

    Title: Double Vector SQuare RooT

  Purpose: This subroutine returns the square root of the double vector
           z[i] = sqrt(x[i])  for i=0,...,N-1

 
     Call: double *dvsqrt(double *z, double *x, int N)

   Output: Returns pointer to z[] array
*/
double *dvsqrt(double *z, double *x, int N)
{
  double *zptr = z ;

  if(N>0)
    {
      do {
	*z++ = sqrt(*x++) ;
      } while(--N > 0) ;
    }
  return(zptr) ;
}
/**<dzvsqrt>**************************************************

    Title: Double Complex Vector SQuare RooT

  Purpose: This subroutine returns the square root of the double complex vector
           z[i] = sqrt(x[i])  for i=0,...,N-1
 
     Call: dcomplex *dzvsqrt( dcomplex *z, dcomplex *x, int N )

   Output: Returns pointer to z[] array

    Notes: Special attention is given to avoiding unnecessary overflow
*/
dcomplex *dzvsqrt( dcomplex *z, dcomplex *x, int N )
{
        register double mag,temp ;
        dcomplex ztemp ;
        register dcomplex *zz = z ;

        while(N-->0) {
                ztemp = *x++ ;
                if((mag = dzabs(ztemp)) == 0) {
                        *zz++ = dzcplx( 0., 0. ) ;

                } else if(ztemp.dr > 0) {
                        temp = sqrt(0.5*(mag + ztemp.dr)) ;
                        *zz++ = dzcplx( temp, 0.5*ztemp.di/temp ) ;

                } else {
                        temp = sqrt(0.5*(mag - ztemp.dr)) ;
                        if(ztemp.di < 0)
                                temp = -temp ;
                        *zz++ = dzcplx( 0.5*ztemp.di/temp , temp ) ;
                }
        }
        return(z) ;
}

/**<dvsum>**************************************************

    Title: Double Vector SUM

  Purpose: Sums elements of vector x[0],...,x[n-1]
 
     Call: double dvsum( double *x, int n )

    Notes: Internal calculation is double precision,
           and the sum is returned in double precision.
*/
double dvsum( double *x, int n )
{
        double sum=0;

        if(n>0)
        {
                do {
                        sum += *x++ ;
                } while(--n > 0) ;
        }
        return(sum);
}

/**<dvmean>**************************************************

    Title: Double Vector MEAN

  Purpose: Compute mean of double vector.
 
     Call: double dvmean( double *x, int n )

*/
double dvmean( double *x, int n )
{
  double sum, mean;
  double dn = n;

  sum  = dvsum(x,n);
  mean = ((double) sum)/dn;

  return(mean);
}

/**<vvar>**************************************************

    Title: Float Vector VARiance

  Purpose: Compute sample variance of float vector.
 
     Call: float vvar( float *x, int n)

*/
float vvar( float *x, int n )
{
       float dn = n-1 ;
       float mean,var,temp  ;
       float sum=0.;

       mean = vmean(x,n);

       if(n>0)
	 {
	   do {
	     temp = *x++ - mean ;
	     sum += temp*temp ;
	   } while(--n > 0) ;
	 }
  
       var = sum / dn;
       return(var);
}
/**<dvvar>**************************************************

    Title: Double Vector VARiance

  Purpose: Compute sample variance of double vector.
 
     Call: double dvvar( double *x, int n)

*/
double dvvar( double *x, int n )
{
       double dn = n-1 ;
       double mean,var,temp  ;
       double sum=0.;

       mean = dvmean(x,n);

       if(n>0)
	 {
	   do {
	     temp = *x++ - mean ;
	     sum += temp*temp ;
	   } while(--n > 0) ;
	 }
  
       var = sum / dn;
       return(var);
}

/**<vstd>**************************************************

    Title: Float Vector Standard Deviation

  Purpose: Compute sample standard deviation of float vector.
 
     Call: float vstd( float *x, int n)

*/
float vstd( float *x, int n )
{
  return sqrtf(vvar(x,n)) ;
}

/**<dvstd>**************************************************

    Title: Double Vector Standard Deviation

  Purpose: Compute sample standard deviation of double vector.
 
     Call: double dvstd( double *x, int n)

*/
double dvstd( double *x, int n )
{
  return sqrt(dvvar(x,n)) ;
}

/**<dvneg.c>****************************************************

     Title: Double Vector NEGate

   Purpose: Computes z[i] = -x[i]  for i=0,...,N-1

      Call: double *dvneg(double *z, double *x, int n)

     Input: x = input vector
	    n = length(x)

    Output: Returns pointer to z[] array
*/
double *dvneg(double *z, double *x, int N)
{
  double *zptr = z ;

  if(N>0)
    {
      do {
	*z++ = -(*x++) ;
      } while( --N > 0 ) ;
    }
  return(zptr) ;

}

/**<vindx>**************************************************

     Title: Float Vector INDeXing

   Purpose: Creating index table for sorting arrays
            into ascending order

      Call: void vindx(float *x,int *indx,int n)

     Usage: To sort an array x, then, do:
        
            vindx(x,indx,n);

	    Then print elements of x in ascending order 

	    for(i=0; i< n; i++) {
	       printf("x[indx[%i]] = %f\n",i,x[indx[i]]);
	    }
*/
int fbc( const void *ap , const void *bp )   /* for qsort */
{
   float *a=(float *)ap,  *b=(float *)bp ;
   if( *a == *b ) return 0  ;
   if( *a <  *b ) return -1 ;
                  return  1 ;
}

int fcompare(const void *i, const void *j)
{
  float **ii = (float **)i;
  float **jj = (float **)j;
  if (**ii < **jj)
    return -1;
  else
    return 1;
}

void vindx(float *x,long int *indx,int n)
{
  int i;

  for (i=0; i<n; i++) indx[i] = (long int)&x[i]; 
  qsort((void *)indx,n,sizeof(float *),fcompare);
  //qsort((void *)indx,n,sizeof(float *),fbc); //DOESN'T WORK!!! SEE tstopt.c
  for (i=0; i<n; i++) indx[i] = ((float *)indx[i]-x);

}

/**<dvindx>**************************************************

     Title: Double Vector INDeXing

   Purpose: Creating index table for sorting arrays
            into ascending order

      Call: void dvindx(double *x,int *indx,int n)

     Usage: To sort an array x, then, do:
        
            dvindx(x,indx,n);

	    Then print elements of x in ascending order 

	    for(i=0; i< n; i++) {
	       printf("x[indx[%i]] = %f\n",i,x[indx[i]]);
	    }
*/
int dbc( const void *ap , const void *bp )   /* for qsort */
{
   double *a=(double *)ap,  *b=(double *)bp ;
   if( *a == *b ) return 0  ;
   if( *a <  *b ) return -1 ;
                  return  1 ;
}

int dcompare(const void *i, const void *j)
{
  double **ii = (double **)i;
  double **jj = (double **)j;
  if (**ii < **jj)
    return -1;
  else
    return 1;
}

void dvindx(double *x, long int *indx, int n)
{
  int i;

  for (i=0; i<n; i++) indx[i] = (long int)&x[i]; 
  qsort((void *)indx,n,sizeof(double *),dcompare);
  //qsort((void *)indx,n,sizeof(double *),dbc); //DOESN'T WORK!!! SEE tstopt.c
  for (i=0; i<n; i++) indx[i] = ((double *)indx[i]-x);

}

/**<vmedian>**************************************************

    Title: Float Vector MEDIAN

  Purpose: Compute median of float vector
 
     Call: float vmedian( float *x, int n )

*/
float vmedian( float *x, int n )
{
  int nmid = (n-1)/2 ;
  float median ;

  long int *iv = (long int *) calloc(n,sizeof(long int)) ;

  vindx(x,iv,n) ;

  if(ODD(n)) {
    median = x[iv[nmid]] ;
  } else {
    median = ( x[iv[nmid]] + x[iv[nmid+1]] ) / 2.;
  }

  free(iv) ;
  return median;
}

/**<dvmedian>**************************************************

    Title: Double Vector MEDIAN

  Purpose: Compute median of double vector
 
     Call: double dvmedian( double *x, int n )

*/
double dvmedian( double *x, int n )
{
  int nmid = (n-1)/2 ;
  double median ;

  long int *iv = (long int *) calloc(n,sizeof(long int)) ;

  dvindx(x,iv,n) ;

  if(ODD(n)) {
    median = x[iv[nmid]] ;
  } else {
    median = ( x[iv[nmid]] + x[iv[nmid+1]] ) / 2.;
  }

  free(iv) ;
  return median;
}

/**<ivany>**************************************************

    Title: Integer Vector ANY

  Purpose: Check to see if any elements of an integer vector "x"
           matches a given value "ival"
 
     Call: int ivany( int *x, int n )

*/
bool ivany( int *x, int ival, int n )
{
  bool ans=0;  

  if(ivmatch(ival,x,n)>0) ans=1;

  return ans;
}

/**<dvany>**************************************************

    Title: Double Vector ANY

  Purpose: Check to see if any elements of an double vector "x"
           matches a given value "dval"
 
     Call: double dvany( double *x, int n )

*/
bool dvany( double *x, double dval, int n )
{
  bool ans=0;  

  if(dvmatch(dval,x,n)>0) ans=1;

  return ans;
}

/**<dvnan>**************************************************

    Title: Double Vector NAN

  Purpose: Check to see if any elements of an double vector "x"
           are "nan"
 
     Call: double dvnan( double *x, int n )

*/
int dvnan( double *x, int n )
{
  int i, count=0;  

  for (i=0; i<n; i++) {
    if isnan(x[i]) count ++ ;
  }

  return count;
}

/**<ivmatch>**************************************************

    Title: Integer Vector MATCH

  Purpose: Count number of elements of an integer vector "x"
           that match a given value "ival"
 
     Call: int ivmatch( int *x, int n )

*/
int ivmatch( int ival, int *x, int n )
{
  int count=0;

  if(n>0)
    {
      do {
	count += (*x++)==ival ;
      } while(--n > 0) ;
    }
  return(count);
}

/**<dvmatch>**************************************************

    Title: Double Vector MATCH

  Purpose: Count number of elements of an double vector "x"
           that match a given value "dval"
 
     Call: double ivmatch( double *x, int n )

*/
int dvmatch( double ival, double *x, int n )
{
  int count=0;

  if(n>0)
    {
      do {
	count += (*x++)==ival ;
      } while(--n > 0) ;
    }
  return(count);
}

/**<ivnonzero>**************************************************

    Title: Integer Vector NON ZERO 

  Purpose: Count non zero elements of an integer vector
 
     Call: int ivnonzero( int *x, int n )

*/
int ivnonzero( int *x, int n )
{
  int count=0;

  if(n>0)
    {
      do {
	count += (*x++)!=0 ;
      } while(--n > 0) ;
    }
  return(count);
}
/* Next two don't work in this directory, but do in grd - fix later */
#if 0
/**<ivfind>**************************************************

    Title: Integer Vector FIND non zero elements

  Purpose: Return array of indeces of non zero elements of an integer vector
 
     Call: int ivfind( &idx, int *x, int n )

   Return: returns length nidx of index array idx[nidx]
*/
int ivfind( int **idx, int *x, int n )
{
  int i,count=0;

  int *pidx = (int *) calloc((n),sizeof(int));

  if(n>0)
    {
      for(i=0; i<n; i++) if(x[i]!=0) pidx[count++] = i;
    }

  pidx = realloc(pidx,count*sizeof(int));
  *idx = pidx;

  return(count);
}
/**<bvfind>**************************************************

    Title: Byte Vector FIND non zero elements

  Purpose: Return array of indeces of non zero elements of a byte vector
 
     Call: int ivfind( &idx, byte *x, int n )

   Return: returns length nidx of index array idx[nidx]
*/
int bvfind( int **idx, byte *x, int n )
{
  int i,count=0;

  int *pidx = (int *) calloc((n),sizeof(int));

  if(n>0)
    {
      for(i=0; i<n; i++) if(x[i]!=0) pidx[count++] = i;
    }

  pidx = realloc(pidx,count*sizeof(int));
  *idx = pidx;

  return(count);
}
#endif

/**<vsum>**************************************************

    Title: Float Vector SUM

  Purpose: Sums elements of vector x[0],...,x[n-1]
 
     Call: float vsum( float *x, int n )

    Notes: Internal calculation is float precision,
           and the sum is returned in float precision.
*/
float vsum( float *x, int n )
{
        float sum=0;

        if(n>0)
        {
                do {
                        sum += *x++ ;
                } while(--n > 0) ;
        }
        return(sum);
}
/**<vmean>**************************************************

    Title: Float Vector MEAN

  Purpose: Compute mean of float vector.
 
     Call: float vmean( float *x, int n )

*/
float vmean( float *x, int n )
{
  float mean;
  float dn = n;

  mean = vsum(x,n)/dn;

  return(mean);
}

/**<zmuli.c>****************************************************

     Title: Complex MULTiplication by i 

   Purpose: Multiply complex number by sqrt(-1)

      Call: complex zmuli(complex a) 

     Input: a = complex arg

    Output: a

*/
complex zmuli(complex a) 
{
	static complex c;
	c.r = -a.i;
	c.i =  a.r;
	return(c);
}

/**<zcnjg.c>****************************************************

     Title: Complex CoNJuGate

   Purpose:  returns complex conjugate of complex arg  

      Call: complex zcnjg(complex a)

     Input: returns a*

*/
complex zcnjg(complex a)
{
	static complex c;
	c.r =  a.r;
	c.i = -a.i;
	return(c);
}

/**<dzcnjg.c>****************************************************

     Title: Double Complex CoNJuGate

   Purpose:  returns dcomplex conjugate of dcomplex arg  

      Call: dcomplex dzcnjg(dcomplex a)

     Input: returns a*

*/
dcomplex dzcnjg(dcomplex a)
{
	static dcomplex c;
	c.dr =  a.dr;
	c.di = -a.di;
	return(c);
}

/**<zadd.c>****************************************************

     Title: Complex ADDition

   Purpose: add two complex numbers

      Call: complex zadd(complex a, complex b)

     Input: returns complex sum

*/
complex zadd(complex a, complex b)
{
	static complex c;
	c.r = a.r + b.r;
	c.i = a.i + b.i;
	return(c);
}
/**<dzadd.c>****************************************************

     Title: Dcomplex ADDition

   Purpose: add two dcomplex numbers

      Call: dcomplex zadd(dcomplex a, dcomplex b)

     Input: returns dcomplex sum

*/
dcomplex dzadd(dcomplex a, dcomplex b)
{
	static dcomplex c;
	c.dr = a.dr + b.dr;
	c.di = a.di + b.di;
	return(c);
}
/**<zmean.c>****************************************************

     Title: Complex MEAN

   Purpose: Mean of two complex numbers

      Call: complex zmean(complex a, complex b)

     Input: returns complex mean

*/
complex zmean(complex a, complex b)
{
	complex sum, mean;

	sum = zadd(a,b);
	mean.r = sum.r/2.;
	mean.i = sum.i/2.;
	
	return(mean);
}

/**<zsub.c>****************************************************

     Title: Complex SUBtraction

   Purpose: Complex difference

      Call: complex zsub(complex a,complex b) 

    Output: a-b

*/
complex zsub(complex a,complex b)
{
	static complex c;
	c.r = a.r - b.r;
	c.i = a.i - b.i;
	return(c);
}

/**<dzsub.c>****************************************************

     Title: Double Complex SUBtraction

   Purpose: Double Complex difference

      Call: dcomplex dzsub(dcomplex a,dcomplex b) 

    Output: a-b

*/
dcomplex dzsub(dcomplex a,dcomplex b)
{
	static dcomplex c;
	c.dr = a.dr - b.dr;
	c.di = a.di - b.di;
	return(c);
}

/**<zmul.c>****************************************************

     Title: Complex MULtiplication 

   Purpose: Multiply two complex numbers

      Call: complex zmul(complex a,complex b)  

    Output: a*b

*/
complex zmul(complex a,complex b)  
{
	static complex c;
        c.r = a.r*b.r - a.i*b.i;
        c.i = a.r*b.i + a.i*b.r;
	return(c);
}
/**<zmulc.c>****************************************************

     Title: Complex Conjugate MULtiplication 

   Purpose: Multiply complex number by complex conjugate 
            of another complex number

      Call: complex zmulc(complex a,complex b)  

    Output: a*conj(b)

*/
complex zmulc(complex a, complex b)  
{
	static complex c;
	c.r = a.r*b.r + a.i*b.i;
	c.i = a.r*b.i - a.i*b.r;
	return(c);
}
/**<zscal.c>****************************************************

     Title: Complex SCALe

   Purpose: Scales a complex number

      Call: complex zscal( double a , complex z )

    Output: a*z

*/
complex zscal( double a , complex z )
{
	static complex c ;
	c.r = a*z.r ;
	c.i = a*z.i ;
	return(c) ;
}
/**<dzscal.c>****************************************************

     Title: Double Complex SCALe

   Purpose: Scales a dcomplex number

      Call: dcomplex dzscal( double a , dcomplex z )

    Output: a*z

*/
dcomplex dzscal( double a , dcomplex z )
{
	static dcomplex c ;
	c.dr = a*z.dr ;
	c.di = a*z.di ;
	return(c) ;
}
/**<factorial.c>****************************************************

     Title: FACTORIAL

   Purpose: Computer N!

      Call: long int factorial(int N)

     Input: N

    Output: N! = N*(N-1)* ... *1

*/
long int factorial(int N)
{
  long int ii = 1;

  if(N==0) return(1);

  do { ii *= N; } while(N-->1) ;

  return(ii) ;
}

/**<vdif.c>****************************************************

     Title: float Vector DIFference

   Purpose: Compute difference

      Call: long int factorial(int N)

     Input: x[n]

    Output: dx[i] = x[i+1] - x[i] ; i=0,..,n-1

     Notes: both x and dx are length n for ease of call
            but only the first n-1 elements of dx are used

*/
float *vdif( float *dx, float *x, int n)
{
  int i,n1 = n - 1;
  
  for(i=0; i<n1; i++) dx[i] = x[i+1] - x[i] ;

  dx[n1] = 0;			/* filler */

  return (dx) ;
}

/**<dvdif.c>****************************************************

     Title: Double Vector DIFference

   Purpose: Compute difference

      Call: long int factorial(int N)

     Input: x[n]

    Output: dx[i] = x[i+1] - x[i] ; i=0,..,n-1

     Notes: both x and dx are length n for ease of call
            but only the first n-1 elements of dx are used

*/
double *dvdif( double *dx, double *x, int n)
{
  int i,n1 = n - 1;
  
  for(i=0; i<n1; i++) dx[i] = x[i+1] - x[i] ;

  dx[n1] = 0;			/* filler */

  return (dx) ;
}

/**<vgrad.c>****************************************************

     Title: Float Vector GRADient

   Purpose: Compute gradient

      Call: float *vgrad( float *gx, float *x, int n)

     Input: x[n]

    Output: grad(x)

*/
float *vgrad( float *gx, float *x, int n)
{
  int i,n1 = n - 1;
  float *h = (float *) calloc(n,sizeof(float)) ;
  
  for(i=0; i<n; i++) h[i] = i;	/* useless now, but maybe not later?
				 (see matlab's grad) */

  if(n>0) {
    gx[0]  = ( x[1]  - x[0]   ) / ( h[1] - h[0]   );
    gx[n1] = ( x[n1] - x[n1-1]) / (h[n1] - h[n1-1]);
  }

  if(n>1) {
    for(i=1; i<n1; i++) gx[i] = (x[i+1] - x[i-1])/(h[i+1] - h[i-1]);
  }

  free(h);
  return (gx) ;
}

/**<dvgrad.c>****************************************************

     Title: Double Vector GRADient

   Purpose: Compute gradient

      Call: double *dvgrad( double *gx, double *x, int n)

     Input: x[n]

    Output: grad(x)

*/
double *dvgrad( double *gx, double *x, int n)
{
  int i,n1 = n - 1;
  double *h = (double *) calloc(n,sizeof(double)) ;
  
  for(i=0; i<n; i++) h[i] = i;	/* useless now, but maybe not later?
				 (see matlab's grad) */

  if(n>0) {
    gx[0]  = ( x[1]  - x[0]   ) / ( h[1] - h[0]   );
    gx[n1] = ( x[n1] - x[n1-1]) / (h[n1] - h[n1-1]);
  }

  if(n>1) {
    for(i=1; i<n1; i++) gx[i] = (x[i+1] - x[i-1])/(h[i+1] - h[i-1]);
  }

  free(h);
  return (gx) ;
}

/**<dzvgrad.c>****************************************************

     Title: Double Complex Vector GRADient

   Purpose: Compute gradient

      Call: dcomplex *dzvgrad( dcomplex *gx, dcomplex *x, int n)

     Input: x[n]

    Output: grad(x)

*/
dcomplex *dzvgrad( dcomplex *gx, dcomplex *x, int n)
{
  int i,n1 = n - 1;
  double *h = (double *) calloc(n,sizeof(double)) ;
  
  for(i=0; i<n; i++) h[i] = i;	/* useless now, but maybe not later?
				 (see matlab's grad) */

  if(n>0) {
    gx[0].dr  = ( x[1].dr  - x[0].dr   ) / ( h[1] - h[0]   );
    gx[n1].dr = ( x[n1].dr - x[n1-1].dr) / (h[n1] - h[n1-1]);

    gx[0].di  = ( x[1].di  - x[0].di   ) / ( h[1] - h[0]   );
    gx[n1].di = ( x[n1].di - x[n1-1].di) / (h[n1] - h[n1-1]);
  }

  if(n>1) {
    for(i=1; i<n1; i++) {
      gx[i].dr = (x[i+1].dr - x[i-1].dr)/(h[i+1] - h[i-1]);
      gx[i].di = (x[i+1].di - x[i-1].di)/(h[i+1] - h[i-1]);
    }
  }

  free(h);
  return (gx) ;
}

/**<unwrap.c>****************************************************

     Title: UNWRAP phase angles in-place 

   Purpose: 

      Call: float *unwrap( float *p, int n, float *maskin )

     Input: 

    Output: 

     Notes: WARNING:  Do NOT use fmodf here!
*/
float *unwrap( float *p, int n, float *maskin )
{
  int i;
  float cutoff = M_PI;		/* could be input ... */
  float dd;
  float *dp = (float *) calloc(n,sizeof(float)) ;
  float *mask = maskin;

  if(maskin==NULL) {
    mask = (float *) calloc(n,sizeof(float)) ;
    for(i=0; i<n; i++) mask[i] = 1.0;
  }

  vdif(dp,p,n);

  for(i=0; i<n; i++) {
    dd = FMOD( (dp[i] + M_PI) , TWOPI) - M_PI ;
    if ( (dd==-M_PI) && (dp[i]>0.) ) dd = M_PI;
    if( fabs(dp[i])<cutoff || !mask[i] ) dp[i] = 0.;
    else dp[i] = dd - dp[i];
  }

  dd = 0.;
  for(i=0; i<(n-1); i++) {
    dd += dp[i];
    p[i+1] += dd;
  }

  free(dp);
  if(maskin==NULL) free(mask);
  return(p) ;
}

/**<vlinsp.c>************************************************************

     Title: Float Vector LINearly SPaced array

   Purpose: Create linearly spaced array of float numbers

      Call: float *vlinsp(float *x,float f1,float f2,int n)

     Input:  x = float x[n] 
            f1 = first point 
            f2 = last point
	     n = number of points

    Output: float *x = [f1,..,f2]

    Author: Lawrence R. Frank

   History: C version of a MATLAB routine
*/
float *vlinsp( float *x, float f1, float f2, int n )
{
  int    i;
  float fn = n;

  for(i=0; i<n-1; i++) x[i] = f1 + i * (f2 - f1)/(fn - 1.);

  x[n-1] = f2;
  
  return(x);
}
/**<dvlinsp.c>************************************************************

     Title: Double Vector LINearly SPaced array

   Purpose: Create linearly spaced array of double numbers

      Call: double *dvlinsp(double *x,double f1,double f2,int n)

     Input:  x = double x[n] 
            f1 = first point 
            f2 = last point
	     n = number of points

    Output: double *x = [f1,..,f2]

    Author: Lawrence R. Frank

   History: C version of a MATLAB routine
*/
double *dvlinsp( double *x, double f1, double f2, int n )
{
  int    i;
  double fn = n;

  for(i=0; i<n-1; i++) x[i] = f1 + i * (f2 - f1)/(fn - 1.);

  x[n-1] = f2;
  
  return(x);
}
/**<svlinsp.c>************************************************************

     Title: Short Vector LINearly SPaced array

   Purpose: Create linearly spaced array of short numbers

      Call: short *svlinsp(short *x,short f1,short f2,int n)

     Input:  x = short x[n] 
            f1 = first point 
            f2 = last point
	     n = number of points

    Output: short *x = [f1,..,f2]

    Author: Lawrence R. Frank

   History: C version of a MATLAB routine
*/
short *svlinsp( short *x, short f1, short f2, int n )
{
  int    i;
  short fn = n;

  for(i=0; i<n-1; i++) x[i] = f1 + i * (int) ((double)(f2 - f1)/((double)(fn - 1.)));

  x[n-1] = f2;
  
  return(x);
}

/**<vtrn2.c>************************************************************

     Title: Float Vector TRaNspose in 2 dimensions

   Purpose: 

     Input: 

*/
float *vtrn2( float *z, int nx, int ny )
{
  int i,j ;
  int nxy = nx*ny;
  float *zt = (float *) malloc( sizeof(float) * (nxy) ) ;

  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      zt[j+i*ny] = z[i+j*nx] ;
    }
  }
  vmove(z,zt,nxy);
  free(zt);
  return (z) ;
}
/**<zvtrn2.c>************************************************************

     Title: Complex Vector TRaNspose in 2 dimensions

   Purpose: 

     Input: 

*/
complex *zvtrn2( complex *z, int nx, int ny )
{
  int i,j ;
  int nxy = nx*ny;
  complex *zt = (complex *) malloc( sizeof(complex) * (nxy) ) ;

  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      zt[j+i*ny].r = z[i+j*nx].r ;
      zt[j+i*ny].i = z[i+j*nx].i ;
    }
  }
  zvmove(z,zt,nxy);
  free(zt);
  return (z) ;
}
/**<zcplx.c>****************************************************

     Title: ComPLeX

   Purpose: Create complex number from real and imag components

      Call: complex zcplx(float r, float i)	

     Input: r = real component
            i = imag component

    Output: z = {c.r,c.i}

*/
complex zcplx(double r, double i)
{
	static complex c;
	c.r = r;
	c.i = i;
	return(c);
}
/**<dzcplx.c>****************************************************

     Title: Double ComPLeX

   Purpose: Create complex number from real and imag components

      Call: complex dzcplx(double r, double i)	

     Input: r = real component
            i = imag component

    Output: dz = {c.r,c.i}

*/
dcomplex dzcplx(double r,double i)
{
	static dcomplex c;
	c.dr = r;
	c.di = i;
	return(c);
}
/**<zexp.c>****************************************************

     Title: Float Complex EXPonential

   Purpose: Computes the dcomplex exponential of a dcomplex number

      Call: dcomplex dzexp(dcomplex z)

     Input: z = double complex number

    Output: exp(z) 

*/
complex zexp(complex z)
{
	float rexp ;

	rexp = exp(z.r) ;
	return(zcplx( rexp*cos(z.i) , rexp*sin(z.i) )) ;
}
/**<dzexp.c>****************************************************

     Title: Double Complex EXPonential

   Purpose: Computes the dcomplex exponential of a dcomplex number

      Call: dcomplex dzexp(dcomplex z)

     Input: z = double complex number

    Output: exp(z) 

*/
dcomplex dzexp(dcomplex dz)
{
	double rexp ;

	rexp = exp(dz.dr) ;
	return(dzcplx( rexp*cos(dz.di) , rexp*sin(dz.di) )) ;
}

/**<zexpI.c>****************************************************

     Title: Double Complex Exonential exp(i*theta)

   Purpose: Compute exp(i*theta)

      Call: dcomplex dzexpI(double theta)

     Input: double theta 

    Output: exp(i*theta)

*/
complex zexpI(float theta)
{
  return(zexp(zcplx(0.,theta))) ;
}

/**<dzexpI.c>****************************************************

     Title: Double Complex Exonential exp(i*theta)

   Purpose: Compute exp(i*theta)

      Call: dcomplex dzexpI(double theta)

     Input: double theta 

    Output: exp(i*theta)

*/
dcomplex dzexpI(double theta)
{
  return(dzexp(dzcplx(0.,theta))) ;
}

/**<dzmul.c>****************************************************

     Title: Double Complex MULtiplication

   Purpose: Computes dcomplex product of 2 dcomplex args 

      Call: dcomplex dzmul(dcomplex a, dcomplex b)

     Input: dcomplex a, b

    Output: a*b

*/
dcomplex dzmul(dcomplex a, dcomplex b)
{
	static dcomplex c;
	c.dr = a.dr*b.dr - a.di*b.di;
	c.di = a.dr*b.di + a.di*b.dr;
	return(c);
}

/**<zrot.c>****************************************************

     Title: Double Complex ROTation

   Purpose: Rotates dcomplex number z by theta

      Call: dcomplex dzmul(dcomplex a, dcomplex b)

     Input: dcomplex a, b

    Output: a*b

*/
complex zrot(complex z, float theta)
{
	return(zmul(z,zexpI(theta)));
}

/**<dzrot.c>****************************************************

     Title: Double Complex ROTation

   Purpose: Rotates dcomplex number z by theta

      Call: dcomplex dzmul(dcomplex a, dcomplex b)

     Input: dcomplex a, b

    Output: a*b

*/
dcomplex dzrot(dcomplex dz, double theta)
{
	return(dzmul(dz,dzexpI(theta)));
}

/**<dzvrot.c>****************************************************

     Title: Double Complex Vector ROTation

   Purpose: Rotates dcomplex vector by theta

      Call: dcomplex *dzvrot(dcomplex *zrot, dcomplex *z, double *theta, int sign, int N)

     Input: z[N], theta[N], N

    Output: z[N]*exp(i*theta*sign)
*/
dcomplex *dzvrot(dcomplex *zrot, dcomplex *z, double *theta, int sign, int N)
{
  dcomplex *zptr = zrot ;

  while(N-->0) {
    *zptr++ = dzrot(*z,*theta * sign) ;
    theta++ ;
    z++ ;
  }

  return(zrot) ;
}

/**<vrotxy.c>****************************************************

     Title: Float Vector ROTation of X and Y components

   Purpose: Rotates complex vector that is expressed as x and y 

      Call: void vrotxy(float *x, float *y, float *phi, int n) 

     Input: x[n] - real channel
            y[n] - imag channel

    Output: x and y components of (x[n] + i y[n])*exp(i*phi)
*/
void vrotxy(float *x, float *y, float phi, int n) 
{
  int   i;
  float c = cos(phi);  
  float s = sin(phi);
  for (i=0; i<n; i++) {
    x[i] = c*x[i] - s*y[i];
    y[i] = s*x[i] + c*y[i];
  }
  return ;
}

/**<zabs.c>****************************************************

     Title: Complex ABSolute value

   Purpose:  returns double prec magnitude  
             of complex argument 
	     spec attn to avoid unnec overflow 

      Call: double zabs(complex a)

     Input: a = complex arg

    Output: abs(a)

*/
double zabs(complex a)
{                          
	double absr, absi, ratio;

	if( (absr=fabs(a.r)) <= (absi=fabs(a.i)) && absi>0.0 ) {
		ratio = absr/absi;
		return(absi * sqrt(1 + ratio*ratio) ); 
	}
	else if (absr>0.0) {
		ratio = absi/absr;
		return(absr * sqrt(1 + ratio*ratio) );
	}
	else return((double)0.0);
}
/**<dzabs.c>****************************************************

     Title: Double Complex ABSolute value

   Purpose:  returns double prec magnitude  
             of double complex argument 
	     spec attn to avoid unnec overflow 

      Call: double dzabs(dcomplex a)

     Input: a = dcomplex arg

    Output: abs(a)

*/
double dzabs(dcomplex a)
{
	double absr, absi, ratio;

	if( (absr=fabs(a.dr)) <= (absi=fabs(a.di)) && absi>0.0 ) {
		ratio = absr/absi;
		return(absi * sqrt(1 + ratio*ratio) ); 
	}
	else if (absr>0.0) {
		ratio = absi/absr;
		return(absr * sqrt(1 + ratio*ratio) );
	}
	else return((double)0.0);
}

/**<zpolr.c>****************************************************

     Title: Complex POLaR conversion

   Purpose:  This routine converts a complex number in rectangular coordinates 
               with z = x + j*y
             to polar coordinates, returning:
		z = radius + j*phase
*/
complex zpolr(complex z)
{
	return(zcplx(zabs(z),zarg(z))) ;
}
/**<dzpolr.c>****************************************************

     Title: Complex POLaR conversion

   Purpose:  This routine converts a complex number in rectangular coordinates 
               with z = x + j*y
             to polar coordinates, returning:
		z = radius + j*phase
*/
dcomplex dzpolr(dcomplex z)
{
	return(dzcplx(dzabs(z),dzarg(z))) ;
}

/**<zvpolr.c>****************************************************

     Title: Complex Vector POLaR conversion

   Purpose: Converts a vector of complex numbers in 
            rectangular coordinates to polar coordinates, 
	    returning:
	    z[i] = radius(x[i]) + j*phase(x[i])

      Call: complex *zvpolr(complex *z, complex *x,int N)

     Input: 

    Output: 

*/
complex *zvpolr(complex *z, complex *x,int N)
{
  complex *zptr = z ;

  while(N-->0) {
    *zptr++ = zcplx(zabs(*x), zarg(*x)) ;
    x++ ;
  }
  return(z) ;
}

/**<dzvpolr.c>****************************************************

     Title: Double complex Vector POLaR conversion

   Purpose: This routine converts a vector of dcomplex numbers 
            in rectangular coordinates to polar coordinates, returning:
		z[i] = radius(x[i]) + j*phase(x[i])

 History:     musicus - original version
          24oct95 lrf - handles abs(x) = 0 correctly

      Call: dcomplex *dzvpolr(dcomplex *z, dcomplex *x, int N)

*/
dcomplex *dzvpolr(dcomplex *z, dcomplex *x, int N)
{
  double angle;
  dcomplex *zptr = z ;

  while(N-->0) {
    *zptr++ = dzcplx(dzabs(*x), dzarg(*x)) ;
    x++ ;
  }
  return(z) ;
}

/**<dzvarg.c>****************************************************

     Title: Double complex Vector ARGument

   Purpose: This routine computes the phase of a vector of dcomplex numbers 
		arg[i] = angel(z[i])

      Call: double *dzvarg(double *arg, dcomplex *z, int N)

 History: 09oct06 lrf - handles abs(x) = 0 correctly

*/
double *dzvarg(double *arg, dcomplex *z, int N)
{
  double *aptr = arg ;

  while(N-->0) {
    *aptr++ = dzarg(*z) ;
    z++ ;
  }
  return(arg) ;
}

/**<zrect.c>****************************************************

     Title: Complex RECTangular conversion

   Purpose: This routines converts a complex number in polar coordinates, 
		z = radius + j*phase
	    into rectangular coordinates
		z = x + j*y

      Call: complex zrect(complex z)

*/
complex zrect(complex z)
{
	return(zcplx(z.r*cos(z.i),z.r*sin(z.i))) ;
}

/**<dzrect.c>****************************************************

     Title: Double Complex RECTangular conversion

   Purpose: This routines converts a complex number in polar coordinates, 
		z = radius + j*phase
	    into rectangular coordinates
		z = x + j*y

      Call: dcomplex dzrect(dcomplex z)

*/
dcomplex dzrect(dcomplex z)
{
	return(dzcplx(z.dr*cos(z.di),z.dr*sin(z.di))) ;
}

/**<zvrect.c>****************************************************

     Title: Complex Vector RECTangular conversion

   Purpose: This routines converts a complex vector in polar coordinates, 
               with x[i] = radius + j*phase
	    into rectangular coordinates
	       z[i] = x + j*y

      Call: complex *zvrect( complex *z, complex *x, int N )

     Input: 

    Output: 
*/
complex *zvrect( complex *z, complex *x, int N )
{
	complex *zptr = z ;

	while(N-->0) {
		*zptr++ = zcplx(x->r * cos(x->i), x->r * sin(x->i)) ;
		x++ ;
	}
	return(z) ;
}

/**<dzvrect.c>****************************************************

     Title: Double Complex Vector RECTangular conversion

   Purpose: This routines converts a double complex vector in polar coordinates, 
               with x[i] = radius + j*phase
	    into rectangular coordinates
	       z[i] = x + j*y

      Call: dcomplex *dzvrect( dcomplex *z, dcomplex *x, int N )

     Input: 

    Output: 
*/
dcomplex *dzvrect(dcomplex *z, dcomplex *x, int N)
{
	dcomplex *zptr = z ;

	while(N-->0) {
		*zptr++ = dzcplx(x->dr * cos(x->di), x->dr * sin(x->di)) ;
		x++ ;
	}
	return(z) ;
}

/**<zvmove.c>*******************************************************

     Title: Complex MOVE

      complex *zvmove(z,x,n)	Moves complex vector x[] into z[] quickly
        complex *z,*x ;		  z[i] = x[i] for i=0,...,n-1 if z<=x
 	int n ;				      for i=n-1,...,0 if z>x
 				Returns pointer to z[]
 
  Works correctly regardless of how x[] and z[] might overlap
*/
complex *zvmove(complex *z, complex *x, int N)
{
	complex *zold = z;
	float *px = (float*) x ;
	float *pz = (float*) z ;

	if(N>0)
	{
		if(pz<=px)
		{
			do {
				*pz++ = *px++ ;
				*pz++ = *px++ ;
			} while( --N > 0) ;
		}
		else
		{
			pz += 2*N ;
			px += 2*N ;
			do {
				*--pz = *--px ;
				*--pz = *--px ;
			} while( --N > 0) ;
		}
	}
	return(zold) ;
}
/*      dcomplex *dzvmove(z,x,n) Moves dcomplex vector x[] into z[] quickly
 *      dcomplex *z,*x ;	  z[i] = x[i] for i=0,...,n-1 if z<=x
 *	int n ;				      for i=n-1,...,0 if z>x
 *				 Returns pointer to z[]
 *
 * Works correctly regardless of how x[] and z[] might overlap
*/
dcomplex *dzvmove(dcomplex *z,dcomplex *x, int N)
{
	dcomplex *zold = z;
	double *px = (double*) x ;
	double *pz = (double*) z ;

	if(N>0)
	{
		if(pz<=px)
		{
			do {
				*pz++ = *px++ ;
				*pz++ = *px++ ;
			} while( --N > 0) ;
		}
		else
		{
			pz += 2*N ;
			px += 2*N ;
			do {
				*--pz = *--px ;
				*--pz = *--px ;
			} while( --N > 0) ;
		}
	}
	return(zold) ;
}

/**<zvcnjg.c>*******************************************************

     Title: Complex CONJuGation

     complex *zvcnjg(z,x,N)		Sets z[i] = conjugate of x[i]
     complex *z,*x ;			Handles in-place calculation specially
     int N ;				  for extra speed
					Returns pointer to z[] array
*/
complex *zvcnjg(complex *z, complex *x, int N)
{
	float *rz = (float*) z ;
	float *rx = (float*) x ;

	if(N>0)
	{
		if(z==x)
		{
			rx++ ;
			do {
				*rx = - *rx ;
				rx += 2 ;
			} while(--N > 0) ;
		}
		else
		{
			do {
				*rz++ = *rx++ ;
				*rz++ = - *rx++ ;
			} while(--N > 0) ;
		}
	}
	return(z) ;
}
/**<dvsubmat.c>*******************************************************

     Title: Double Vector SUBMATrix extraction

   Purpose: Extract submatrix b from larger matrix a.
            b = a(rstart:rend,cstart:cend); 

      Call: dvsubmat(b, a, ma, na, rstart, rend, cstart, cend)
              double *b;       Result matrix (rend-rstart+1)-by-(cend-cstart+1)
	      double *a;       Original matrix ma-by-na 
	      int ma;	       Row size of a 
	      int na;	       Column size of a 
	      int rstart,rend; Row range of submatrix b: rstart:rend 
	      int cstart,cend; Column range of submatrix b: cstart:cend 

    Author: Clay Thompson, MATLAB

   History: 15-Mar-95 lrf conversion to dv routine 
            (stripped out all the MATLAB stuff)
*/

double *dvsubmat(double *b, double *a, int ma, int na, int rstart, int rend, int cstart, int cend)
{
  double *p,*q;			/* Pointers to elements in a and b */
  int	 i,j;			/* Loop counters */
  int	 mb,nb,step;

  /* Size of result array */
  mb = rend - rstart + 1;
  nb = cend - cstart + 1;

  /* Increments and pointers */

  step = ma - mb;
  q    = b;
  p    = a + rstart + cstart*ma;

  /* Copy elements from subsection of a to b */

  for (j=0;j<nb;++j) {
    for (i=0;i<mb;++i) {
      *(q++) = *(p++);
    }
    p += step;
  }
  return(b);
}
/**<dvxcor1.c>*******************************************************

     Title: Double Vector CROSEE CORrelation in 1-dimension

   Purpose: Cross-correlates two vectors x and y and returns
            vector z, either the same size as x, or full size.

      Call: double *dvxcor1(z, x, y, nx, ny, shape)
              double *z	         Result matrix (nx)-by-(1) 
	      double *x;	 Larger vector 
	      double *y;	 Smaller vector 
	      int    nx;	 Length of x
	      int    ny;	 Length of y
	      char   *shape      Output vector size
	                         shape is either 'same' or 'full'

   Returns pointer to z;

    Author: Lawrence Frank

   History: 04-Jan-96

*/
double *dvxcor1(double *z, double *x, double *y, int nx, int ny, char *shape)
{
  dvcor2(z,x,y,1,nx,ny,1,shape);
  return(z);
}

/**<dvcor2.c>*******************************************************

     Title: Double Vector Correlation in 2-dimensions

   Purpose: Correlates two matrices x and y and returns
            vector z, the same size as x.

      Call: double *dvcor2(c, a, b, ma, na, mb, nb, shape)
              double  *c         Result matrix (ma+mb-1)-by-(na+nb-1) 
	      double  *a;	 Larger matrix 
	      double  *b;	 Smaller matrix 
	      int     ma;	 Row size of a 
	      int     na;	 Column size of a 
	      int     mb;	 Row size of b 
	      int     nb;	 Column size of b 
	      char *shape        Output matrix size
	                         shape is either 'same' or 'full'

   Returns pointer to c;

    Author: Lawrence R. Frank 
            (based on Clay Thompson MATLAB routine)

   History: 04-Jan-96 lrf conversion of dvconv2.c

*/
double *dvcor2(double *c, double *a, double *b, int ma, int na, int mb, int nb, char *shape)
{
  void dvcor2b(double *c, double *a, double *b, int ma, int na, int mb, int nb);
  double *x;
  int mc     = ma;
  int nc     = na;
  int rstart = mb/2;
  int rend   = mb/2 + mc - 1;
  int cstart = nb/2;
  int cend   = nb/2 + nc - 1;
  int nx     = (ma+mb-1)*(na+nb-1);

  /* Allocate space for full correlation */
  x = (double *) calloc(nx,sizeof(double));

  /* Perform full correlation */
  dvcor2b(x,a,b,ma,na,mb,nb);

  /* Return desired shape */
  if(strcmp(shape,"full")==0) {
    dvmove(c,x,nx);
  } else if(strcmp(shape,"same")==0) {
    dvsubmat(c, x, ma, na, rstart, rend, cstart, cend);
  } else {
    perror("Illegal shape value!");
    exit(0);
  }
  free(x);
  return(c);
}
/*
     double *c;	         Result matrix (ma+mb-1)-by-(na+nb-1) 
     double *a;	         Larger matrix 
     double *b;	         Smaller matrix 
     int    ma;		 Row size of a 
     int    na;		 Column size of a 
     int    mb;		 Row size of b 
     int    nb;		 Column size of b 
*/
void dvcor2b(double *c, double *a, double *b, int ma, int na, int mb, int nb)
{
  register double *p,*q; /* Pointer to elements in 'a' and 'c' matrices */
  register double w;	 /* Weight (element of 'b' matrix) */
  register int k,l,i,j;
  int          mc,nc;
  double       *r;	 /* Pointer to elements in 'b' matrix */
	
  mc = ma+mb-1;
  nc = na+nb-1;
	
  /* Perform correlation */

  r = b+nb*mb-1;	
  for (j=0; j<nb; ++j) {		/* For each non-zero element in b */
    for (i=0; i<mb; ++i) {
      w = *(r--);			/* Get weight from b matrix */
      if (w != 0.0) {
	p = c + i + j*mc;	        /* Start at first column of a in c. */
	for (l=0, q=a; l<na; l++) {	/* For each column of a ... */
	  for (k=0; k<ma; k++) {	
	    *(p++) += *(q++) * w;      /* multiply by weight and add. */
	  }
	  p += mb - 1;	/* Jump to next column position of a in c */
	}
      } /* end if */
    } 
  }
}

/**<vcor.c>****************************************************

     Title: Float Vector CORrelation

   Purpose: Correlates two float vectors quickly

      Call: float vcor( float *x, float *y, int n )

     Input: x,y = input float vectors
	    n = length(x)

    Output: z = x[0]*cnjg(y[0]) +...+ x[N-1]*cnjg(y[N-1])

     Notes: Internal calculation is float precision,
            and the sum is returned in float precision.

   History: Musicus, circa 1989
*/
float vcor(float *x, float *y, int n)
{
  float sum=0.;

  if(n>0) {

      do {
	sum += *x++ * *y ;
	y++ ;
      } while(--n > 0) ;

  }

  return(sum);
}

/**<dvcor.c>****************************************************

     Title: Double Vector CORrelation

   Purpose: Correlates two double vectors quickly

      Call: double dvcor( double *x, double *y, int n )

     Input: x,y = input double vectors
	    n = length(x)

    Output: z = x[0]*cnjg(y[0]) +...+ x[N-1]*cnjg(y[N-1])

     Notes: Internal calculation is double precision,
            and the sum is returned in double precision.

   History: Musicus, circa 1989
*/
double dvcor(double *x, double *y, int n)
{
  double sum=0.;

  if(n>0) {

      do {
	sum += *x++ * *y ;
	y++ ;
      } while(--n > 0) ;

  }

  return(sum);
}

/**<zvcor.c>****************************************************

     Title: Complex Vector CORrelation

   Purpose: Correlates two complex vectors

      Call: complex zvcor( complex *x, complex *y, int n )

     Input: x,y = input complex vectors
	    n = length(x)

    Output: z = x[0]*cnjg(y[0]) +...+ x[N-1]*cnjg(y[N-1])

     Notes: Handles case where x=y specially, taking half the time

   History: Musicus, circa 1989
*/
complex zvcor( complex *x, complex *y, int n )
{
	float *px = (float*)x ;
	float *py = (float*)y ;

	double sumr = 0.;
	double sumi = 0.;
	complex ztemp ;

	if(n>0)
	{
		if(px==py)
		{
			do {
				sumr += *px++ * *py++ ;
				sumr += *px++ * *py++ ;
			} while(--n > 0) ;
		}
		else
		{
			do {
				sumr += *px * *py ;
				sumi -= *px++ * *(py+1) ;
				sumi += *px * *py++ ;
				sumr += *px++ * *py++ ;
			} while(--n > 0) ;
		}
	}

	ztemp.r = sumr ;
	ztemp.i = sumi ;
	return(ztemp) ;
}

/**<dzcor.c>****************************************************

     Title: Double Complex CORrelation

   Purpose: Correlates two complex numbers

      Call: dcomplex dzcor( dcomplex a, dcomplex b )

     Input: a,b = dcomplex numbers

    Output: z = x*cnjg(y) 

   History: lrf 8-Oct-07
*/
dcomplex dzcor( dcomplex a, dcomplex b )
{

  return dzmul(a,dzcnjg(b)) ;

}
/**<dzvcor.c>****************************************************

     Title: Double Complex Vector CORrelation

   Purpose: Correlates two complex vectors

      Call: dcomplex dzvcor( dcomplex *x, dcomplex *y, int n )

     Input: x,y = input dcomplex vectors
	    n = length(x)

    Output: z = x[0]*cnjg(y[0]) +...+ x[N-1]*cnjg(y[N-1])

     Notes: Handles case where x=y specially, taking half the time

   History: Musicus, circa 1989
*/
dcomplex dzvcor( dcomplex *x, dcomplex *y, int n )
{
	double *px = (double*)x ;
	double *py = (double*)y ;

	double sumr = 0.;
	double sumi = 0.;
	dcomplex dztemp ;

	if(n>0)
	{
		if(px==py)
		{
			do {
				sumr += *px++ * *py++ ;
				sumr += *px++ * *py++ ;
			} while(--n > 0) ;
		}
		else
		{
			do {
				sumr += *px * *py ;
				sumi -= *px++ * *(py+1) ;
				sumi += *px * *py++ ;
				sumr += *px++ * *py++ ;
			} while(--n > 0) ;
		}
	}

	dztemp.dr = sumr ;
	dztemp.di = sumi ;
	return(dztemp) ;
}

/**<zvconv.c>****************************************************

     Title: Complex Vector CONVolution

   Purpose: Convolves two complex vectors

      Call: complex zvconv( complex *x, complex *y, int n )

     Input: x,y = input complex vectors
	    n = length(x)

    Output: zvconv(x,y+k,N) = x[0]*cnjg(y[k]) +...+ x[N-1]*cnjg(y[k-N+1])

   History: Musicus, circa 1989
*/
complex zvconv( complex *x, complex *y, int n )
{
	float *px = (float*)x ;
	float *py = (float*)y ;

	double sumr = 0.;
	double sumi = 0.;
	complex ztemp ;

	if(n>0)
	{
		py+=2 ;		/* position just beyond array */
		do {
			sumi -= *px * *--py ;
			sumr += *px++ * *--py ;
			sumi += *px * *py ;
			sumr += *px++ * *(py+1) ;
		} while(--n > 0) ;
	}

	ztemp.r = sumr ;
	ztemp.i = sumi ;
	return(ztemp) ;
}

/**<dzvconv.c>****************************************************

     Title: Double Complex Vector CONVolution

   Purpose: Convolves two dcomplex vectors

      Call: dcomplex dzvconv( dcomplex *x, dcomplex *y, int n )

     Input: x,y = input dcomplex vectors
	    n = length(x)

    Output: dzvconv(x,y+k,N) = x[0]*cnjg(y[k]) +...+ x[N-1]*cnjg(y[k-N+1])

   History: Musicus, circa 1989
*/
dcomplex dzvconv( dcomplex *x, dcomplex *y, int n )
{
	double *px = (double*)x ;
	double *py = (double*)y ;

	double sumr = 0.;
	double sumi = 0.;
	dcomplex dztemp ;

	if(n>0)
	{
		py+=2 ;		/* position just beyond array */
		do {
			sumi -= *px * *--py ;
			sumr += *px++ * *--py ;
			sumi += *px * *py ;
			sumr += *px++ * *(py+1) ;
		} while(--n > 0) ;
	}

	dztemp.dr = sumr ;
	dztemp.di = sumi ;
	return(dztemp) ;
}

/**<zvrev.c>****************************************************

     Title: Complex Vector REVersal

   Purpose: Reverses vector: z[i]=x[n-1-i] for i=0,...,n-1

      Call: complex zvrev( complex *x, complex *y, int n )

     Input: x,y = input complex vectors
	    n = length(x)

    Output: Returns pointer to z[]

     Notes: Works correctly if x[] and z[] vectors do not overlap at all, 
            or if they overlap completely

   History: Musicus, circa 1989
*/
complex *zvrev( complex *z, complex *x, int N )
{
	complex *zold = z;
	complex *ze = z+N ;
	complex *xe = x+N ;
	complex temp1,temp2;

	if(N>0)
	{
		N = (N+1)>>1 ;
		do {
			temp1 = *x++ ;
			temp2 = *--xe ;
			*z++  = temp2 ;
			*--ze = temp1 ;
		} while( --N > 0) ;
	}
	return(zold) ;
}

/* Recursive split technique
   Split first and last halves into reals and imaginaries
   Then merge the two lists */

static void split(complex *z, int N)
{
	int Ndn,Nup ;
	float *upptr,*dnptr ;
	float temp ;

	if(N==2)		/* Handle specially */
	{
		temp = z[0].i ;
		z[0].i = z[1].r ;
		z[1].r = temp ;
	}
	else if(N>2)
	{
		Ndn = N/2 ;
		Nup = (N+1)/2 ;
		split(z,Nup) ;		/* split first half */
		split(z+Nup,Ndn) ;	/* split second half */
		upptr = (float*)(z+Nup) ;
		dnptr = ((float*)z)+Nup ;
		if(N&1)
		{		/* Handle odd lengths */
			temp = *--upptr ;	/* save last imaginary
						of first half */
			do	/* Exchange 1st half imag, 2nd half reals */
			{
				*upptr++ = *dnptr ;
				*dnptr++ = *upptr ;
			} while(--Ndn > 0) ;
			*upptr = temp ;		/* Last imag of 1st half */
		}
		else
		{		/* Handle even lengths */
			do	/* Exchange 1st half imag, 2nd half reals */
			{
				temp = *dnptr ;
				*dnptr++ = *upptr ;
				*upptr++ = temp ;
			} while(--Nup > 0) ;
		}
	}
	return ;
}

/* This routine fills the complex array from real and imaginary parts
    using a recursive conquer and divide algorithm.
  Handles N=1,2 specially.  Anything larger, it splits into two halves,
    roughly sorts out the two halves, then calls itself recursively on
    each half.
*/
static void unsplit(float *z , int N)
{
	float *px ;
	float *py ;
	float tmp ;
	int N2dn,N2up ;

	if(N<=1) return ;

	if(N&1)	/* odd */
	{
		N2dn = N>>1 ;
		N2up = N2dn+1 ;
		px = z+N ;	/* just above reals */
		py = px+N2up ;	/* halfway through imags */
		tmp = *(py-1) ;	/* save middle imag */
		do {
			*--py = *--px ;	/* move real up */
			*px = *(py-1) ;	/* move imag down */
		} while(--N2dn > 0) ;
		*--py = tmp ;	/* move last imag */
	}
	else	/* even */
	{
		N2dn = N>>1 ;
		px = z+N2dn ;
		py = z+N ;
		do {		/* exchange half the reals and imags */
			tmp = *px ;
			*px++ = *py ;
			*py++ = tmp ;
		} while(--N2dn > 0) ;
	}

	/* recursively split bottom and top halves */
	N2dn = N>>1 ;
	N2up = (N+1)>>1 ;
	if(N2up>1)
		unsplit(z,N2up) ;
	if(N2dn>1)
		unsplit(z+(N2up<<1),N2dn) ;

	return ;
}
/*      float *dvtov(z,x,n)	Converts double vector x[] into float z[]
 *      float *z ;		  z[i] = x[i] for i=0,...,n-1 if z<x
 *	double *x ;			      for i=n-1,...,0 if z>=x
 *	int n ;			Returns pointer to z[]
 *
 * Works correctly if z[] and x[] overlap provided that z[] is not wholly
 *   contained in the interior of x[] - i.e. z<=x or &x[N]<=&z[N]
*/
float *dvtov(float *z,double *x, int N)
{
	float *zold = z;

	if(N>0)
	{
		if( z <= (float*)x )
		{
			do {
				*z++ = *x++ ;
			} while( --N > 0) ;
		}
		else
		{
			z += N ;
			x += N ;
			do {
				*--z = *--x ;
			} while( --N > 0) ;
		}
	}
	return(zold) ;
}

/*      dcomplex *zvtodzv(z,x,n) Converts complex vector x[] into dcomplex z[]
 *      dcomplex *z ;		  z[i] = x[i] for i=0,...,n-1 if z<x
 *	complex *x ;			      for i=n-1,...,0 if z>=x
 *	int n ;			 Returns pointer to z[]
 *
 * Works correctly if z[] and x[] overlap provided that x[] is not wholly
 *   contained in the interior of z[] - i.e. x<=z or &z[N]<=&x[N]
*/
dcomplex *zvtodzv(dcomplex *z,complex *x,int N)
{
	dcomplex *zold = z;

	if(N>0)
	{
		if( (complex*)z < x )
		{
			do {
				z->dr = x->r ;
				z->di = x->i ;
				z++ ; x++ ;
			} while( --N > 0) ;
		}
		else
		{
			z += N ;
			x += N ;
			do {
				z-- ; x-- ;
				z->di = x->i ;
				z->dr = x->r ;
			} while( --N > 0) ;
		}
	}
	return(zold) ;
}

/**<zvtori.c>****************************************************

     Title: Complex Vector TO Real and Imaginary

   Purpose: Splits complex vector into real and imaginary parts:
            for(i=0 ; i<N ; i++) { real[i]=z[i].r ; imag[i]=z[i].i ; }

      Call: complex *zvtori(float *real, float *imag, complex *z, int N)

    Output: Returns pointer to z[] array

     Notes: If real=NULL or imag=NULL, then zvtori discards 
            the real or imag part
	    If both real=NULL and imag=NULL then zvtori 
	    splits the complex vector in place using a 
	    conquer and divide algorithm (time NlogN instead of N)

   History: Musicus, circa 1989
*/
complex *zvtori(float *real, float *imag, complex *z, int N)
{
	float *zz = (float*)z ;

	if(N>0)
	{
		if(real)
		{
			if(imag)
			{
				do
				{
					*real++ = *zz++ ;
					*imag++ = *zz++ ;
				} while(--N > 0) ;
			}
			else
			{
				do
				{
					*real++ = *zz++ ;
					zz++ ;
				} while(--N > 0) ;
			}
		}
		else
		{
			if(imag)
			{
				do
				{
					zz++ ;
					*imag++ = *zz++ ;
				} while(--N > 0) ;
			}
			else
			{
				split(z,N) ;
			}
		}
	}
	return(z) ;
}

/* Recursive split technique
   Split first and last halves into reals and imaginaries
   Then merge the two lists
*/

static void dsplit(dcomplex *z, int N)
{
	int Ndn,Nup ;
	double *upptr,*dnptr ;
	double temp ;

	if(N==2)		/* Handle specially */
	{
		temp = z[0].di ;
		z[0].di = z[1].dr ;
		z[1].dr = temp ;
	}
	else if(N>2)
	{
		Ndn = N/2 ;
		Nup = (N+1)/2 ;
		dsplit(z,Nup) ;		/* split first half */
		dsplit(z+Nup,Ndn) ;	/* split second half */
		upptr = (double*)(z+Nup) ;
		dnptr = ((double*)z)+Nup ;
		if(N&1)
		{		/* Handle odd lengths */
			temp = *--upptr ;	/* save last imaginary
						of first half */
			do	/* Exchange 1st half imag, 2nd half reals */
			{
				*upptr++ = *dnptr ;
				*dnptr++ = *upptr ;
			} while(--Ndn > 0) ;
			*upptr = temp ;		/* Last imag of 1st half */
		}
		else
		{		/* Handle even lengths */
			do	/* Exchange 1st half imag, 2nd half reals */
			{
				temp = *dnptr ;
				*dnptr++ = *upptr ;
				*upptr++ = temp ;
			} while(--Nup > 0) ;
		}
	}
	return ;
}
/* This routine fills the dcomplex array from real and imaginary parts
    using a recursive conquer and divide algorithm.
  Handles N=1,2 specially.  Anything larger, it splits into two halves,
    roughly sorts out the two halves, then calls itself recursively on
    each half.
*/
static void dunsplit(double *z, int N)
{
	double *px ;
	double *py ;
	double tmp ;
	int N2dn,N2up ;

	if(N<=1) return ;

	if(N&1)	/* odd */
	{
		N2dn = N>>1 ;
		N2up = N2dn+1 ;
		px = z+N ;	/* just above reals */
		py = px+N2up ;	/* halfway through imags */
		tmp = *(py-1) ;	/* save middle imag */
		do {
			*--py = *--px ;	/* move real up */
			*px = *(py-1) ;	/* move imag down */
		} while(--N2dn > 0) ;
		*--py = tmp ;	/* move last imag */
	}
	else	/* even */
	{
		N2dn = N>>1 ;
		px = z+N2dn ;
		py = z+N ;
		do {		/* exchange half the reals and imags */
			tmp = *px ;
			*px++ = *py ;
			*py++ = tmp ;
		} while(--N2dn > 0) ;
	}

/* recursively split bottom and top halves */
	N2dn = N>>1 ;
	N2up = (N+1)>>1 ;
	if(N2up>1)
		dunsplit(z,N2up) ;
	if(N2dn>1)
		dunsplit(z+(N2up<<1),N2dn) ;

	return ;
}

/**<dzvtori.c>****************************************************

     Title: Double Complex Vector TO Real and Imaginary

   Purpose: Splits dcomplex vector into real and imaginary parts:
            for(i=0 ; i<N ; i++) { real[i]=z[i].r ; imag[i]=z[i].i ; }

      Call: dcomplex *dzvtori(double *real, double *imag, dcomplex *z, int N)

    Output: Returns pointer to z[] array

     Notes: If real=NULL or imag=NULL, then dzvtori discards 
            the real or imag part
	    If both real=NULL and imag=NULL then dzvtori 
	    splits the dcomplex vector in place using a 
	    conquer and divide algorithm (time NlogN instead of N)

   History: Musicus, circa 1989
*/
dcomplex *dzvtori(double *real, double *imag, dcomplex *z, int N)
{
	double *zz = (double*)z ;

	if(N>0)
	{
		if(real)
		{
			if(imag)
			{
				do
				{
					*real++ = *zz++ ;
					*imag++ = *zz++ ;
				} while(--N > 0) ;
			}
			else
			{
				do
				{
					*real++ = *zz++ ;
					zz++ ;
				} while(--N > 0) ;
			}
		}
		else
		{
			if(imag)
			{
				do
				{
					zz++ ;
					*imag++ = *zz++ ;
				} while(--N > 0) ;
			}
			else
			{
				dsplit(z,N) ;
			}
		}
	}
	return(z) ;
}

/**<zvmul.c>****************************************************

     Title: Complex Vector MULtiplication

   Purpose: Computes z[i]=x[i]*y[i]  for i=0,...,N-1

      Call: complex zvmul( complex *z, complex *x, complex *y, int n )

    Output: Returns pointer to z[]

   History: Musicus, circa 1989
            02_05_04 lrf - make in-place
*/
complex *zvmul( complex *z, complex *x, complex *y, int n )
{
	float *pz = (float*)z ;
	complex ztmp;

	if(n>0)
	{
		do {
		        ztmp.r = x->r*y->r - x->i*y->i ;
			ztmp.i = x->r*y->i + x->i*y->r ;
			*pz++  = ztmp.r ;
			*pz++  = ztmp.i ;
			x++ ; y++ ;
		} while( --n > 0) ;
	}
	return(z) ;
}

 /**<dvmul.c>****************************************************

      Title: Double Vector MULtiplication

    Purpose: Computes z[i]=x[i]*y[i]  for i=0,...,N-1

       Call: double zvmul( double *z, double *x, double *y, int n )

     Output: Returns pointer to z[]

    History: Musicus, circa 1989
	     02_05_04 lrf - make in-place
 */
double *dvmul(double *z,double *x,double *y,int n)
{
	double *zold = z ;

	if(n>0)
	{
		do {
			*z++ = *x++ * *y ;
			y++ ;
		} while( --n > 0) ;
	}
	return(zold) ;
}
/**<dzvmul.c>****************************************************

     Title: Doubble Complex Vector MULtiplication

   Purpose: Computes z[i]=x[i]*y[i]  for i=0,...,N-1

      Call: dcomplex dzvmul( dcomplex *z, dcomplex *x, dcomplex *y, int n )

    Output: Returns pointer to z[]

   History: Musicus, circa 1989
            02_05_04 lrf - make in-place
*/
dcomplex *dzvmul( dcomplex *z, dcomplex *x, dcomplex *y, int n )
{
	double *pz = (double*)z ;

	if(n>0)
	{
		do {
			*pz++ = x->dr*y->dr - x->di*y->di ;
			*pz++ = x->dr*y->di + x->di*y->dr ;
			x++ ; y++ ;
		} while( --n > 0) ;
	}
	return(z) ;
}

/**<zvmulc.c>****************************************************

     Title: Complex Vector Conjugate MULtiplication

   Purpose: Computes z[i]=x[i]*cnjg(y[i])  for i=0,...,N-1

      Call: complex zvmulc( complex *z, complex *x, complex *y, int n )

    Output: Returns pointer to z[]

   History: Musicus, circa 1989
            02_05_04 lrf - make in-place
*/
complex *zvmulc( complex *z, complex *x, complex *y, int n )
{
	float   *pz = (float*)z ;
	complex ztmp;

	if(n>0)
	{
		do {
			ztmp.r = x->r*y->r + x->i*y->i ;
			ztmp.i = x->i*y->r - x->r*y->i ;
			*pz++  = ztmp.r ;
			*pz++  = ztmp.i ;
			x++ ; y++ ;
		} while( --n > 0) ;
	}
	return(z) ;
}

/**<dzvmulc.c>****************************************************

     Title: Complex Vector Conjugate MULtiplication

   Purpose: Computes z[i]=x[i]*cnjg(y[i])  for i=0,...,N-1

      Call: complex dzvmulc( complex *z, complex *x, complex *y, int n )

    Output: Returns pointer to z[]

   History: Musicus, circa 1989
            02_05_04 lrf - make in-place
*/
dcomplex *dzvmulc( dcomplex *z, dcomplex *x, dcomplex *y, int n )
{
	float   *pz = (float*)z ;
	dcomplex ztmp;

	if(n>0)
	{
		do {
			ztmp.dr = x->dr*y->dr + x->di*y->di ;
			ztmp.di = x->di*y->dr - x->dr*y->di ;
			*pz++  = ztmp.dr ;
			*pz++  = ztmp.di ;
			x++ ; y++ ;
		} while( --n > 0) ;
	}
	return(z) ;
}

/**<vmnmx>************************************************************

    Title: Float Vector MiN/MaX

  Purpose: Find the minimum and maximum values,
           and their respective indices, of a float vector.
 
     Call: float vmnmx(float *x, int n, int *nmin, int *nmax)

    Input: x = float vector to examine
           n = length(x)

   Output:  *nmin = index of mininum value in x[]
            *nmax = index of maxinum value in x[]

	    Also value of function is the maximum of the
	    absolute values of the vector points.
	    If nmin=0 or nmax=0, then does not return these pointer values
	    If n<=0, returns immediately with value 0 and *nmin = *nmax = 0.

  History: Musicus, circa 1980
*/
float vmnmx(float *x, int n, int *nmin, int *nmax)
{
#define MABS(x)		(((x)<0) ? -(x) : (x))
#define MMAX(x,y)	(((x)<(y)) ? (y) : (x))

	float *xold = x ;
	float xmin,xmax ;
	int nnmin,nnmax ;
	int count ;

	if(n<=0)
	{
		if(nmin) *nmin = 0 ;
		if(nmax) *nmax = 0 ;
		return(0.) ;
	}
	nnmin = nnmax = 0 ;
	xmin = xmax = *x++ ;
	n-- ;
	count = 1 ;
	if(n>0)
	{
		do {
			if(*x < xmin) {
				nnmin = count++ ;
				xmin = *x++ ;
			} else if(*x > xmax) {
				nnmax = count++ ;
				xmax = *x++ ;
			} else {
				count++ ;
				x++ ;
			}
		} while( --n > 0 ) ;
	}
	if(nmin) *nmin = nnmin ;
	if(nmax) *nmax = nnmax ;
	xmax = MABS(xmax) ;
	xmin = MABS(xmin) ;
	return(MMAX(xmax,xmin)) ;
}

/**<vran>**************************************************************

    Title: Float Vector RANge

  Purpose: Find the minimum and maximum values,
           and their respective indices, of a float vector.
 
    Input: a = float vector
           n = length(x)

   Output:  amin (float) = min(a)
           lamin   (int) = location of amin in a[] 
	    amax (float) = max(a)
           lamax   (int) = location of amax in a[] 

    Notes: call must use addresses!
           vran(a,&amin,&lamin,&amax,&lamax,n);
*/
void vran(float *a, float *amin, int *lamin, float *amax, int *lamax, int n)
{
  vmnmx(a,n,lamin,lamax);

  *amin = *(a + *lamin);
  *amax = *(a + *lamax); 
}

/**<print_vran.c>****************************************************

     Title: PRINT Float Vector RANge

   Purpose: Report float vector min,max and their indeces

      Call: void print_vran(float *x, int nsamp, char *text)

     Input: x = float vector
	    n = length(x)
	    text = informational text

    Output: none

    History: lrf 11/2005
*/
void print_vran(float *x, int nsamp, char *text)
{
  float xmin,xmax;
  int lxmin,lxmax;

  vran(x,&xmin,&lxmin,&xmax,&lxmax,nsamp);

  fprintf(stderr,"%s: min = %f at index %i\n",text,xmin,lxmin);
  fprintf(stderr,"%s: max = %f at index %i\n",text,xmax,lxmax);
}

/**<ivmnmx>************************************************************

    Title: Int Vector MiN/MaX

  Purpose: Find the minimum and maximum values,
           and their respective indices, of a int vector.
 
     Call: int ivmnmx(int *x, int n, int *nmin, int *nmax)

    Input: x = int vector to examine
           n = length(x)

   Output:  *nmin = index of mininum value in x[]
            *nmax = index of maxinum value in x[]

	    Also value of function is the maximum of the
	    absolute values of the vector points.
	    If nmin=0 or nmax=0, then does not return these pointer values
	    If n<=0, returns immediately with value 0 and *nmin = *nmax = 0.

  History: lrf, Sept 11, 2006 (from vmnmx)
*/
int ivmnmx(int *x, int n, int *nmin, int *nmax)
{
#define MABS(x)		(((x)<0) ? -(x) : (x))
#define MMAX(x,y)	(((x)<(y)) ? (y) : (x))

	int *xold = x ;
	int xmin,xmax ;
	int nnmin,nnmax ;
	int count ;

	if(n<=0)
	{
		if(nmin) *nmin = 0 ;
		if(nmax) *nmax = 0 ;
		return(0) ;
	}
	nnmin = nnmax = 0 ;
	xmin = xmax = *x++ ;
	n-- ;
	count = 1 ;
	if(n>0)
	{
		do {
			if(*x < xmin) {
				nnmin = count++ ;
				xmin = *x++ ;
			} else if(*x > xmax) {
				nnmax = count++ ;
				xmax = *x++ ;
			} else {
				count++ ;
				x++ ;
			}
		} while( --n > 0 ) ;
	}
	if(nmin) *nmin = nnmin ;
	if(nmax) *nmax = nnmax ;
	xmax = MABS(xmax) ;
	xmin = MABS(xmin) ;
	return(MMAX(xmax,xmin)) ;
}

/**<ivran>**************************************************************

    Title: Int Vector RANge

  Purpose: Find the minimum and maximum values,
           and their respective indices, of a int vector.
 
    Input: a = int vector
           n = length(x)

   Output:  amin (int) = min(a)
           lamin   (int) = location of amin in a[] 
	    amax (int) = max(a)
           lamax   (int) = location of amax in a[] 

    Notes: call must use addresses!
           ivran(a,&amin,&lamin,&amax,&lamax,n);
*/
void ivran(int *a, int *amin, int *lamin, int *amax, int *lamax, int n)
{
  ivmnmx(a,n,lamin,lamax);

  *amin = *(a + *lamin);
  *amax = *(a + *lamax); 
}

/**<print_ivran.c>****************************************************

     Title: PRINT Int Vector RANge

   Purpose: Report int vector min,max and their indeces

      Call: void print_ivran(int *x, int nsamp, char *text)

     Input: x = int vector
	    n = length(x)
	    text = informational text

    Output: none

    History: lrf 11/2005
*/
void print_ivran(int *x, int nsamp, char *text)
{
  int xmin,xmax;
  int lxmin,lxmax;

  ivran(x,&xmin,&lxmin,&xmax,&lxmax,nsamp);

  fprintf(stderr,"%s: min = %i at index %i\n",text,xmin,lxmin);
  fprintf(stderr,"%s: max = %i at index %i\n",text,xmax,lxmax);
}

/**<dvmnmx>************************************************************

    Title: Double Vector MiN/MaX

  Purpose: Find the minimum and maximum values,
           and their respective indices, of a double vector.
 
     Call: double dvmnmx(double *x, int n, int *nmin, int *nmax)

    Input: x = double vector to examine
           n = length(x)

   Output:  *nmin = index of mininum value in x[]
            *nmax = index of maxinum value in x[]

	    Also value of function is the maximum of the
	    absolute values of the vector points.
	    If nmin=0 or nmax=0, then does not return these pointer values
	    If n<=0, returns immediately with value 0 and *nmin = *nmax = 0.

  History: Musicus, circa 1980
*/
double dvmnmx(double *x, int n, int *nmin, int *nmax)
{
#define MABS(x)		(((x)<0) ? -(x) : (x))
#define MMAX(x,y)	(((x)<(y)) ? (y) : (x))

	double *xold = x ;
	double xmin,xmax ;
	int nnmin,nnmax ;
	int count ;

	if(n<=0)
	{
		if(nmin) *nmin = 0 ;
		if(nmax) *nmax = 0 ;
		return(0.) ;
	}
	nnmin = nnmax = 0 ;
	xmin = xmax = *x++ ;
	n-- ;
	count = 1 ;
	if(n>0)
	{
		do {
			if(*x < xmin) {
				nnmin = count++ ;
				xmin = *x++ ;
			} else if(*x > xmax) {
				nnmax = count++ ;
				xmax = *x++ ;
			} else {
				count++ ;
				x++ ;
			}
		} while( --n > 0 ) ;
	}
	if(nmin) *nmin = nnmin ;
	if(nmax) *nmax = nnmax ;
	xmax = MABS(xmax) ;
	xmin = MABS(xmin) ;
	return(MMAX(xmax,xmin)) ;
}

/**<dvran>**************************************************************

    Title: Double Vector RANge

  Purpose: Find the minimum and maximum values,
           and their respective indices, of a double vector.
 
    Input: a = double vector
           n = length(x)

   Output:  amin (double) = min(a)
           lamin   (int) = location of amin in a[] 
	    amax (double) = max(a)
           lamax   (int) = location of amax in a[] 

    Notes: call must use addresses!
           dvran(a,&amin,&lamin,&amax,&lamax,n);
*/
void dvran(double *a, double *amin, int *lamin, double *amax, int *lamax, int n)
{
  dvmnmx(a,n,lamin,lamax);

  *amin = *(a + *lamin);
  *amax = *(a + *lamax); 
}


/**<print_dvran.c>****************************************************

     Title: PRINT Double Vector RANge

   Purpose: Report double vector min,max and their indeces

      Call: void print_dvran(double *x,int n,char *text)

     Input: x = double vector
	    n = length(x)
	    text = informational text

    Output: none

    History: lrf 11/2005
*/
void print_dvran(double *x, int n, char *text)
{
  double xmin,xmax;
  int lxmin,lxmax;

  dvran(x,&xmin,&lxmin,&xmax,&lxmax,n);

  fprintf(stderr,"%s: min = %f at index %i\n",text,xmin,lxmin);
  fprintf(stderr,"%s: max = %f at index %i\n",text,xmax,lxmax);
}

/**<zvmnmx>************************************************************

    Title: Complex Vector MiN/MaX

  Purpose:  Returns the indices (*nmin and *nmax) of the 
            minimum and maximum absolute value in the complex vector x[].

     Call: double vmnmx(complex *x, int n, int *nmin, int *nmax)

    Input: x = complex vector to examine
           n = length(x)

   Output:  *nmin = index of mininum value in x[]
            *nmax = index of maxinum value in x[]

	    Also value of function is the maximum of the
	    absolute values of the vector points.
	    If nmin=0 or nmax=0, then does not return these pointer values
	    If n<=0, returns immediately with value 0 and *nmin = *nmax = 0.

  History: Musicus, circa 1980
*/
double zvmnmx(complex *z, int N, int *nmin, int *nmax)
{
	float *zptr = (float*) z ;
	double absr, absi, ratio, biggie, smallie;
	int NN ;

	if(nmin) *nmin = 0 ;
	if(nmax) *nmax = 0 ;
	if(N<=0)
	{
		return(0.) ;
	}

	biggie = -1 ; /* flag first entry */
	for(NN=0 ; NN<N ; NN++)
	{
		absr = fabs(*zptr++) ;
		absi = fabs(*zptr++) ;
		if( absr <= absi && absi>0.0 ) {
			ratio = absr/absi;
			absr = absi * sqrt(1 + ratio*ratio); 
		}
		else if (absr>0.0) {
			ratio = absi/absr;
			absr = absr * sqrt(1 + ratio*ratio);
		}
		else {
			absr = 0.0;
		}
		if(biggie<0) {
			biggie = smallie = absr ;
		}
		else if(absr>biggie) {
			biggie = absr ;
			if(nmax) *nmax = NN ;
		}
		else if(absr<smallie) {
			smallie = absr ;
			if(nmin) *nmin = NN ;
		}
	}
	return(biggie) ;
}
/**<dzvmnmx>************************************************************

    Title: Double Complex Vector MiN/MaX

  Purpose:  Returns the indices (*nmin and *nmax) of the 
            minimum and maximum absolute value in the complex vector x[].

     Call: double vmnmx(dcomplex *x, int n, int *nmin, int *nmax)

    Input: x = dcomplex vector to examine
           n = length(x)

   Output:  *nmin = index of mininum value in x[]
            *nmax = index of maxinum value in x[]

	    Also value of function is the maximum of the
	    absolute values of the vector points.
	    If nmin=0 or nmax=0, then does not return these pointer values
	    If n<=0, returns immediately with value 0 and *nmin = *nmax = 0.

  History: Musicus, circa 1980
*/
double dzvmnmx(dcomplex *z,int N,int *nmin,int *nmax)
{
	double *zptr = (double*) z ;
	double absr, absi, ratio, biggie, smallie;
	int NN ;

	if(nmin) *nmin = 0 ;
	if(nmax) *nmax = 0 ;
	if(N<=0)
	{
		return(0.) ;
	}

	biggie = -1 ; /* flag first entry */
	for(NN=0 ; NN<N ; NN++)
	{
		absr = fabs(*zptr++) ;
		absi = fabs(*zptr++) ;
		if( absr <= absi && absi>0.0 ) {
			ratio = absr/absi;
			absr = absi * sqrt(1 + ratio*ratio); 
		}
		else if (absr>0.0) {
			ratio = absi/absr;
			absr = absr * sqrt(1 + ratio*ratio);
		}
		else {
			absr = 0.0;
		}
		if(biggie<0) {
			biggie = smallie = absr ;
		}
		else if(absr>biggie) {
			biggie = absr ;
			if(nmax) *nmax = NN ;
		}
		else if(absr<smallie) {
			smallie = absr ;
			if(nmin) *nmin = NN ;
		}
	}
	return(biggie) ;
}

/**<zvran>**************************************************************

    Title: Complex Vector RANge

  Purpose: Find the minimum and maximum values,
           and their respective indices, of a complex vector.
 
    Input: a = complex vector
           n = length(x)

   Output:  amin (float) = min(a)
           lamin   (int) = location of amin in a[] 
	    amax (float) = max(a)
           lamax   (int) = location of amax in a[] 

    Notes: call must use addresses!
           zvran(a,&amin,&lamin,&amax,&lamax,n);
*/
void zvran(complex *a, float *amin, int *lamin, float *amax, int *lamax, int n)
{
  complex zmin, zmax ;

  zvmnmx(a,n,lamin,lamax);

  zmin = *(a + *lamin);
  zmax = *(a + *lamax); 

  *amin = ZABS(zmin);
  *amax = ZABS(zmax);
}

/**<dzvran>**************************************************************

    Title: Complex Vector RANge

  Purpose: Find the minimum and maximum values,
           and their respective indices, of a complex vector.
 
    Input: a = complex vector
           n = length(x)

   Output:  amin (float) = min(a)
           lamin   (int) = location of amin in a[] 
	    amax (float) = max(a)
           lamax   (int) = location of amax in a[] 

    Notes: call must use addresses!
           zvran(a,&amin,&lamin,&amax,&lamax,n);
*/
void dzvran(dcomplex *a, double *amin, int *lamin, double *amax, int *lamax, int n)
{
  dcomplex zmin, zmax ;

  dzvmnmx(a,n,lamin,lamax);

  zmin = *(a + *lamin);
  zmax = *(a + *lamax); 

  *amin = dzabs(zmin);
  *amax = dzabs(zmax);
}

/**<print_zvran.c>****************************************************

     Title: PRINT Complex Vector RANge

   Purpose: Report complex vector min,max and their indeces

      Call: void print_zvran(complex *z,int nsamp, char *text);

     Input: z = complex vector
	    n = length(z)

    Output: none

    History: lrf 11/2005
*/
void print_zvran(complex *z,int nsamp, char *text)
{
  float zmin,zmax;
  int lzmin,lzmax;

  zvran(z,&zmin,&lzmin,&zmax,&lzmax,nsamp);

  fprintf(stdout,"%s:min = %f at index %i\n", text,zmin,lzmin); fflush(stdout);
  fprintf(stdout,"%s:max = %f at index %i\n", text,zmax,lzmax); fflush(stdout);
}

/**<print_dzvran.c>****************************************************

     Title: PRINT Complex Vector RANge

   Purpose: Report complex vector min,max and their indeces

      Call: void print_zvran(complex *z,int nsamp, char *text);

     Input: z = complex vector
	    n = length(z)

    Output: none

    History: lrf 11/2005
*/
void print_dzvran(dcomplex *z,int nsamp, char *text)
{
  double dzmin,dzmax;
  int lzmin,lzmax;

  dzvran(z,&dzmin,&lzmin,&dzmax,&lzmax,nsamp);

  fprintf(stdout,"%s: min = %f at index %i\n",text,dzmin,lzmin); fflush(stdout);
  fprintf(stdout,"%s: max = %f at index %i\n",text,dzmax,lzmax); fflush(stdout);
}
/**<vscal>**************************************************

    Title: Float Vector SCALing

  Purpose:  Scales float vector x[] by a and b giving z[]
            z[i] = a*x[i] + b for i=0,...,n-1 if z<=x
                              for i=n-1,...,0 if z>x

     Call: float *vscal(z,a,x,b,n) 
           float *z,a,*x,b ;	  
	   int n ;			

    Notes: Works correctly regardless of how x[] and z[] might overlap
           Handles cases of a=1,-1 specially
*/

float *vscal(float *z, float a, float *x, float b, int N)
{
	float *zold = z;

	if(N>0)
	{
		if(z<=x)
		{
		  if(a==1.)
		    {
			do {
				*z++ = *x++ + b ;
			} while( --N > 0) ;
		      }
		  else if(a== -1.)
		    {
			do {
				*z++ = - *x++ + b ;
			} while( --N > 0) ;
		      }
		  else
		    {
			do {
				*z++ = *x++ * a + b ;
			} while( --N > 0) ;
		      }
		}
		else
		{
		        z += N ;
			x += N ;
			if(a==1.)
			  {
			    do {
				*--z = *--x + b;
			      } while( --N > 0) ;
			  }
			else if(a== -1.)
			  {
			    do {
				*--z = - *--x + b;
			      } while( --N > 0) ;
			  }
			else
			  {
			    do {
				*--z = *--x * a + b;
			      } while( --N > 0) ;
			  }
		      }
	      }
	return(zold) ;
}

/**<dvscal>**************************************************

    Title: Double Vector SCALing

  Purpose:  Scales double vector x[] by a and b giving z[]
            z[i] = a*x[i] + b for i=0,...,n-1 if z<=x
                              for i=n-1,...,0 if z>x

     Call: double *dvscal(z,a,x,b,n) 
           double *z,a,*x,b ;	  
	   int n ;			

    Notes: Works correctly regardless of how x[] and z[] might overlap
           Handles cases of a=1,-1 specially
*/

double *dvscal(double *z, double a, double *x, double b, int N)
{
	double *zold = z;

	if(N>0)
	{
		if(z<=x)
		{
		  if(a==1.)
		    {
			do {
				*z++ = *x++ + b ;
			} while( --N > 0) ;
		      }
		  else if(a== -1.)
		    {
			do {
				*z++ = - *x++ + b ;
			} while( --N > 0) ;
		      }
		  else
		    {
			do {
				*z++ = *x++ * a + b ;
			} while( --N > 0) ;
		      }
		}
		else
		{
		        z += N ;
			x += N ;
			if(a==1.)
			  {
			    do {
				*--z = *--x + b;
			      } while( --N > 0) ;
			  }
			else if(a== -1.)
			  {
			    do {
				*--z = - *--x + b;
			      } while( --N > 0) ;
			  }
			else
			  {
			    do {
				*--z = *--x * a + b;
			      } while( --N > 0) ;
			  }
		      }
	      }
	return(zold) ;
}

/**<vclip.c>****************************************************

     Title: Vector CLIP

   Purpose: Clips an n point float array x[] to a minimum 
            value of xmin, and a maximum value of xmax, 
	    puts result in z[]

      Call: int vclip( float *z, float *x, int n,
                       float xmin, float xmax)

    Output:  Returns -1 if xmin > xmax, 
             otherwise returns the number of points clipped
 
    History: Musicus, circa 1980
*/
int vclip( float *z, float *x, int n,
	   float xmin, float xmax)
{
	float *zold = z ;
	int count = 0 ;

	if(xmin>xmax)
	{
		if(n>0)
		{
			do {
				*z++ = *x++ ;
			} while( --n > 0) ;
		}
		return(-1) ;
	}
	if(n>0)
	{
		do {
			if(*x < xmin) {
				*z++ = xmin ;
				x++ ;
				count++ ;
			} else if(*x > xmax) {
				*z++ = xmax ;
				x++ ;
				count++ ;
			} else {
				*z++ = *x++ ;
			}
		} while( --n > 0) ;
	}
	return(count) ;
}

/**<vexp.c>****************************************************

     Title: Float Vector EXPonential

   Purpose: Computes z[i] = exp(x[i])  for i=0,...,N-1

      Call: float *vexp(float *z, float *x, int n)

     Input: x = input vector
	    n = length(x)

    Output: Returns pointer to z[] array
*/
float *vexp(float *z, float *x, int N)
{
	float *zptr = z ;

	if(N>0)
	{
		do {
			*z++ = exp(*x++) ;
		} while( --N > 0 ) ;
	}
	return(zptr) ;
}

/**<zvexp.c>****************************************************

     Title: Complex Vector EXPonential

   Purpose: Computes the complex exponential of a complex vector

      Call: complex *zvexp(complex *z, complex *x, int n)

     Input: x = input vector
	    n = length(x)

    Output: z[i] = exp(x[i])

    Returns pointer to z[] array
*/
complex *zvexp(complex *z, complex *x, int n)
{
  float *zptr = (float*)z ;
  double rexp ;

  while(n-->0) {
    rexp = exp(x->r) ;
    *zptr++ = rexp*cos(x->i) ;
    *zptr++ = rexp*sin(x->i) ;
    x++ ;
  }
  return(z) ;
}

/**<zdiv.c>****************************************************

     Title: Complex DIVision

   Purpose: Divide complex args

      Call: complex *zdiv( complex *z, complex *x, complex *y, int N)

    Output: Returns complex quotient a/b 

   History: Musicus, circa 1989

     Notes: Uses special care to avoid unnecessary overflow
*/
complex zdiv(complex a, complex b)
{
	static complex c;
	float  ratio, den;

	if( fabs(b.r) <= fabs(b.i) ) {
		ratio = b.r/b.i;
		den = (ratio*ratio + 1.) * b.i;
		c.r = (a.r*ratio + a.i) / den;
		c.i = (a.i*ratio - a.r) / den;
	}
	else {
		ratio = b.i/b.r;
		den = (ratio*ratio + 1.) * b.r;
		c.r = (a.i*ratio + a.r) / den;
		c.i = (a.i - a.r*ratio) / den;
	}
	return(c);
}

/**<dzdiv.c>****************************************************

     Title: Dcomplex DIVision

   Purpose: Divide dcomplex args

      Call: dcomplex *zdiv( dcomplex *z, dcomplex *x, dcomplex *y, int N)

    Output: Returns dcomplex quotient a/b 

   History: Musicus, circa 1989

     Notes: Uses special care to avoid unnecessary overflow
*/
dcomplex dzdiv(dcomplex a, dcomplex b)
{
	static dcomplex c;
	double  ratio, den;

	if( fabs(b.dr) <= fabs(b.di) ) {
		ratio = b.dr/b.di;
		den = (ratio*ratio + 1.) * b.di;
		c.dr = (a.dr*ratio + a.di) / den;
		c.di = (a.di*ratio - a.dr) / den;
	}
	else {
		ratio = b.di/b.dr;
		den = (ratio*ratio + 1.) * b.dr;
		c.dr = (a.di*ratio + a.dr) / den;
		c.di = (a.di - a.dr*ratio) / den;
	}
	return(c);
}

/**<zvdiv.c>****************************************************

     Title: Complex Vector DIVision

   Purpose: Divide complex vector by another
            z[i] = x[i] / y[i]  for i=0,...,N-1

      Call: complex *zvdiv( complex *z, complex *x, complex *y, int N)

    Output: Returns pointer to z[] array

   History: Musicus, circa 1989
            02_05_04 lrf - make in-place

     Notes: Uses special care to avoid unnecessary overflow
*/
complex *zvdiv( complex *z, complex *x, complex *y, int N)
{
	float *pz = (float*)z ;
	float  ratio, den;
	complex ztmp;

#define mabs(x)	(((x)>=0)?(x):(-(x)))

	if(N>0)
	{
		do {
			if( mabs(y->r) <= mabs(y->i) )
			{
				ratio = y->r/y->i;
				den = (ratio*ratio + 1.) * y->i;
				ztmp.r = (x->r*ratio + x->i) / den; /* real */
				ztmp.i = (x->i*ratio - x->r) / den; /* imag */
				*pz++  = ztmp.r ;
				*pz++  = ztmp.i ;
				*x++ ; *y++ ;
			}
			else
			{
				ratio = y->i/y->r;
				den = (ratio*ratio + 1.) * y->r;
				ztmp.r = (x->i*ratio + x->r) / den; /* real */
				ztmp.i = (x->i - x->r*ratio) / den; /* imag */
				*pz++  = ztmp.r ;
				*pz++  = ztmp.i ;
				*x++ ; *y++ ;
			}
		} while(--N > 0) ;
	}
	return(z);
}

/**<dzvdiv.c>****************************************************

     Title: Double Complex Vector DIVision

   Purpose: Divide dcomplex vector by another
            z[i] = x[i] / y[i]  for i=0,...,N-1

      Call: dcomplex *dzvdiv( dcomplex *z, dcomplex *x, dcomplex *y, int N)

    Output: Returns pointer to z[] array

   History: Musicus, circa 1989
            02_05_04 lrf - make in-place

     Notes: Uses special care to avoid unnecessary overflow
*/
dcomplex *dzvdiv( dcomplex *z, dcomplex *x, dcomplex *y, int N)
{
	double *pz = (double*)z ;
	double  ratio, den;
	dcomplex ztmp;

#define mabs(x)	(((x)>=0)?(x):(-(x)))

	if(N>0)
	{
		do {
			if( mabs(y->dr) <= mabs(y->di) )
			{
				ratio = y->dr/y->di;
				den = (ratio*ratio + 1.) * y->di;
				ztmp.dr = (x->dr*ratio + x->di) / den; /* real */
				ztmp.di = (x->di*ratio - x->dr) / den; /* imag */
				*pz++  = ztmp.dr ;
				*pz++  = ztmp.di ;
				*x++ ; *y++ ;
			}
			else
			{
				ratio = y->di/y->dr;
				den = (ratio*ratio + 1.) * y->dr;
				ztmp.dr = (x->di*ratio + x->dr) / den; /* real */
				ztmp.di = (x->di - x->dr*ratio) / den; /* imag */
				*pz++  = ztmp.dr ;
				*pz++  = ztmp.di ;
				*x++ ; *y++ ;
			}
		} while(--N > 0) ;
	}
	return(z);
}
/*	float *vabs(z,x,n)	Sets z[i] = |x[i]| for i=0,...,n-1
 *	float *x,*z ;		returns pointer to z[]
 *	int n ;
*/
float *vabs(float *z, float *x, int n)
{
	float *zold = z;

	if(n>0)
	{
		do {
			*z++ = (*x>0) ? *x : - *x ;
			x++ ;
		} while( --n > 0) ;
	}
	return(zold) ;
}

/*	double *dvabs(z,x,n)	Sets z[i] = |x[i]| for i=0,...,n-1
 *	double *x,*z ;		returns pointer to z[]
 *	int n ;
*/
double *dvabs(double *z, double *x, int n)
{
	double *zold = z;

	if(n>0)
	{
		do {
			*z++ = (*x>0) ? *x : - *x ;
			x++ ;
		} while( --n > 0) ;
	}
	return(zold) ;
}

/**<zvabs.c>****************************************************

     Title: Complex Vector ABSolute value

   Purpose: Compute magnitude of complex vector
            z[i] = |x[i]|  for i=0,...,N-1

      Call: complex *zvabs(float *z, complex *x,int N)

    Return: Returns pointer to z[]
    
     Notes: Uses special care to avoid overflow
*/
float *zvabs( float *z, complex *x, int N)
{
	double absr, absi, ratio;
	float *xptr = (float*) x ;
	float *oldz = z ;

	if(N>0)
	{
	    do {
		absr = fabs(*xptr++) ;
		absi = fabs(*xptr++) ;
		if( absr <= absi && absi>0.0 ) {
			ratio = absr/absi;
			*z++ = absi * sqrt(1 + ratio*ratio); 
		}
		else if (absr>0.0) {
			ratio = absi/absr;
			*z++ = absr * sqrt(1 + ratio*ratio);
		}
		else {
			*z++ = 0.0;
		}
	    } while(--N > 0) ;
	}
	return(oldz) ;
}

/**<dzvabs.c>****************************************************

     Title: Double Complex Vector ABSolute value

   Purpose: Compute magnitude of complex vector
            z[i] = |x[i]|  for i=0,...,N-1

      Call: complex *dzvabs(double *z, complex *x,int N)

    Return: Returns pointer to z[]
    
     Notes: Uses special care to avoid overflow
*/
double *dzvabs( double *z, dcomplex *x, int N)
{
	double absr, absi, ratio;
	double *xptr = (double*) x ;
	double *oldz = z ;

	if(N>0)
	{
	    do {
		absr = fabs(*xptr++) ;
		absi = fabs(*xptr++) ;
		if( absr <= absi && absi>0.0 ) {
			ratio = absr/absi;
			*z++ = absi * sqrt(1 + ratio*ratio); 
		}
		else if (absr>0.0) {
			ratio = absi/absr;
			*z++ = absr * sqrt(1 + ratio*ratio);
		}
		else {
			*z++ = 0.0;
		}
	    } while(--N > 0) ;
	}
	return(oldz) ;
}

/* complex *zvcplx(z,real,imag,N)	Combines real and imaginary arrays
   complex z[] ;			into complex array
   float real[],imag[] ;		  for(i=0 ; i<N ; i++) {
   int N ;				    z[i].r=real[i] ; z[i].i=imag[i] ; }

   Returns pointer to z[] array
   If real=NULL or imag=NULL, then zvcplx fills in with zeroes
   If both real=NULL and imag=NULL then zvcplx fills the complex vector
   in place using a conquer and divide algorithm (time NlogN instead of N)
   (assumes EQUIVALENCE (z[0],real[0]),((float*)z+N,imag[0]) )
*/
complex *zvcplx( complex *z, float *real, float *imag, int N )
{

#define NULL0 0

	float *pz = (float*)z ;
	float *px ;
	float *py ;

	if(N>0)
	{
		if(real==(float*)NULL0)
		{
			if(imag==(float*)NULL0)
			{
				unsplit((float*)z,N) ;
			}
			else
			{
				px = imag ;
				do {
					*pz++ = 0. ;
					*pz++ = *px++ ;
				} while(--N > 0) ;
			}
		}
		else
		{
			px = real ;
			if(imag==(float*)NULL0)
			{
				do {
					*pz++ = *px++ ;
					*pz++ = 0. ;
				} while(--N > 0) ;
			}
			else
			{
				py = imag ;
				do {
					*pz++ = *px++ ;
					*pz++ = *py++ ;
				} while(--N > 0) ;
			}
		}
	}
	return(z) ;
}

/* dcomplex *dzvcplx(z,real,imag,N)	Combines real and imaginary arrays
   dcomplex z[] ;			into dcomplex array
   double real[],imag[] ;		  for(i=0 ; i<N ; i++) {
   int N ;				    z[i].dr=real[i]; z[i].di=imag[i] ;}

 Returns pointer to z[] array
 If real=NULL or imag=NULL, then dzvcplx fills in with zeroes
 If both real=NULL and imag=NULL then dzvcplx fills the complex vector
   in place using a conquer and divide algorithm (time NlogN instead of N)
   (assumes EQUIVALENCE (z[0],real[0]),((double*)z+N,imag[0]) )
*/
#define NULL0	0

dcomplex *dzvcplx(dcomplex *z, double *real, double *imag, int N)
{
	double *pz = (double*)z ;
	double *px ;
	double *py ;

	if(N>0)
	{
		if(real==(double*)NULL0)
		{
			if(imag==(double*)NULL0)
			{
				dunsplit((double*)z,N) ;
			}
			else
			{
				px = imag ;
				do {
					*pz++ = 0. ;
					*pz++ = *px++ ;
				} while(--N > 0) ;
			}
		}
		else
		{
			px = real ;
			if(imag==(double*)NULL0)
			{
				do {
					*pz++ = *px++ ;
					*pz++ = 0. ;
				} while(--N > 0) ;
			}
			else
			{
				py = imag ;
				do {
					*pz++ = *px++ ;
					*pz++ = *py++ ;
				} while(--N > 0) ;
			}
		}
	}
	return(z) ;
}

/*      float *vmove(z,x,n)	Moves vector x[] into z[] quickly
 *      float *z,*x ;		  z[i] = x[i] for i=0,...,n-1 if z<=x
 *	int n ;				      for i=n-1,...,0 if z>x
 *				Returns pointer to z[]
 *
 * Works correctly regardless of how x[] and z[] might overlap
*/
float *vmove(float *z, float *x,int N)
{
	float *zold = z;

	if(N>0)
	{
		if(z<=x)
		{
			do {
				*z++ = *x++ ;
			} while( --N > 0) ;
		}
		else
		{
			z += N ;
			x += N ;
			do {
				*--z = *--x ;
			} while( --N > 0) ;
		}
	}
	return(zold) ;
}

/*      double *dvmove(z,x,n)	Moves vector x[] into z[] quickly
 *      double *z,*x ;		  z[i] = x[i] for i=0,...,n-1 if z<=x
 *	int n ;				      for i=n-1,...,0 if z>x
 *				Returns pointer to z[]
 *
 * Works correctly regardless of how x[] and z[] might overlap
*/
double *dvmove(double *z, double *x,int N)
{
	double *zold = z;

	if(N>0)
	{
		if(z<=x)
		{
			do {
				*z++ = *x++ ;
			} while( --N > 0) ;
		}
		else
		{
			z += N ;
			x += N ;
			do {
				*--z = *--x ;
			} while( --N > 0) ;
		}
	}
	return(zold) ;
}

/**<vclr.c>****************************************************

     Title: Float Vector CLeaR

   Purpose: Clears float vector: x[i]=0 for i=0,...,n-1

      Call: float *vclr( float *x, int n )

    Output: Returns pointer to x[]
*/
float *vclr( float *x, int n )
{
	float *xold = x ;

	if(n>0)
	{
		do {
			*x++ = 0 ;
		} while( --n>0) ;
	}
	return(xold) ;
}

/**<dvclr.c>****************************************************

     Title: Double Vector CLeaR

   Purpose: Clears double vector: x[i]=0 for i=0,...,n-1

      Call: double *dvclr( double *x, int n )

    Output: Returns pointer to x[]
*/
double *dvclr( double *x, int n )
{
	double *xold = x ;

	if(n>0)
	{
		do {
			*x++ = 0 ;
		} while( --n>0) ;
	}
	return(xold) ;
}

/**<vfill.c>****************************************************

     Title: Float Vector FILL

   Purpose: Fills float vector: x[i]=a for i=0,...,n-1

      Call: float *vfill( float *x, float a, int n )

    Output: Returns pointer to x[]
*/
float *vfill( float *x, float a, int n )
{
	float *xold = x ;

	if(n>0)
	{
		do {
			*x++ = a ;
		} while( --n>0) ;
	}
	return(xold) ;
}

/**<dvfill.c>****************************************************

     Title: Double Vector FILL

   Purpose: Fills double vector: x[i]=a for i=0,...,n-1

      Call: double *vfill( double *x, double a, int n )

    Output: Returns pointer to x[]
*/
double *dvfill( double *x, double a, int n )
{
	double *xold = x ;

	if(n>0)
	{
		do {
			*x++ = a ;
		} while( --n>0) ;
	}
	return(xold) ;
}

/**<zvclr.c>****************************************************

     Title: Complex Vector CLeaR

   Purpose: Clears complex vector: x[i]=0 for i=0,...,N-1

      Call: complex *zvclr(complex *x , int N )

    Output: Returns pointer to x[]
*/
complex *zvclr(complex *x , int N )
{
	float *rx = (float*) x ;

	if(N>0)
	{
		do {
			*rx++ = 0 ;
			*rx++ = 0 ;
		} while(--N > 0) ;
	}
	return(x) ;
}
/**<dzvclr.c>****************************************************

     Title: Dcomplex Vector CLeaR

   Purpose: Clears dcomplex vector: x[i]=0 for i=0,...,N-1

      Call: dcomplex *dzvclr(dcomplex *x , int N )

    Output: Returns pointer to x[]
*/
dcomplex *dzvclr(dcomplex *x , int N )
{
	double *rx = (double*) x ;

	if(N>0)
	{
		do {
			*rx++ = 0 ;
			*rx++ = 0 ;
		} while(--N > 0) ;
	}
	return(x) ;
}

/*	float *vadd(z,x,a,y,N)	Adds scaled float vector to another vector:
 *	float *z,*x,a,*y ;	  z[] = x[] + a*y[]
 *	int N ;			where x,y,z are float arrays of length N, and
 *				a is a double (float) scalar
 *
 * Handles cases of a=1,-1 specially
 * Returns pointer to z[] array
*/
float *vadd(float *z, float *x, float a, float *y, int N)
{
	float *zold = z;

	if(N>0)
	{
		if(a==1.)
		{
			do {
				*z++ = *x++ + *y ;
				y++ ;
			} while(--N > 0) ;
		}
		else if(a== -1.)
		{
			do {
				*z++ = *x++ - *y ;
				y++ ;
			} while(--N > 0) ;
		}
		else
		{
			do {
				*z++ = *x++ + a * *y++ ;
			} while(--N > 0) ;
		}
	}
	return(zold) ;
}

/*	double *dvadd(z,x,a,y,N) Adds scaled double vector to another vector:
 *	double *z,*x,a,*y ;	  z[] = x[] + a*y[]
 *	int N ;			 where x,y,z are double arrays of length N, and
 *				 a is a double scalar
 *
 * Handles cases of a=1,-1 specially
 * Returns pointer to z[] array
 * (Software hacked to improve compiler output)
*/
double *dvadd(double *z, double *x, double a, double *y, int N)
{
	double *zold = z;

	if(N>0)
	{
		if(a==1.)
		{
			do {
				*z++ = *x++ + *y ;
				y++ ;
			} while(--N > 0) ;
		}
		else if(a== -1.)
		{
			do {
				*z++ = *x++ - *y ;
				y++ ;
			} while(--N > 0) ;
		}
		else
		{
			do {
				*z++ = *x++ + *y++ * a ;
			} while(--N > 0) ;
		}
	}
	return(zold) ;
}

/**<vacc.c>****************************************************

     Title: Float Vector Accumulate

   Purpose: z[] += a*y[]

      Call: float *vacc(z,a,y,N) 

    Output: Returns pointer to z[]

  Comments: Handles cases of a=1,-1 specially
            Returns pointer to z[] array
*/
float *vacc(float *z, float a, float *y, int N)
{
	float *zold = z;

	if(N>0)
	{
		if(a==1.)
		{
			do {
				*z++ += *y ;
				y++ ;
			} while(--N > 0) ;
		}
		else if(a== -1.)
		{
			do {
				*z++ -= *y ;
				y++ ;
			} while(--N > 0) ;
		}
		else
		{
			do {
				*z++ += *y++ * a ;
			} while(--N > 0) ;
		}
	}
	return(zold) ;
}
/**<dvacc.c>****************************************************

     Title: Double Vector Accumulate

   Purpose: z[] += a*y[]

      Call: double *dvacc(z,a,y,N) 

    Output: Returns pointer to z[]

  Comments: Handles cases of a=1,-1 specially
            Returns pointer to z[] array
*/
double *dvacc(double *z, double a, double *y, int N)
{
	double *zold = z;

	if(N>0)
	{
		if(a==1.)
		{
			do {
				*z++ += *y ;
				y++ ;
			} while(--N > 0) ;
		}
		else if(a== -1.)
		{
			do {
				*z++ -= *y ;
				y++ ;
			} while(--N > 0) ;
		}
		else
		{
			do {
				*z++ += *y++ * a ;
			} while(--N > 0) ;
		}
	}
	return(zold) ;
}

/**<zvadd.c>****************************************************

     Title: Complex Vector ADDition

   Purpose: Add scaled complex vector to another
            z[i] = x[i] + a*y[i]  for i=0,...,N-1

      Call: complex *zvadd(complex *z, complex *x, complex a, 
                           complex *y, int N)

    Output: Returns pointer to z[] array

     Notes: For extra speed, real, imaginary, and +1, -1 scalars handled specially

   History: Musicus, circa 1989

*/
complex *zvadd(complex *z, complex *x, complex a, complex *y, int N)
{
	float *rz = (float*) z ;
	float *rx = (float*) x ;
	float *ry = (float*) y ;

	if(N>0)
	{
		if(a.i == 0)
		{
			if(a.r == 1)
			{
				do {
					*rz++ = *rx++ + *ry++ ;
					*rz++ = *rx++ + *ry++ ;
				} while(--N>0) ;
			}
			else if(a.r == -1)
			{
				do {
					*rz++ = *rx++ - *ry++ ;
					*rz++ = *rx++ - *ry++ ;
				} while(--N>0) ;
			}
			else
			{
				do {
					*rz++ = *rx++ + a.r * *ry++ ;
					*rz++ = *rx++ + a.r * *ry++ ;
				} while(--N>0) ;
			}
		}
		else if(a.r == 0)
		{
			if(a.i == 1)
			{
				do {
					*rz++ = *rx++ - y->i ;
					*rz++ = *rx++ + (y++)->r ;
				} while(--N>0) ;
			}
			else if(a.i == -1)
			{
				do {
					*rz++ = *rx++ + y->i ;
					*rz++ = *rx++ - (y++)->r ;
				} while(--N>0) ;
			}
			else
			{
				do {
					*rz++ = *rx++ - a.i * y->i ;
					*rz++ = *rx++ + a.i * (y++)->r ;
				} while(--N>0) ;
			}
		}
		else
		{
			do {
				*rz++ = *rx++ + a.r * y->r - a.i * y->i ;
				*rz++ = *rx++ + a.r * y->i + a.i * y->r ;
				y++ ;
			} while(--N>0) ;
		}
	}
	return(z) ;
}
/**<dzvadd.c>****************************************************

     Title: Double Complex Vector ADDition

   Purpose: Add scaled complex vector to another
            z[i] = x[i] + a*y[i]  for i=0,...,N-1

      Call: dcomplex *dzvadd(dcomplex *z, dcomplex *x, dcomplex a, dcomplex *y, int N)

    Output: Returns pointer to z[] array

     Notes: For extra speed, real, imaginary, and +1, -1 scalars handled specially

   History: Musicus, circa 1989

*/
dcomplex *dzvadd(dcomplex *z, dcomplex *x, dcomplex a, dcomplex *y, int N)
{
	double *rz = (double*) z ;
	double *rx = (double*) x ;
	double *ry = (double*) y ;

	if(N>0)
	{
		if(a.di == 0)
		{
			if(a.dr == 1)
			{
				do {
					*rz++ = *rx++ + *ry++ ;
					*rz++ = *rx++ + *ry++ ;
				} while(--N>0) ;
			}
			else if(a.dr == -1)
			{
				do {
					*rz++ = *rx++ - *ry++ ;
					*rz++ = *rx++ - *ry++ ;
				} while(--N>0) ;
			}
			else
			{
				do {
					*rz++ = *rx++ + a.dr * *ry++ ;
					*rz++ = *rx++ + a.dr * *ry++ ;
				} while(--N>0) ;
			}
		}
		else if(a.dr == 0)
		{
			if(a.di == 1)
			{
				do {
					*rz++ = *rx++ - y->di ;
					*rz++ = *rx++ + (y++)->dr ;
				} while(--N>0) ;
			}
			else if(a.di == -1)
			{
				do {
					*rz++ = *rx++ + y->di ;
					*rz++ = *rx++ - (y++)->dr ;
				} while(--N>0) ;
			}
			else
			{
				do {
					*rz++ = *rx++ - a.di * y->di ;
					*rz++ = *rx++ + a.di * (y++)->dr ;
				} while(--N>0) ;
			}
		}
		else
		{
			do {
				*rz++ = *rx++ + (a.dr * y->dr - a.di * y->di) ;
				*rz++ = *rx++ + (a.dr * y->di + a.di * y->dr) ;
				y++ ;
			} while(--N>0) ;
		}
	}
	return(z) ;
}

/*	complex zvsum(x,n)	Sums elements of vector x[0],...,x[n-1]
 *	complex *x ;		Internal calculation is double precision,
 *	int n ;
*/
complex zvsum(complex *x, int n)
{
	complex sum;
	double sumr = 0;
	double sumi = 0;
	float *px = (float*)x;

	if(n>0)
	{
		do {
			sumr += *px++ ;
			sumi += *px++ ;
		} while(--n > 0) ;
	}
	sum.r = sumr ;
	sum.i = sumi ;
	return(sum) ;
}

/**<zvmean>**************************************************

    Title: Complex Vector MEAN

  Purpose: Compute mean of complex vector.
 
     Call: complex zvmean( complex *x, int n)

*/
complex zvmean( complex *x, int n)
{
  complex sum, mean;

  sum = zvsum(x,n);
  mean.r = sum.r/n;
  mean.i = sum.i/n;

  return(mean);
}

/**<zvstd>**************************************************

    Title: Complex Vector STD

  Purpose: Compute std of complex vector.
 
     Call: double zvstd( complex *x, int n )

     Note: Returns double (_not_ complex)!
*/
double zvstd( complex *x, int n )
{
       complex  mean, temp, sum ;
       double   std;
       double   dn = n-1 ;
       float    *dx = (float*) x;

       sum.r = 0.;
       sum.i = 0.;
       mean = zvmean(x,n);

       if(n>0)
	 {
	   do {
	     temp.r = *dx++ - mean.r ;
	     temp.i = *dx++ - mean.i ;
	     sum.r += temp.r*temp.r ;
	     sum.i += temp.i*temp.i ;
	   } while(--n > 0) ;
	 }
  
       sum.r = sqrt(sum.r / dn);
       sum.i = sqrt(sum.i / dn);

       std = zabs(sum);

       return(std);
}
/**<zvscal.c>****************************************************

     Title: Complex Vector SCALe

   Purpose: Scales complex vector x[] by a giving z[]

            z[i] = a*x[i] for i=0,...,n-1 if z<=x
	                  for i=n-1,...,0 if z>x

      Call: complex *zvscal(complex *z, complex a, complex *x, int N)

    Output: Returns pointer to x[]

     Notes: Works correctly regardless of how x[] and z[] might overlap
            Twice as fast if a is real or pure imaginary
*/
complex *zvscal(complex *z, complex a, complex *x, int N)
{
	float *pz,*px ;

	if(N>0)
	{
		if(z <= x)
		{
			pz = (float*)z ;
			if(a.i == 0)
			{
				px = (float*)x ;
				do {
					*pz++ = a.r * *px++ ;
					*pz++ = a.r * *px++ ;
				} while( --N > 0) ;
			}
			else if(a.r == 0)
			{
				do {
					*pz++ = -a.i * x->i ;
					*pz++ = a.i * x->r ;
					x++ ;
				} while( --N > 0) ;
			}
			else
			{
				do {
					*pz++ = a.r * x->r - a.i * x->i ;
					*pz++ = a.r * x->i + a.i * x->r ;
					x++ ;
				} while( --N > 0) ;
			}
		}
		else
		{
			pz = (float*)(z+N) ;
			x += N ;
			if(a.i == 0)
			{
				px = (float*)x ;
				do {
					*--pz = a.r * *--px ;
					*--pz = a.r * *--px ;
				} while( --N > 0) ;
			}
			else if(a.r == 0)
			{
				do {
					--x ;
					*--pz = a.i * x->r ;
					*--pz = -a.i * x->i ;
				} while( --N > 0) ;
			}
			else
			{
				do {
					--x ;
					*--pz = a.r*x->i + a.i*x->r ;
					*--pz = a.r*x->r - a.i*x->i ;
				} while( --N > 0) ;
			}
		}
	}
	return(z) ;
}

/**<dzvscal.c>****************************************************

     Title: Double Complex Vector SCALe

   Purpose: Scales dcomplex vector x[] by a giving z[]

            z[i] = a*x[i] for i=0,...,n-1 if z<=x
	                  for i=n-1,...,0 if z>x

      Call: dcomplex *dzvscal(dcomplex *z, dcomplex a, dcomplex *x, int N)

    Output: Returns pointer to x[]

     Notes: Works correctly regardless of how x[] and z[] might overlap
            Twice as fast if a is real or pure imaginary
*/
dcomplex *dzvscal(dcomplex *z, dcomplex a, dcomplex *x, int N)
{
	double *pz,*px ;

	if(N>0)
	{
		if(z <= x)
		{
			pz = (double*)z ;
			if(a.di == 0)
			{
				px = (double*)x ;
				do {
					*pz++ = a.dr * *px++ ;
					*pz++ = a.dr * *px++ ;
				} while( --N > 0) ;
			}
			else if(a.dr == 0)
			{
				do {
					*pz++ = -a.di * x->di ;
					*pz++ = a.di * x->dr ;
					x++ ;
				} while( --N > 0) ;
			}
			else
			{
				do {
					*pz++ = a.dr * x->dr - a.di * x->di ;
					*pz++ = a.dr * x->di + a.di * x->dr ;
					x++ ;
				} while( --N > 0) ;
			}
		}
		else
		{
			pz = (double*)(z+N) ;
			x += N ;
			if(a.di == 0)
			{
				px = (double*)x ;
				do {
					*--pz = a.dr * *--px ;
					*--pz = a.dr * *--px ;
				} while( --N > 0) ;
			}
			else if(a.dr == 0)
			{
				do {
					--x ;
					*--pz = a.di * x->dr ;
					*--pz = -a.di * x->di ;
				} while( --N > 0) ;
			}
			else
			{
				do {
					--x ;
					*--pz = a.dr*x->di + a.di*x->dr ;
					*--pz = a.dr*x->dr - a.di*x->di ;
				} while( --N > 0) ;
			}
		}
	}
	return(z) ;
}

/**<zvzap.c>************************************************************

     Title: Float Complex Vector ZAP in 1 dimension

   Purpose: Nyquist modulation in 1D

     Input: complex z[n]

*/
complex *zvzap( complex *z , int n )
{
  int i ;
  int i0 = 1;

  for (i=i0; i<n; i+=2) {
    z[i].r *= -1. ;
    z[i].i *= -1. ;
  }

  return (z) ;
}

/**<dzvzap.c>************************************************************

     Title: Double Complex Vector ZAP in 1 dimension

   Purpose: Nyquist modulation in 1D

     Input: dcomplex z[n]

*/
dcomplex *dzvzap( dcomplex *z , int n )
{
  int i ;
  int i0 = 1;

  for (i=i0; i<n; i+=2) {
    z[i].dr *= -1. ;
    z[i].di *= -1. ;
  }

  return (z) ;
}

/**<zvzap2rc.c>************************************************************

     Title: Float Complex Vector ZAP in along 1 dimension,
            either Row or Column, in a 1-dim implementation of a 2D-array

   Purpose: 

     Input: complex z[nr*nc]
            int nr = row dim
            int nc = col dim
	    int rowcol: 0 = zap rows
	                1 = zap cols
*/
complex *zvzap2rc( complex *z, int nr, int nc, int rowcol )
{
  int i, j, jj;
  int ixy;

  if(rowcol==0) {		/* zap rows */

    for (j = 0; j < nr; j++) {
      for (i = 0; i < nc; i+=2) {
	ixy =  i +  j * nc ;
	z[ixy].r *= -1. ;
	z[ixy].i *= -1. ;
      }
    }

  } else if (rowcol==1) {	/* zap columns */

    for (j = 0; j < nr; j+=2) {
      for (i = 0; i < nc; i++) {
	ixy = i +  j * nc ;
	z[ixy].r *= -1. ;
	z[ixy].i *= -1. ;
      }
    }

  }

  return (z) ;
}

/**<zvzap2.c>************************************************************

     Title: Float Complex Vector ZAP in 2 dimensions

   Purpose: Nyquist modulation in 2D

     Input: complex z[nr*nc]

*/
complex *zvzap2( complex *z, int nr, int nc )
{
  int i,  j ;
  int i0, j0 = 0;
  int idx;

  int debug = 0;

  for (j = j0; j < nr; j++) {
    if(debug) if(j>4) exit(0);
    i0 = 1 - j % 2 ;
    if(debug) printf("\n");
    for (i = i0; i < nc; i+=2) {
      if(debug) printf("+++++++++++++++++++ (%i,%i) = -1\n",i,j);
      idx = i + j*nc ;
      z[idx].r *= -1. ;
      z[idx].i *= -1. ;
    }
  }

  return (z) ;
}

/**<dzvzap2.c>************************************************************

     Title: Double Complex Vector ZAP in 2 dimensions

   Purpose: Nyquist modulation in 2D

     Input: dcomplex z[nr*nc]

*/
dcomplex *dzvzap2( dcomplex *z, int nx, int ny )
{
  int i,j,jj ;
  int i0, j0 = 0;
  int idx;

  for (j = j0; j < ny; j++) {
    jj = j*nx;
    i0 = 1 - j % 2 ;
    for (i = i0; i < nx; i+=2) {
      idx = i + jj ;
      z[idx].dr *= -1. ;
      z[idx].di *= -1. ;
    }
  }

  return (z) ;
}

/**<dvzap2.c>************************************************************

     Title: Double Vector ZAP in 2 dimensions

   Purpose: Nyquist modulation in 2D

     Input: double z[nr*nc]

*/
double *dvzap2( double *z, int nx, int ny )
{
  int i,j,jj ;
  int i0, j0 = 0;
  int idx;

  for (j = j0; j < ny; j++) {
    jj = j*nx;
    i0 = 1 - j % 2 ;
    for (i = i0; i < nx; i+=2) {
      idx = i + jj ;
      z[idx] *= -1. ;
    }
  }

  return (z) ;
}

/**<zvzap3.c>************************************************************

     Title: Float Complex Vector ZAP in 3 dimensions

   Purpose: Nyquist modulation in 3D

     Input: complex z[nx*ny*nz]

*/
complex *zvzap3_old( complex *z, int nx, int ny, int nz )
{
  int i,j,k ;
  int i0,j0,k0;
  int idx;

  for (k = 0; k < nz; k++) {
    j0 = k % 2 ;
    for (j = j0; j < ny; j++) {
      i0 = (j + k) % 2 ;
      for (i = i0; i < nx; i+=2) {
	idx = k + (j + i*ny)*nz;
	z[idx].r *= -1. ;
	z[idx].i *= -1. ;
      }
    }
  }

  return (z) ;
}

/**<zvzap3.c>************************************************************

     Title: Complex Vector ZAP in 3 dimensions

   Purpose: Nyquist modulation in 3D

     Input: complex z[nx*ny*nz]

*/
complex *zvzap3( complex *z, int nx, int ny, int nz )
{
  int i,k, kk ;
  int nxy = nx*ny;
  complex a = zcplx(-1.,-1.);

  for (k = 0; k < nz; k++) {
    kk = k*nxy;
    zvzap2(&z[kk],nx,ny);
    if(k%2) for (i = 0; i < nxy; i++) z[kk+i] = zscal(-1.,z[kk+i]);
  }

  return (z) ;
}

/**<dvzap3.c>************************************************************

     Title: Double Vector ZAP in 3 dimensions

   Purpose: Nyquist modulation in 3D

     Input: double z[nx*ny*nz]

*/
double *dvzap3( double *z, int nx, int ny, int nz )
{
  int i,k, kk ;
  int nxy = nx*ny;
  double a = -1.;

  for (k = 0; k < nz; k++) {
    kk = k*nxy;
    dvzap2(&z[kk],nx,ny);
    if(k%2) for (i = 0; i < nxy; i++) z[kk+i] *= -1.;
  }

  return (z) ;
}

/**<dzvzap3.c>************************************************************

     Title: Double Complex Vector ZAP in 3 dimensions

   Purpose: Nyquist modulation in 3D

     Input: dcomplex z[nx*ny*nz]

*/
dcomplex *dzvzap3( dcomplex *z, int nx, int ny, int nz )
{
  int i,k, kk ;
  int nxy = nx*ny;
  dcomplex a = dzcplx(-1.,-1.);

  for (k = 0; k < nz; k++) {
    kk = k*nxy;
    dzvzap2(&z[kk],nx,ny);
    if(k%2) for (i = 0; i < nxy; i++) z[kk+i] = dzscal(-1.,z[kk+i]);
  }

  return (z) ;
}

/**<dzvzap3_old.c>************************************************************

     Title: Double Complex Vector ZAP in 3 dimensions

   Purpose: Nyquist modulation in 3D

     Input: dcomplex z[nx*ny*nz]

*/
dcomplex *dzvzap3_old( dcomplex *z, int nx, int ny, int nz )
{
  int i,j,k,jj ;
  int i0,j0,k0;
  int idx;

  int nyz = ny*nz;

  for (k = 0; k < nz; k++) {
    j0 = k % 2 ;
    for (j = j0; j < ny; j++) {
      i0 = (j + k) % 2 ;
      jj = j*nz;
      for (i = i0; i < nx; i+=2) {
	idx = k + jj + i*nyz;
	z[idx].dr *= -1. ;
	z[idx].di *= -1. ;
      }
    }
  }

  return (z) ;
}
