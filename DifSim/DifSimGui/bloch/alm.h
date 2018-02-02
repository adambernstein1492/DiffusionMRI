/*<alm.h>***************************************************************/

#ifndef _LRF_ALM_HEADER_
#define _LRF_ALM_HEADER_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/* lrf includes */
#include "sphfun.h"
#include "alloc.h"
#include "vector.h"
#include "matrix.h"
#include "ioutils.h"
#include "coords.h"
#include "grads.h"

#ifndef TYPEDEF_complex
#define TYPEDEF_complex
typedef struct { float r, i ; } complex ;
#endif

#ifndef TYPEDEF_dcomplex
#define TYPEDEF_dcomplex
typedef struct { double r, i ; } dcomplex ;
#endif

/**************** Alm routines ****************************/

int fac2(int l);
void alm2olm(complex **Alm, complex **Olm, int Lmax);
float **alm2shape( float **shape, complex **Alm, 
		   int Lmin, int Lmax,
		   float *az, float *po, int naz, int npo ) ;

float   **rlm_alloc(int Lmax) ;
complex **alm_alloc(int Lmax) ;
void alm_free(complex **alm, int Lmax) ;
float   **alm2tensor( float **Dap, complex **Alm, int Lmax );
float    *alm2egy( float *egyL, complex **Alm, int Lmax ) ;
complex  *sumalm( complex *sumL, complex **Alm, int Lmax ) ;
complex  pickalm( complex **Alm, int Lpick, int Mpick, int Lmax ) ;
complex **zvtoAlm( complex **Alm, complex *zv, int Lmax );
complex *Almtozv( complex *zv, complex **Alm, int Lmax ) ;
float alm2gra( complex **Alm, int Lmax ) ;
float alm2gfa( complex **Alm, int Lmax ) ;
float alm2fai( complex **Alm, int Lmax ) ;

/* 10-Aug-03 */
int numAlm( int Lmax ) ;
complex almcorr( complex **Alm1, complex **Alm2, int Lmax ) ;
complex **loadalm( complex **AlmTo, complex **AlmFrom, int Lmax ) ;
#define LOAD_ALM(almto,almfrom,lmax) loadalm((almto),(almfrom),(lmax)) ;
complex ***almr_realloc(complex ***almr, int Lmax, int npts) ;

#define CHECK_L(lval,lmax) (((lval)<=(lmax)) ? (1) : (0))
#define CHECK_M(mval,lmax) (((ABS(mval))<=(lmax)) ? (1) : (0))

#define getMrowFromAlm(mrow,Alm,Lval) \
do{ int j; \
      for(j=-(Lval); j<=(Lval); j++) { \
	mrow[j+(Lval)].r = (Alm)[(Lval)][j].r; \
	mrow[j+(Lval)].i = (Alm)[(Lval)][j].i; \
      } \
} while(0) ;

#define printAlm(Alm,Lmax) \
do{ int i,j,nl = (Lmax)+1; \
    for(i=0; i<nl; i++) { \
      printf("\n"); \
      for(j=-i; j<=i; j++) { \
	printf("Alm[%i][%2i] = (%13.6g,%13.6g)\n",\
                i,j,(Alm)[i][j].r,(Alm)[i][j].i); \
      } \
    }\
} while(0) ;

#define fillAlm(Alm,Lmax) \
do{ int i,j,nl = (Lmax)+1; \
    for(i=0; i<nl; i++) { \
      for(j=-i; j<=i; j++) { \
	(Alm)[i][j].r = i; \
        (Alm)[i][j].i = j; \
      } \
    }\
} while(0) ;

#define writeAlm(fp,Alm,Lmax,nb) \
do{ int i,j,nl = Lmax+1; \
   fprintf(fp,"%d %d\n",nb,nb); \
   fprintf(fp,"%d %d\n",Lmax,Lmax); \
   for(i=0; i<nl; i++) { \
     for(j=-i; j<=i; j++) { \
       if(isnan(Alm[i][j].r)||isnan(Alm[i][j].i)) { \
	 printf("writeAlm: Zeroing nan found at (%i,%i)\n",i,j); \
	 fprintf(fp,"%13.6g %13.6g\n",0.,0.) ;\
       } else {\
	 fprintf(fp,"%13.6g %13.6g\n",\
		 Alm[i][j].r,Alm[i][j].i); \
       }\
     }\
   }\
} while(0) ;

/* 20-Sep-03 */
#define scaleAlm(Alm,scal,Lmax) \
do{ int i,j,nl = (Lmax)+1; \
    for(i=0; i<nl; i++) { \
      for(j=-i; j<=i; j++) { \
	(Alm)[i][j].r *= (scal); \
        (Alm)[i][j].i *= (scal); \
      } \
    }\
} while(0) ;

#define printRlm(Rlm,Lmax) \
do{ int i,j,nl = (Lmax)+1; \
    for(i=0; i<nl; i++) { \
      printf("\n"); \
      for(j=-i; j<=i; j++) { \
	printf("Rlm[%i][%2i] = %13.6g\n",\
                i,j,(Rlm)[i][j]); \
      } \
    }\
} while(0) ;

#define fillRlm(Rlm,Lmax,val) \
do{ int i,j,nl = (Lmax)+1; \
    for(i=0; i<nl; i++) { \
      for(j=-i; j<=i; j++) { \
	(Rlm)[i][j].r = (val); \
      } \
    }\
} while(0) ;

#define clearAlm(Alm,Lmax) fillAlm((Alm),(Lmax),(0.))
#define unitAlm(Alm,Lmax)  fillAlm((Alm),(Lmax),(1.))

#define clearRlm(Rlm,Lmax) fillRlm((Rlm),(Lmax),(0.))
#define unitRlm(Rlm,Lmax)  fillRlm((Rlm),(Lmax),(1.))


/**************** Ylm routines ****************************/

complex ***ylm_alloc(int Lmax, int npts) ;
void ylm_free(complex ***ylm, int Lmax) ;

#define printYlm(Ylm,Lmax,nn) \
do{ int i,j,k,nl = Lmax+1; \
    for(i=0; i<nl; i++) { \
     printf("\n"); \
     for(j=-i; j<=i; j++) { \
      for(k=0; k<nn; k++) { \
	printf("Ylm[%i][%2i][%i] = (%13.6g,%13.6g)\n",\
                i,j,k,Ylm[i][j][k].r,Ylm[i][j][k].i); \
      } \
     }\
    }\
} while(0) ;

#define fillYlm(Ylm,Lmax,nn) \
do{   int i,j,k,nl = Lmax+1; \
    for(k=0; k<nn; k++) { \
    for(i=0; i<nl; i++) { \
      for(j=-i; j<=i; j++) { \
	Ylm[i][j][k].r = i*k; \
	Ylm[i][j][k].i = j*k; \
      } \
    }\
    }\
} while(0) ;

/**************** Wlmr routines ****************************/

complex ****wlmr_alloc(int Lmax, int nrad, int nvert) ;

#define printWlmr(Wlmr,Lmax,nv,nr) \
do{ int i,j,k,l,nl = Lmax+1; \
    for(k=0; k<nr; k++) { \
     for(i=0; i<nl; i++) { \
      printf("\n"); \
      for(j=-i; j<=i; j++) { \
       for(l=0; l<nv; l++) { \
	printf("Wlmr[%i][%2i][%i][%i] = (%13.6g,%13.6g)\n",\
                k,i,j,l,Wlmr[k][i][j][l].r,Wlmr[k][i][j][l].i); \
       } \
      } \
     }\
    }\
} while(0) ;


/**************** Routines to replace Almr routines **************************/

complex ***almN_alloc(int Lmax, int npts) ;
complex ****almNM_alloc(int Lmax, int npts, int mpts) ;

/**************** Almr routines ****************************/

/* complex ***almr_alloc(int Lmax, int npts) ;*//*Obsoleted by almN_alloc*/
float   ***rlmr_alloc(int Lmax, int npts);

#define printAlmN(Almr,Lmax,nr) \
do{ int i,j,k,nl = Lmax+1; \
    for(k=0; k<nr; k++) { \
      for(i=0; i<nl; i++) { \
      printf("\n"); \
        for(j=-i; j<=i; j++) { \
	  printf("Almr[%i][%2i][%i] = (%13.6g,%13.6g)\n",\
                  k,i,j,Almr[k][i][j].r,Almr[k][i][j].i); \
        } \
      } \
    }\
} while(0) ;

#define writeAlmr(fp,Almr,Lmax,nb,nr) \
do{ int i,j,k,nl = Lmax+1,count=0; \
   fprintf(fp,"%d %d\n",nb,nb); \
   fprintf(fp,"%d %d\n",nr,nr); \
   fprintf(fp,"%d %d\n",Lmax,Lmax); \
   for(k=0; k<nr; k++) { \
     for(i=0; i<nl; i++) { \
       for(j=-i; j<=i; j++) { \
         if(isnan(Almr[k][i][j].r)||isnan(Almr[k][i][j].i)) { \
	   printf("writeAlmr: Zeroing nan found at (%i,%i,%i)\n",i,j,k); \
	   fprintf(fp,"%13.6g %13.6g\n",0.,0.) ;\
         } else {\
	   fprintf(fp,"%13.6g %13.6g\n",\
		   Almr[k][i][j].r,Almr[k][i][j].i); \
         }\
       }\
     }\
   }\
} while(0) ;

#endif
