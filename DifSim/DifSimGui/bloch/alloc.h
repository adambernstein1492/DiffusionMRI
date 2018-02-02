/*<alloc.h>***************************************************************/

#ifndef _ALLOC_HEADER_
#define _ALLOC_HEADER_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <err.h>
#include <errno.h>

#ifndef TYPEDEF_complex
#define TYPEDEF_complex
typedef struct { float r, i ; } complex ;
#endif

#ifndef TYPEDEF_dcomplex
#define TYPEDEF_dcomplex
typedef struct { double dr, di ; } dcomplex ;
#endif

#ifndef TYPEDEF_byte
#define TYPEDEF_byte
typedef unsigned char byte ;
#endif

void f3_free( float ***ft, int ns, int nr, int nc ) ;
void zt_free( complex ***zt, int ns, int nr, int nc );
void z4_free( complex ****z4, int nm, int ns, int nr, int nc ) ;
void z5_free( complex *****z5, int nq, int nm, int ns, int nr, int nc ) ;
void dzt_free( dcomplex ***dzt, int ns, int nr, int nc );

#define ft_free f3_free

void dz4_free( dcomplex ****dz4, int nm, int ns, int nr, int nc ) ;
void dz5_free( dcomplex *****dz5, int nq, int nm, int ns, int nr, int nc ) ;

void fm_free( float **fm, int nrow, int ncol ) ;
void dm_free( double **dm, int nrow, int ncol ) ;
void im_free( int **im, int nrow, int ncol ) ;
void sm_free( short **sm, int nrow, int ncol ) ;
void zm_free( complex **zm, int nrow, int ncol ) ;
void dzm_free( dcomplex **dzm, int nrow, int ncol ) ;

void _mm_free( char **mm, int nr, int nc ) ;
void _tt_free( char ***tt, int ns, int nr, int nc ) ;
void _tupleN_free( char ****xx, int ns, int nr, int nc, int ncomp ) ;

#define mm_free(mm,nr,nc)    _mm_free((char **)(mm),(nr),(nc)) 
#define tt_free(tt,ns,nr,nc) _tt_free((char ***)(tt),(ns),(nr),(nc))
#define tupleN_free(xx,ns,nr,nc,ncomp) _tupleN_free((char ****)(xx),(ns),(nr),(nc),(ncomp)) 
#define tuple3_free(xx,ns,nr,nc) _tupleN_free((char ****)(xx),(ns),(nr),(nc),(3)) 

char    **cm_alloc(int nr, int nc) ;

short   **sm_alloc( int nr, int nc );
short  ***st_alloc( int ns, int nr, int nc ) ;

int     **im_alloc( int nr, int nc );
int    ***it_alloc( int ns, int nr, int nc ) ;

float   **fm_alloc(int nr, int nc) ;
float  ***f3_alloc( int ns, int nr, int nc ) ;
float ****ftupleN_alloc( int ns, int nr, int nc, int ncomp ) ;
#define ft_alloc f3_alloc

double   **dm_alloc(int nr, int nc) ;
double ****dtupleN_alloc( int ns, int nr, int nc, int ncomp ) ;
double   ***d3_alloc( int ns, int nr, int nc ) ;
double  ****d4_alloc( int nm, int ns, int nr, int nc ) ;
double *****d5_alloc( int nq, int nm, int ns, int nr, int nc ) ;
#define dt_alloc d3_alloc

#define d2_alloc dm_alloc

void dm_free( double **dm, int nrow, int ncol ) ;
void d3_free( double ***dt, int ns, int nr, int nc ) ;
void d4_free( double ****d4, int nm, int ns, int nr, int nc ) ;
void d5_free( double *****d5, int nq, int nm, int ns, int nr, int nc ) ;

#define dt_free d3_free

complex   **zm_alloc(int nr, int nc) ;
complex  ***zt_alloc( int ns, int nr, int nc ) ;
complex ****z4_alloc( int nm, int ns, int nr, int nc ) ;
complex *****z5_alloc( int nq, int nm, int ns, int nr, int nc ) ;
complex ****ztupleN_alloc( int ns, int nr, int nc, int ncomp ) ;

dcomplex **dzm_alloc(int nr, int nc) ;
dcomplex  ***dzt_alloc( int ns, int nr, int nc ) ;
dcomplex ****dz4_alloc( int nm, int ns, int nr, int nc ) ;
dcomplex *****dz5_alloc( int nq, int nm, int ns, int nr, int nc ) ;

#define z3_alloc zt_alloc
#define z3_free zt_free
#define dz2_alloc dzm_alloc
#define dz3_alloc dzt_alloc
#define dz3_free dzt_free
#define dz2_free dzm_free

#define ftuple3_alloc(ns,nr,nc) (ftupleN_alloc(ns,nr,nc,(3)))
#define dtuple3_alloc(ns,nr,nc) (dtupleN_alloc(ns,nr,nc,(3)))

char *cv_alloc(int n) ;
short *sv_alloc(int n) ;
byte *bv_alloc(int n) ;
int *iv_alloc(int n) ;
float *fv_alloc(int n) ;
double *dv_alloc(int n) ;
complex *zv_alloc(int n) ;
dcomplex *dzv_alloc(int n) ;

#define dzv_realloc(v,n) ((dcomplex *) realloc( (void *)(v) , (n)*sizeof(dcomplex)))
#define dwv_alloc(n) (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(n));

#endif /* _ALLOC_HEADER_ */
