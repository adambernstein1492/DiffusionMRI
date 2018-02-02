/* matrix.h */

#ifndef _LRF_MATRIX_HEADER_
#define _LRF_MATRIX_HEADER_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/*#include "mrilib.h"*/
#include "alloc.h"

#ifndef TYPEDEF_byte
#define TYPEDEF_byte
typedef unsigned char byte ;
#endif

#ifndef TYPEDEF_complex
#define TYPEDEF_complex
typedef struct { float r, i ; } complex ;
#endif

double **dmadd(double **z,double **x,double a,double **y,int nr,int nc) ;
double *dvmmul(double *d,double **f, double *a,int m,int n) ;

float **fmmul( float **z, float **x, float **y,
	       int nrx, int ncx, int ncy) ;

float **fmtrn( float **z, float**x, int nr, int nc);

#define DMVMV33(x,b,y) dmvmv((x),(b),(y),3,3)
#define FMMUL33(z,x,y) fmmul((z),(x),(y),3,3,3)
#define ROT33(z,x,y) FMMUL33((z),(x),(y))

#define LOAD_MAT33_DIAG(A,a11,a22,a33) \
        LOAD_MAT33(A,a11,0.0,0.0,0.0,a22,0.0,0.0,0.0,a33) 

#define LOAD_MAT33(A,a11,a12,a13,a21,a22,a23,a31,a32,a33) \
 ( (A)[0][0] = (a11) , (A)[0][1] = (a12) ,      \
   (A)[0][2] = (a13) , (A)[1][0] = (a21) ,      \
   (A)[1][1] = (a22) , (A)[1][2] = (a23) ,      \
   (A)[2][0] = (a31) , (A)[2][1] = (a32) , (A)[2][2] = (a33) )

#define PRINT_MAT33(str,A)                                  \
   fprintf(stderr,                                         \
          "%10.10s: [ %13.6g %13.6g %13.6g ]\n"            \
       "            [ %13.6g %13.6g %13.6g ]\n"            \
       "            [ %13.6g %13.6g %13.6g ]\n" ,          \
     str , (A)[0][0] , (A)[0][1] , (A)[0][2] , \
           (A)[1][0] , (A)[1][1] , (A)[1][2] , \
           (A)[2][0] , (A)[2][1] , (A)[2][2]  )

#define PRINT_FMAT(v,n,m,label) \
do{ int i; \
    for(i=0; i<(n); i++) { \
      fprintf(stderr,"%s[%i] = (",(label),i);	\
      for(j=0; j<(m); j++) { \
        fprintf(stderr,"%7.4g",(v)[i][j]);	\
        if(j<(m-1)) fprintf(stderr,", ");		\
      }\
      fprintf(stderr,")\n");	\
    }\
} while(0) 

#define PRINT_FM3(v,n,label) PRINT_FMAT((v),(n),3,(label))

/*****************************************/

#define COPY_MAT(A,B)                                       \
 ( (A)[0][0] = (B)[0][0] , \
   (A)[1][0] = (B)[1][0] , \
   (A)[2][0] = (B)[2][0] , \
   (A)[0][1] = (B)[0][1] , \
   (A)[1][1] = (B)[1][1] , \
   (A)[2][1] = (B)[2][1] , \
   (A)[0][2] = (B)[0][2] , \
   (A)[1][2] = (B)[1][2] , \
   (A)[2][2] = (B)[2][2] )

#define LOAD_ROTGEN_MAT(A,th,ff,aa,bb)             \
 ( (A)[aa][aa] = (A)[bb][bb] = cos((th)) , \
   (A)[aa][bb] = -sin((th)) ,                   \
   (A)[bb][aa] = -(A)[aa][bb] ,            \
   (A)[ff][ff] = 1.0 ,                         \
   (A)[aa][ff] = (A)[bb][ff] = (A)[ff][aa] = (A)[ff][bb] = 0.0 )

#define LOAD_DIAG_MAT(A,x,y,z) \
 ( (A)[0][0] = (x) , \
   (A)[1][1] = (y) , \
   (A)[2][2] = (z) , \
   (A)[0][1] = (A)[0][2] = (A)[1][0] = \
   (A)[1][2] = (A)[2][0] = (A)[2][1] = 0.0 )

#define LOAD_ROTX_MAT(A,th) LOAD_ROTGEN_MAT(A,th,0,1,2)
#define LOAD_ROTY_MAT(A,th) LOAD_ROTGEN_MAT(A,th,1,2,0)
#define LOAD_ROTZ_MAT(A,th) LOAD_ROTGEN_MAT(A,th,2,0,1)

#define LOAD_ROT_MAT(A,th,i)                    \
  do{ switch( (i) ){                            \
        case 0: LOAD_ROTX_MAT(A,th)   ; break ; \
        case 1: LOAD_ROTY_MAT(A,th)   ; break ; \
        case 2: LOAD_ROTZ_MAT(A,th)   ; break ; \
       default: LOAD_DIAG_MAT(A,1,1,1); break ; \
      } } while(0)

/* vector-vector multiply */

#define DOTVEC(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

/* matrix-vector multiply */

#define MATVEC(y,A,x) \
  ( (y)[0] = (A)[0][0] * (x)[0]  \
            +(A)[0][1] * (x)[1]  \
            +(A)[0][2] * (x)[2] ,\
    (y)[1] = (A)[1][0] * (x)[0]  \
            +(A)[1][1] * (x)[1]  \
            +(A)[1][2] * (x)[2] ,\
    (y)[2] = (A)[2][0] * (x)[0]  \
            +(A)[2][1] * (x)[1]  \
            +(A)[2][2] * (x)[2] )

/* vector-matrix-vector multiply */

#define VECMATVEC(val,y,A,x) \
do{ double vtmp[3]; \
    MATVEC(vtmp,A,x) ; \
    (val) = DOTVEC((y),vtmp); } while(0)

/* matrix-matrix multiply */

#define ROW_DOT_COL(A,B,i,j) (  (A)[i][0] * (B)[0][j] \
                              + (A)[i][1] * (B)[1][j] \
                              + (A)[i][2] * (B)[2][j] )

#define MAT_SCALAR(AS,A,c)                          \
  ( AS[0][0] = (c)*(A)[0][0] ,  \
    AS[1][0] = (c)*(A)[1][0] ,  \
    AS[2][0] = (c)*(A)[2][0] ,  \
    AS[0][1] = (c)*(A)[0][1] ,  \
    AS[1][1] = (c)*(A)[1][1] ,  \
    AS[2][1] = (c)*(A)[2][1] ,  \
    AS[0][2] = (c)*(A)[0][2] ,  \
    AS[1][2] = (c)*(A)[1][2] ,  \
    AS[2][2] = (c)*(A)[2][2] )

#define MAT_MUL(AB,A,B)                   \
  ( AB[0][0] = ROW_DOT_COL((A),(B),0,0) , \
    AB[1][0] = ROW_DOT_COL((A),(B),1,0) , \
    AB[2][0] = ROW_DOT_COL((A),(B),2,0) , \
    AB[0][1] = ROW_DOT_COL((A),(B),0,1) , \
    AB[1][1] = ROW_DOT_COL((A),(B),1,1) , \
    AB[2][1] = ROW_DOT_COL((A),(B),2,1) , \
    AB[0][2] = ROW_DOT_COL((A),(B),0,2) , \
    AB[1][2] = ROW_DOT_COL((A),(B),1,2) , \
    AB[2][2] = ROW_DOT_COL((A),(B),2,2)  )

#define TRANSPOSE_MAT(AT,A) \
 ( AT[0][0] = (A)[0][0] , \
   AT[1][0] = (A)[0][1] , \
   AT[2][0] = (A)[0][2] , \
   AT[0][1] = (A)[1][0] , \
   AT[1][1] = (A)[1][1] , \
   AT[2][1] = (A)[1][2] , \
   AT[0][2] = (A)[2][0] , \
   AT[1][2] = (A)[2][1] , \
   AT[2][2] = (A)[2][2] )

/*****************************************/

float   fmvmv(float *x,float **b,float *y,int m,int n);
float **fmtrn33( float **AT, float **A);
float **fmrot33( float **rotx, float **x, float **rotmat);
float **fmmul33( float **AB, float **A, float **B);
float row_dot_col33( float **A, float **B, int i, int j);
#define fmvmv33(x,b,y) fmvmv((x),(b),(y),3,3)

double dmvmv(double *x, double **b, double *y, int m,int n) ;

#endif
