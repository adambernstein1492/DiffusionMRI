/*<grads.h>***************************************************************/

#ifndef _LRF_GRADS_HEADER_
#define _LRF_GRADS_HEADER_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/* Tessellation polytopes */

static char *POLY_TYPES[] = {"Octahedral","Tetrahedral","Icosahedral"} ;

#define POLY_OCT 0
#define POLY_TET 1
#define POLY_ICO 2

static int POLY_TYPE = POLY_ICO; /* default = icosahedral */

/*--- Gradient definitions -------------------------------------------------*/

typedef struct BMATRIX {	/* bmatrix */
  float  *bvals;		/* bvalues */
  float **bvecs;		/* bvectors */
  float **bmatrix;		/* copath bmatrix (bvals * averaged bvecs) */
} BMATRIX ;

typedef struct DUP {		/* Duplicate structure */
  int  n;			/* number of duplicate */
  int *v;			/* vector of duplicate indeces */
} DUP ;

typedef struct Grad {		/* Gradient structure */

  int   polytope ;		/* polytope to tesselate */
  int   ntess ;			/* tesselation level */

  float *gx ;			/* gx[nvert] - x grad */
  float *gy ;			/* gy[nvert] - y grad */
  float *gz ;			/* gz[nvert] - z grad */
  double **gxyz ;		/* gxyz[nvert][3] = (gx[nvert],gy[nvert],gz[nvert]) */
  int   nvert ;			/* number of grad directions */

  int   nb0 ;			/* number of true b=0 images (BOTH ste and hyp b=0)*/
  int   nb0_ste ;		/* number of ste b=0 images (ste has b=0 but hyp b!=0) */
  int   nb0_hyp ;		/* number of hyp b=0 images (hyp has b=0 but ste b!=0) */

  int   nbv_mix ;		/* number of diffusion weighted images (both ste and hyp on) */
  int   nbv_ste ;		/* number of ste diffusion weighted images (no hyp diff weighting) */
  int   nbv_hyp ;		/* number of hyp diffusion weighted images (no ste diff weighting) */

  int   nb_total ;		/* total number images: nb0 + nbv_ste + nbv_hyp + nbv_mix */

  int   ndup_mix ;		/* number of duplicate diffusion weighted images (both ste and hyp on) */
  int   ndup_ste ;		/* number of duplicate ste diffusion weighted images (no hyp diff weighting) */
  int   ndup_hyp ;		/* number of duplicate hyp diffusion weighted images (no ste diff weighting) */

  int   *indx_b0 ;		/* locations of nb0 images */
  int   *indx_b0_ste ;		/* locations of nb0_ste images */
  int   *indx_b0_hyp ;		/* locations of nb0_hyp images */

  int   *indx_bv_mix ;		/* locations of nbv_ste and nbv_hyp images */
  int   *indx_bv_ste ;		/* locations of nbv_ste images */
  int   *indx_bv_hyp ;		/* locations of nbv_hyp images */

  DUP   **dup_mix ;		/* number of duplicates of indx_bv */
  DUP   **dup_ste ;		/* number of duplicates of indx_bv_ste */
  DUP   **dup_hyp ;		/* number of duplicates of indx_bv_hyp */

  int   *tris ;			/* tris[ntris] - triangles for graphics */
  int   ntris ;			/* number of triangles */

} Grad;

/* pp=0 don't print either grads or tris
   pp=1 print grads
   pp=2 print tris
   pp=3 print grads and tris */

#define PRINT_GRAD(grd,pmode) \
do{ int i;                 \
 if((grd) != NULL) {        \
   fprintf(stderr,           \
	   "polytope = %d\n" \
	   "ntess = %d\n"    \
	   "nvert = %d\n"    \
	   "ntris = %d\n",   \
	   (grd)->polytope,(grd)->ntess,(grd)->nvert,(grd)->ntris) ; \
   if((pmode)==1 || (pmode)==3) {       \
     for(i=0; i<(grd)->nvert; i++) {  \
       fprintf(stderr,"grad[%i] = (%f,%f,%f)\n", \
	       i,(grd)->gx[i],(grd)->gy[i],(grd)->gz[i]); \
     } \
   } \
   if((pmode)==2 || (pmode)==3) { \
     for(i=0; i<(grd)->ntris; i++) { \
       fprintf(stderr, "tris[%i] = %i\n",i,(grd)->tris[i]); \
     } \
   } \
 } else { \
   fprintf(stderr,"PRINT_GRAD:  Grad is NULL!\n") ;\
 }\
} while(0);

int sphtes( Grad *grad ) ;

#endif
