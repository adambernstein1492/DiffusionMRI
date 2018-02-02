/*<sphfun.h>***************************************************************/

#ifndef _LRF_SPH_HEADER_
#define _LRF_SPH_HEADER_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

/* lrf includes */
#include "alm.h"
#include "gsl_utils.h"
#include "cs.h"
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

/**************** SHD routines ******************************/

typedef struct SHD {           /* Gradient structure */

  int Lmin, Lmax ;		/* min,max L for recon */
  int      nvert ;		/* number of grad directions */
  complex ***Ylm ;		/* spherical harmonics */
  float     *wt  ;		/* voronoi weights */

  float    *rad  ;		/* radius (really just used in conversion) */
  float    *pol  ;		/* polar angles */
  float    *azi  ;		/* azi angles */

} SHD ;

/* Spherical Harmonics functions */

float **ylm_scale( float **Ylm_scale, int Lmin, int Lmax) ;

complex *SphericalHarmonicY(complex *Ylm, int l, int m, 
			    float *polr, float *azim, int npts) ;

complex ***SphericalHarmonicsYlm( complex ***Ylm, int Lmin, int Lmax,
				  float *polr, float *azim, int npts ) ;

complex **SpharmTransformV1( complex **Alm, float *data , 
			     int Lmin, int Lmax, 
			     float *pol, float *azi, int nvert) ;

complex **SpharmTransform( complex **Alm, float *data, float **gvecs, SHD *shd) ;
complex **SpharmTransformOpt( complex **Alm, float *data, SHD *shd) ;

#define YLM_COEF(l,m) \
(((2*(l)+1)/(4.*PI))*((float)factorial((l)-(m))/(float)factorial((l)+(m))))

/* END SHD Stuff */

/**************** SWD routines ******************************/

float ***jl_scale( float ***jlscal, int Lmax, float *jrad, int nrad ) ;

complex ***SphWaveTransform( complex ***Almr, float **data , 
			     int Lmin, int Lmax, int nvert,
			     float *rad, int nrad,
			     complex ****Wlmr, float *wt );

void o_dot(complex ***plmr, float *data,
	   int Lmin, int Lmax, int nvert,
	   float *rad, int nrad,
	   complex ****alpha_llmr, complex ***Ylm, float *wt, float t);


complex ****SphericalWavesWlmr( complex ****Wlmr, 
				float *jrad, int nrad,
				int   Lmin, int   Lmax,
				float *polr, float *azim, int nvert );

float **SphericalBesseljl( float **jl, int Lmin, int Lmax,
			   float *rad, int nrad );

/* END SWD Stuff */

/* ************************************ */

int sphlin( float *x, float *y, float *z, int npolr, int nazim );
int sphtes( Grad *grad ) ;
int sphrnd( float *x, float *y, float *z, int n ) ;
float *sbesselj(float *sbj, int l, float *x, int n ) ;

void *print_Omega( Omega *omega ) ;
Omega *alloc_omega( int nazim, int npolr ) ;
Omega *gen_omega( int *nazim, int *npolr, int setn, int rev ) ;
void get_angles_from_omega( Omega *omega , 
			    float *azim, float *polr,
			    int nazim, int npolr) ;

#endif
