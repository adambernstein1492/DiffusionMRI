/* gsl_utils.h */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>

/************************** vector routines *******************************/

gsl_vector *dv2gsl( gsl_vector *gvec, double *x, int nn ) ;
    double *gsl2dv( double *x, gsl_vector *gvec, int nn ) ;

/************************** matrix routines *******************************/

gsl_matrix *dm2gsl( gsl_matrix *gmat, double **x, int nr, int nc ) ;
   double **gsl2dm( double **x, gsl_matrix *gmat, int nr, int nc ) ;
