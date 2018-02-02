/**<gsl_utils.c>********************************************************

      Title: GSL UTILitieS

    Purpose: Manipulation/conversion of gsl routines/structures

*/
#include "gsl_utils.h"

/************************** vector routines *******************************/

gsl_vector *dv2gsl( gsl_vector *gvec, double *x, int nn ) 
{
  int ii ;

  for (ii=0; ii<nn; ii++) gsl_vector_set(gvec,ii,x[ii]) ;

  return (gvec) ;
}

double *gsl2dv( double *x, gsl_vector *gvec, int nn ) 
{
  int ii ;

  for (ii=0; ii<nn; ii++) x[ii] = gsl_vector_get(gvec,ii) ;

  return (x) ;
}

/************************** matrix routines *******************************/

gsl_matrix *dm2gsl( gsl_matrix *gmat, double **x, int nr, int nc ) 
{
  int ir,ic ;

  for (ir=0; ir<nr; ir++) {
    for (ic=0; ic<nc; ic++) {
      gsl_matrix_set(gmat,ir,ic,x[ir][ic]) ;
    }
  }

  return (gmat) ;

}

double **gsl2dm( double **x, gsl_matrix *gmat, int nr, int nc ) 
{
  int ir,ic ;

  for (ir=0; ir<nr; ir++) {
    for (ic=0; ic<nc; ic++) {
      x[ir][ic] = gsl_matrix_get(gmat,ir,ic) ;
    }
  }

  return (x) ;
}
