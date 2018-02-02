/* matrix.c */
#include "matrix.h"

/**<dmadd>**************************************************

    Title: Double Matrix ADDition

  Purpose: Compute sum of two double matrices.
 
     Call: double **dmadd(z,x,a,y,nr,nc)
                double **z,**x,a,**y;
		int    nr,nc;

    Input:  z[nr][nc] = double matrix
            x[nr][nc] = double matrix
	            a = double scalar
            y[nc][nc] = double matrix
	           nr = number of rows
	           nc = number of cols

   Output: double **z[nr][nc] = x[nr][nc] + a*y[nr][nc]

    Notes: Matrices must be of same size
           Handles cases of a=1,-1 specially
	   Returns pointer to z[] array
*/
double **dmadd(double **z,double **x,double a,double **y,int nr,int nc)
{
        register int    r,c;
	double   **zold = z;

	if(nr>0 && nc>0)
	{
		if(a==1.)
		{
		        for(r=0; r<nr; r++) {
			  for(c=0; c<nc; c++) {
			    z[r][c] = x[r][c] + y[r][c];
			  }
			}
		}
		else if(a== -1.)
		{
		        for(r=0; r<nr; r++) {
			  for(c=0; c<nc; c++) {
			    z[r][c] = x[r][c] - y[r][c];
			  }
			}
		}
		else
		{
		        for(r=0; r<nr; r++) {
			  for(c=0; c<nc; c++) {
			    z[r][c] = x[r][c] + a * y[r][c];
			  }
			}
		}
	}
	return(zold) ;
}

/**<dvmmul>*******************************************************

    Title: Double Vector-Matrix MULtiplication

  Purpose: Compute product of double vector and double matrix.
 
     Call: double *dvmmul(d,f,a,m,n)
                double *d,**f,*a;
	        int    m,n;

   Returns: double d[n] = f[n][m] a[m] 

*/
double *dvmmul(double *d,double **f, double *a,int m,int n)
{
        register int    i,j;
	register double sum;

	for(j=0; j< n; j++) {
	  sum = 0.0;
	  for(i=0; i< m; i++) {
	    sum += f[j][i] * a[i];
	  }
	  d[j] = sum;
	}
	return(d);
}

float row_dot_col33( float **A, float **B, int i, int j)
{
  float suma ;

  suma = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j] ;

  return(suma);

}
/* matrix-matrix multiply */

float **fmmul33( float **AB, float **A, float **B)
{

  AB[0][0] = row_dot_col33(A,B,0,0) ;
  AB[1][0] = row_dot_col33(A,B,1,0) ;
  AB[2][0] = row_dot_col33(A,B,2,0) ;
  AB[0][1] = row_dot_col33(A,B,0,1) ;
  AB[1][1] = row_dot_col33(A,B,1,1) ;
  AB[2][1] = row_dot_col33(A,B,2,1) ;
  AB[0][2] = row_dot_col33(A,B,0,2) ;
  AB[1][2] = row_dot_col33(A,B,1,2) ;
  AB[2][2] = row_dot_col33(A,B,2,2) ;

  return(AB);
}

/* matrix transpose */

float **fmtrn33( float **AT, float **A)
{

  AT[0][0] = A[0][0] ;
  AT[1][0] = A[0][1] ;
  AT[2][0] = A[0][2] ;
  AT[0][1] = A[1][0] ;
  AT[1][1] = A[1][1] ;
  AT[2][1] = A[1][2] ;
  AT[0][2] = A[2][0] ;
  AT[1][2] = A[2][1] ;
  AT[2][2] = A[2][2] ;

  return(AT);
}

/* Rotate x by rotmat to get rotx */

float **fmrot33( float **rotx, float **x, float **rotmat)
{
  float **tmpmat1 = fm_alloc(3,3);
  float **tmpmat2 = fm_alloc(3,3);

  fmtrn33(tmpmat1,rotmat);  
  fmmul33(tmpmat2,tmpmat1,x);
  fmmul33(rotx,tmpmat2,rotmat);

  fm_free(tmpmat1,3,3);
  fm_free(tmpmat2,3,3);

  return(rotx);
}

/**<fmmul>**************************************************

    Title: Float Matrix MULtiplication

  Purpose: Compute product of two float matrices.
 
     Call: float **fmmul( float **z, float **x, float **y,
		int nrx, int ncx, int ncy)

    Input:  z[nrx][ncy] = float matrix
            x[nrx][ncx] = float matrix
            y[ncx][ncy] = float matrix
	            nrx = number of rows in x
	            ncx = number of cols in x and rows in y
	            ncy = number of cols in y

   Output: float **z[nrx][ncy] = x[nrx][ncx] * y[ncx][ncy]

    Notes: Inner dimensions of x and y must be equal

*/
float **fmmul( float **z, float **x, float **y,
		int nrx, int ncx, int ncy)
{
        register int    i,j,k;
	register float sum;
	float   **zold = z;

	for(k=0; k< nrx; k++) {
	  for(j=0; j< ncy; j++) {
	    sum = 0.0;
	    for(i=0; i< ncx; i++) {
	      sum += x[k][i] * y[i][j];
	    }
	    z[k][j] = sum;
	  }
	}

	return(zold);
}
/**<fmtrn.c>************************************************************

     Title: Float Matrix TRANSpose

   Purpose: Compute transpose of float matrix

      Call: float **fmtrn( float **z, float**x, int nr, int nc)

     Input: z[nc][nr] = float matrix
            x[nr][nc] = input matrix
		   nr = number of rows in x
		   nc = number of cols in x

    Output: z[nc][nr] = transpose of x
            x[nr][nc] = (unchanged)
                   nr = (unchanged)
                   nc = (unchanged)

*/
float **fmtrn( float **z, float**x, int nr, int nc)
{
  int   ir,ic;
  float **zold = z;

  for (ir=0; ir<nr; ir ++) {
    for (ic=0; ic<nc; ic ++) {
      z[ic][ir] = x[ir][ic];
    }
  }
  return(zold);
}

/**<fmvmv.c>********************************************************

      Title: Float Vector Matrix Vector multiply

    Purpose: X'*B*Y

       Call: float fmvmv(x,b,y,m,n)
                  float *x[m],**b[m][n],*y[n];
		  int    m,n;

*/
float fmvmv(float *x,float **b,float *y,int m,int n)
{
  float suma = 0.0;
  int    i,j;

  for(i=0; i<m; i++) {
    for(j=0; j<n; j++) {
      suma += x[i] * b[i][j] * y[j];
    }
  }

  return(suma);
}

/**<dmvmv.c>********************************************************

      Title: Double Matrix Vector Matrix Vector multiply

    Purpose: X'*B*Y

       Call: double dmvmv(x,b,y,m,n)
                  double *x[m],**b[m][n],*y[n];
		  int    m,n;

*/
double dmvmv(double *x, double **b, double *y, int m,int n)
{
  double suma = 0.0;
  int    i,j;

  for(i=0; i<m; i++) {
    for(j=0; j<n; j++) {
      suma += x[i] * b[i][j] * y[j];
    }
  }

  return(suma);
}
