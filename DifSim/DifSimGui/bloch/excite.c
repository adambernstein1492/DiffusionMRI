#include "bloch.h"

/**<rotation_matrix_rot_3d.c>***********************************************
 
    Title: rotation_matrix_rot_3d
 
  Purpose: Matrix for 3D rotations in terms of rotations angles

     Call: double **rotation_matrix_rot_3d(double **R, double theta, double phi, double xi)
 
    Input: R[3][3] = rotation matrix
             theta = first rotation angle
               phi = second rotation angle
                xi = third rotation angle

   Returns: double **R[3][3]
*/
double **rotation_matrix_rot_3d(double **R, double theta, double phi, double xi)
{
  double nx,ny,nz;
  double nx2,ny2,nz2;
  double sx,hsx,hsx2;

  /* Compute rotation axes */

  nx = sin(theta) * cos(phi);
  ny = sin(theta) * sin(phi);
  nz = cos(theta);

  nx2  = nx * nx;  
  ny2  = ny * ny;
  nz2  = nz * nz;

  /* Compute trig relations */

  sx   = sin(xi);
  hsx  = sin(.5*xi);
  hsx2 = 2. * hsx * hsx;

  /* Compute rotation matrix components */

  R[0][0] =  1. -  (ny2 + nz2) * hsx2;
  R[0][1] =  nz * sx + nx * ny * hsx2;
  R[0][2] = -ny * sx + nz * nx * hsx2;

  R[1][0] = -nz * sx + nx * ny * hsx2;
  R[1][1] =  1. -  (nz2 + nx2) * hsx2;
  R[1][2] =  nx * sx + ny * nz * hsx2;

  R[2][0] =  ny * sx + nz * nx * hsx2;
  R[2][1] = -nx * sx + ny * nz * hsx2;
  R[2][2] =  1. -  (nx2 + ny2) * hsx2;

  /* Return rotation matrix */

  return(R);
}
/**<rotation_matrix_fm.c>***********************************************
 
    Title: rotation_matrix_fm
 
  Purpose: Matrix for rotations to fm frame

     Call: double **rotation_matrix_fm(R,phi)
                double **R,phi;
 
    Input: R[3][3] = rotation matrix
               phi = signal phase

   Returns: double **R[3][3]
*/
double **rotation_matrix_fm( double **R, double phi )
{
  double cp = cos(phi);
  double sp = sin(phi);

  /* Compute rotation matrix components */

  R[0][0] =  cp;
  R[0][1] =  sp;
  R[0][2] = 0.0;

  R[1][0] = -sp;
  R[1][1] =  cp;
  R[1][2] = 0.0;

  R[2][0] = 0.0;
  R[2][1] = 0.0;
  R[2][2] = 1.0;

  /* Return rotation matrix */

  return(R);
}
/**<relaxation_matrix_rot_3d.c>***********************************************
 
    Title: relaxation_matrix_rot_3d
 
  Purpose: Matrix for 3D relaxation in terms of rotations

     Call: double **relaxation_matrix_rot_3d( double **E, 
                                              double theta, double phi, double xi, 
					      double we, double t1, double t2)
 
    Input: E[3][3] = rotation matrix
             theta = first rotation angle
               phi = second rotation angle
                xi = third rotation angle

   Returns: double **E[3][3]
*/
double **relaxation_matrix_rot_3d( double **E, 
				   double theta, double phi, double xi, 
				   double we, double t1, double t2)
{
  double t1sq,t2sq,T;
  double nx,ny,nz;
  double cx,sx,hsx,hsx2;

  /* Compute rotation axes */

  nx = sin(theta) * cos(phi);
  ny = sin(theta) * sin(phi);
  nz = cos(theta);

  /* Compute relaxation relations */

  t1  *= we;
  t2  *= we;
  t1sq = t1 * t1;
  t2sq = t2 * t2;
  T    = (t1 * t2)/(t1 + t2);

  /* Compute trig relations */

  cx   = cos(xi);
  sx   = sin(xi);
  hsx  = sin(.5*xi);
  hsx2 = 2. * hsx * hsx;

  /* Compute matrix components */

  E[0][0] =  hsx2 / t2sq - sx / t2;
  E[0][1] = -2. * nz * hsx2 / t2;
  E[0][2] =  ny * hsx2 / T;

  E[1][0] = -E[0][1];
  E[1][1] =  hsx2 / t2sq - sx / t2;
  E[1][2] = -nx * hsx2 / T;

  E[2][0] = -E[0][2];
  E[2][1] = -E[1][2];
  E[2][2] = hsx2 / t1sq - sx / t1;

  /* Return relaxation matrix */
  return(E);
}

/**<propagate_rot.c>********************************************************
 
    Title: propagate_rot
 
  Purpose: Evolution matrix for spin magnetization - rotation operator method

     Call: double *propagate_rot(double *M, double **R, double **E, 
                                 double w1, double p1, double dw, 
				 double t1, double t2, double dt)
 
 
    Input:    M[3] = magnetization vector
           R[3][3] = rotation matrix
           E[3][3] = relaxation matrix
                w1 = RF angular frequency
                p1 = RF phase
	        dw = Offset angular frequency
	        t1 = spin-lattice relaxation time
	        t2 = spin-spin relaxation time
  	        dt = time increment

  Routines called: rotation_matrix_rot_3d()
                   relaxation_matrix_rot_3d()
		   minfinity_rot()

   Returns: double *M[3]
*/
double *propagate_rot(double *M, double **R, double **E, 
		      double w1, double p1, double dw, 
		      double t1, double t2, double dt)
{
  char *program_name= "propagate_rot";
  double Mnew[3],Minf[3];
  double we,theta,phi,xi,sign;
  bool debug_local = 0;

  if (debug_local) 
    fprintf(stderr,"********** %s ************\n",program_name);

  /* Effective angular frequency */

  if(dw==0) {
    sign = -1.;
  } else {
    sign = dw/fabs(dw);
  }
  we = -sign*sqrt(w1*w1 + dw*dw); 
  if (debug_local) fprintf(stderr,"(w1,dw,we) = (%lf,%lf,%lf)\n",w1,dw,we);

  /* Rotation angles */

  theta = we==0.?0:asin(w1/we);  /* or -> theta = acos(dw/we); */
  phi   = p1;
  xi    = we * dt;

  if (debug_local) fprintf(stderr,"(theta,phi,xi) = (%lf,%lf,%lf)\n",theta,phi,xi);
  
  /* Ideal rotation matrix */

  rotation_matrix_rot_3d(R,theta,phi,xi);
  if (debug_local) PRINT_MAT33("R",R);

  /* Relaxation matrix */

  relaxation_matrix_rot_3d(E,theta,phi,xi,we,t1,t2);
  if (debug_local) PRINT_MAT33("E",E);

  /* Full propagation matrix R + E */

  dmadd(R,R,1.,E,3,3);
  if (debug_local) PRINT_MAT33("R",R);

  /* Steady-state magnetization */

  minfinity_rot(Minf,theta,phi,xi,we,t1,t2);
  dvfill(Minf,0.,3);		/* TEMP AHB!!! */
  if (debug_local) PRINT_VEC3("\tMinf =",Minf);

  /*if (dvnan(Minf,3)) {
    perror("Minf is nan!");
    exit(-1);
    }*/

  if (1) {	     /* AHB until we find source of nan's in minfinity_rot! */
    int i;
    bool uhoh = 0;
    for (i=0; i<3; i++) {
      if isnan(Minf[i]) {
	  Minf[i] = 0;
	  uhoh = 1;
	}
    }
    if (uhoh) {
      fprintf(stderr,"propagate_rot: Minf is nan!  Setting to:\n");
      PRINT_VEC3("\tMinf =",Minf);
    }
  }

  /* Propagate the magnetization vector */

  dvadd(M,M,-1.,Minf,3);
  if (debug_local) PRINT_VEC3("M",M);
  dvmmul(Mnew,R,M,3,3);
  if (debug_local) PRINT_VEC3("Mnew",Mnew);
  dvadd(M,Mnew,1.,Minf,3); 
  if (debug_local) PRINT_VEC3("M",M);

  return(M);
}
/**<minfinity_rot.c>***********************************************
 
    Title: M INFINITY - ROTation method
 
  Purpose: Calculate equilibrium magnetization in rotation method

     Call: double *minfinity_rot(double *M, double theta, double phi, 
                                 double xi, double we, double t1, double t2)
 
    Input: M[3] = rotation matrix
          theta = first rotation angle
            phi = second rotation angle
             we = effective angular frequency
             t1 = longitudinal relaxation time
             t2 = transverse relaxation time

   Returns: double *M[3]

      Note: At zero frequency offset, for t2 = Inf, the actual
            value of t2 (as defined in bloch.c, for example)
	    will cause severe rounding errors in the last
	    step of this routine.  That is why Inf = BIG
	    rather than HUGE (=MAXFLOAT) in bloch.c.
*/
#define M0 1.0

double *minfinity_rot(double *M, 
		      double theta, double phi, double xi,
		      double we, double t1, double t2)
{
  double nx,ny,nz;
  double nx2,ny2,nz2;
  double mx,my,mz;
  double t2sq;
  double d;

  bool debug_local= 0;

  /* Compute rotation axes */

  nx = sin(theta) * cos(phi);
  ny = sin(theta) * sin(phi);
  nz = cos(theta);

  if (debug_local) fprintf(stderr,"(nx,ny,nz) = (%lf,%lf,%lf)\n",nx,ny,nz);

  nx2 = nx * nx;
  ny2 = ny * ny;
  nz2 = nz * nz;

  t1   *= we;
  t2   *= we;
  t2sq  = t2 * t2;

  /* Compute matrix components */

  d = 1. + (nx2 + ny2) * t1 * t2 + nz2 * t2sq;

  mx = M0 * t2 * (-ny + nx * nz * t2);
  my = M0 * t2 * (nx + ny * nz * t2);
  mz = M0 * (1. + nz2 * t2sq); 

  if (debug_local) fprintf(stderr,"(mx,my,mz) = (%lf,%lf,%lf)\n",mx,my,mz);

  M[0] = mx / d;
  M[1] = my / d;
  M[2] = mz / d;

  /* Return magnetization */

  return(M);
}
