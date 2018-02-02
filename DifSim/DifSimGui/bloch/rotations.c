/* rotations.c */

#include "rotations.h"

/**<RotationMatrix3D>**************************************************

    Title: RotationMatrix3D(polr,azim,roll)

  Purpose: Create 3D Rotation Matrix
 
     Call: double **RotationMatrix3D(double **R, double polr, double azim, double roll)

    Input: polr = polar angle
           azim = azimuthal angle
           roll = roll angle

   Output: R[3][3]

   theta = polar angle (radians)
   phi   = azimuthal angle (radians)
   psi   = roll (radians)

*/
double **RotationMatrix3D(double **R, double polr, double azim, double roll)
{
  double tht = polr;
  double phi = azim;
  double psi = roll;

  double Ctht = cos(tht) ;
  double Cphi = cos(phi) ;
  double Cpsi = cos(psi) ;

  double Stht = sin(tht) ;
  double Sphi = sin(phi) ;
  double Spsi = sin(psi) ;

  R[0][0] =  Cpsi*Cphi - Ctht*Spsi*Sphi ;
  R[1][0] = -Spsi*Cphi - Ctht*Cpsi*Sphi ;
  R[2][0] =  Stht*Sphi ;

  R[0][1] =  Ctht*Spsi*Cphi + Cpsi*Sphi ;
  R[1][1] =  Ctht*Cpsi*Cphi - Spsi*Sphi ;
  R[2][1] = -Stht*Cphi ;

  R[0][2] =  Stht*Spsi ;
  R[1][2] =  Stht*Cpsi ;
  R[2][2] =  Ctht ;

  return(R);
}

/**<CylRotationMatrix3D>**************************************************

    Title: CylRotationMatrix3D(elev,azim)

  Purpose: CylRotationMatrix3D(elev,azim)
 
     Call:

    Input: 

   Output: 

  elev = elevation (radians)
  azim = azimuth (radians)

  Cylindrical symmetry
  Hsu,Mori, Eqn [20] 
  (theta = elev, phi = azim)

*/
double **CylRotationMatrix3D(double **R, double elev, double azim)
{
  double thet = elev;
  double phi  = azim;

  double Cphi = cos(phi);  
  double Sphi = sin(phi);
  double Ctht = cos(thet); 
  double Stht = sin(thet);

  R[0][0] =  Sphi ;
  R[1][0] =  Cphi * Ctht ;
  R[2][0] =  Cphi * Stht ;

  R[0][1] = -Cphi ;
  R[1][1] =  Sphi * Ctht ;
  R[2][1] =  Sphi * Stht ;

  R[0][2] =  0. ;
  R[1][2] = -Stht;
  R[2][2] =  Ctht;

  return(R);
}
