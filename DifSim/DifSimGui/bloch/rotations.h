/* rotations.h */

/* Generic includes */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/fcntl.h> 

double **RotationMatrix3D(double **R, double polr, double azim, double roll) ;
double **CylRotationMatrix3D(double **R, double elev, double azim) ;
