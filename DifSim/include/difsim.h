/* difsim.h */
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include <time.h>
#include <sys/time.h>
#include <sys/times.h>

/* shd */
//#include "shd.h"
int sphlin( float *x, float *y, float *z, int npolr, int nazim );
int sphtes( float *gx, float *gy, float *gz, int *tris, int polytope, int ntess ) ;
int sphrnd( float *x, float *y, float *z, int n ) ;

#include <X11/Intrinsic.h> // Boolean
/* vtk */
#include "vtkInteractorStyleTerrain.h"
#include "vtkDecasteljauWidget.h"
#include "vtkSphereSource.h"
#include "vtkInteractorStylePointPicker.h"
#include "vtkXGtkRenderWindowInteractor.h"
#include <vtkCommand.h>
#include <vtkCallbackCommand.h>
#include <vtkXRenderWindowInteractor.h>
#include <vtkCellArray.h>
#include <vtkLight.h>
#include <vtkCylinder.h>
#include <vtkGlyph3D.h>
#include <vtkPointSource.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkCylinderSource.h>
#include <vtkDataSetMapper.h>
#include <vtkConvexPointSet.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <vtkMapper.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkOutlineFilter.h>
#include <vtkCubeSource.h>
#include <vtkRenderer.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkCamera.h>
#include "vtkMyExtractPolyDataGeometry.h"
#include "colordefs.h"	      /* Color definitions for vtk graphics */

/* Physical constants */

/* Diffusion coefficients in mm^2/s */

/* The Diffusivity of White Matter (approx!)
   Pierpaoli + Basser, MRM 36, 893-906 (1996)
   Pierpaoli, Jezzard, Basser, Barnett, diChiro, 
   Radiology 201, 637-648 (1996).
   Gray matter diffusivity (a guess!)*/

#define Dwm   .75 * .001	/* white matter diffusion coeff */
#define Dgm   .75 * .01		/* gray matter diffusion coeff */
#define Tpm   .25 		/* tissue permeability (0<=Tpm<=1 */
#define CM2UM .0001		/* convertion factor for cm to um */
#define BUNDLE_NAME_LEN 32
#define BUNDLE_HEADER "DIFSIM2_BUNDLE"

#define SPH_LEN(x) ((x).pos[0] * (x).pos[0] +\
                    (x).pos[1] * (x).pos[1] +\
                    (x).pos[2] * (x).pos[2])
#define SPH_SUB(x,y,z) {\
(z).pos[0] = (x).pos[0] - (y).pos[0];\
(z).pos[1] = (x).pos[1] - (y).pos[1];\
(z).pos[2] = (x).pos[2] - (y).pos[2]; }

float Gamma = 26752.0 ; /* units of rad/(sec-Gauss) */

#ifndef BIGVAL
#  define BIGVAL 1.e+10
#endif

/* Logical */

#if 0
#define FALSE 0
#define TRUE  1
#endif

#ifndef YES
#define YES 1 
#endif
#ifndef NO
#define NO 0 
#endif
#ifndef MAYBE
#define MAYBE -1
#endif

/* Math */

#ifndef EVEN
#define EVEN(x) (((x) % 2) ? (0) : (1))
#endif
#ifndef ODD
#define ODD(x)  (((x) % 2) ? (1) : (0))
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define NULL_IVAL -666 
#define NULL_FVAL -666. 
#define EXIST(a)  ( (a) != NULL_FVAL )
#define BETWEEN(x,a,b)  (((x) >= (a)) && ((x)<=(b)))

#ifndef MAX
#  define MAX(a,b) (((a)<(b)) ? (b) : (a))
#endif

#ifndef MIN
#  define MIN(a,b) (((a)>(b)) ? (b) : (a))
#endif

#ifndef ABS
#  define ABS(x)  (((x)<0) ? -(x) : (x))
#endif

#define ISVALID_VOX_CLASS(mm) ( !strcmp((mm),voxel_class[0] ) || \
				!strcmp((mm),voxel_class[1] ) || \
				!strcmp((mm),voxel_class[2] ) || \
				!strcmp((mm),voxel_class[3] ) )

static char *voxel_class[] = { "NONE" , "ISO" , "SINGLE" , "MULTI" } ;
#define NUM_VOX_CLASS (sizeof(voxel_class)/sizeof(char *))

#define PRINT_FPT(fv,name) \
  printf("%s(x,y,z)         = (%13.6g,%13.6g,%13.6g)\n", name,fv[0],fv[1],fv[2]);

#define PRINT_FV(v,n,label) \
do{ int i; \
    for(i=0; i<(n); i++) { \
        printf("%s[%i] = %13.6g\n",(label),i,(float)(v)[i]); \
    }\
} while(0) 

/**************** tesselation stuff from shd.h ****************************/

#define MAX_TESS_PTS 10000

/* Sampling type */

static char *SAMPLE_TYPES[] = \
{"Polytope tessellation on sphere","Linear spherical","Random spherical"} ;

#define SAMP_TESS 0
#define SAMP_SPHL 1
#define SAMP_SRND 2

static int SAMPLE_TYPE = SAMP_TESS; /* default = polytope tessellation */

/* Tessellation polytopes */

static char *POLY_TYPES[] = {"Octahedral","Tetrahedral","Icosahedral"} ;

#define POLY_OCT 0
#define POLY_TET 1
#define POLY_ICO 2

static int POLY_TYPE = POLY_ICO; /* default = icosahedral */

/**************** end: tesselation stuff ****************************/

/* lrf - structs */

typedef struct {
  char  *vclass;
  int   ir, ic, iz ;
  short rai_mask ;
  float rai ;
  float eval1 ;
  float *evec1 ;
  float *couple ;
  float *scpl ;
} Voxel ;

typedef struct {
  int   vcnt ;
  int   nvox ;
  int   nrows, ncols ;
  int   nslices ;
  float min_rai, max_rai ;
  float rai_thresh ;
  Voxel **data ;
} Map;

/* Vector */

#include "vector.h"

/* lrf function prototypes */

void instructions(),send_help();
int  parse_args(int argc, char **argv) ;
int  initMap( int nrows , int ncols, int nslices ) ;
static int addVoxel( Map *map, 
		     int ir , int ic , int isl ,
		     float rai , float eval1 , float *evec1 ,
		     char *vclass ) ;
static void freeMap( Map * map );
static void coupling_coeff( Map *map ) ;

/* Allocation prototypes */

char  **allocate2D( int nrows, int ncols, int element_size ) ;
char ***allocate3D( int nrows, int ncols, int nslcs, int element_size ) ;
void free2D( char  **a , int nrows ) ;
void free3D( char ***a , int nslcs, int nrows ) ;

double *dvrndn(double *x,int n);
long int clock_seed() ;

/* PComplex means "pseudo-complex" 
   vtk has no complex numbers (!) */ 

#define PComplex2Complex(zv,pzv,n) \
do{ int i; \
    for(i=0; i<(n); i++) { \
        (zv)[i].r = (pzv)[2*i] ; \
        (zv)[i].i = (pzv)[2*i+1] ; \
    }\
} while(0) 

#define Complex2PComplex(pzv,zv,n) \
do{ int i; \
    for(i=0; i<(n); i++) { \
        (pzv)[2*i]   = (zv)[i].r ; \
        (pzv)[2*i+1] = (zv)[i].i; \
    }\
} while(0) 

/* Change the sign of every element of a vector */
#define NEGATE(x,n) \
do{ int i; \
    for(i=0; i<(n); i++) { \
      (x)[i] = -(x)[i] ; \
    }\
} while(0) 

