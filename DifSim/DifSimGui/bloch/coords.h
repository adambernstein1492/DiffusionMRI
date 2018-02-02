/**<coords.h>******************************************************/

#ifndef _COORD_HEADER_
#define _COORD_HEADER_

#include <stdio.h>
#include "math.h"

/************************* Useful stuff *****************************/

#ifndef PI
#  define PI 3.141592653589793238462643
#endif

#ifndef TWOPI
#  define TWOPI 6.283185307179586476925286     /* 2 * PI */
#endif

#ifndef PI2
#  define PI2 1.57079632679490     /* PI/2 */
#endif

#ifndef MAX
#   define MAX(a,b) (((a)<(b)) ? (b) : (a))
#endif

#ifndef MIN
#   define MIN(a,b) (((a)>(b)) ? (b) : (a))
#endif

#ifndef ABS
#  define ABS(x)  (((x)<0) ? -(x) : (x))
#endif

#ifndef MABS
#define MABS(x)   (((x)<0) ? -(x) : (x))
#endif


#ifndef DEG2RAD
#define DEG2RAD (M_PI/180.)
#endif

#ifndef RAD2DEG
#define RAD2DEG (180./M_PI)
#endif

/************************* Structures *****************************/

typedef struct {	/* Spherical point */
  float rad, pol, azi ; 
} spoint;

typedef struct {	/* Cartesian point */
  float  x, y, z;
} cpoint;

typedef struct {	/* Index 3 point */
    int  ix, iy, iz;
} indx3;

typedef struct {	/* Index 2 point */
    int  ix, iy ;
} indx2;

typedef struct {
  cpoint    pt[3];	/* Vertices of triangle */
  float    area;	/* Unused; might be used for adaptive subdivision */
} triangle;

typedef struct {
    int       nface;	/* # of triangles in object */
    triangle  *face;	/* Triangles */
    int       nvert;	/* # of vertices in object */
    cpoint    *vert;	/* Vertices */
    char      *type;	/* Type of polygon */
} object;

typedef struct {		/* Apr 6, 2004 */
    int       nazim;	/* # of azimuthal */
    int       npolr;	/* # of polr */
    int       nomega;	/* # of angles */
    spoint   *spts;	/* spherical points */
} Omega;

/* New allocations 24-Aug-03 */

#define sph_alloc(n) ((spoint *) calloc( (n) ,  sizeof(spoint) ))

/* printindx3(indx3 ipt,char *name) */
#define printindx3(ipt,name) \
  printf("indx3 %s = {%i,%i,%i}\n",name,ipt.ix,ipt.iy,ipt.iz);

/* printindx2(indx2 ipt,char *name) */
#define printindx2(ipt,name) \
  printf("indx2 %s = {%i,%i}\n",name,ipt.ix,ipt.iy);

/* printcpt(cpoint cpt,char *name) */
#define printcpt(cpt,name) \
  printf("cpt %s(x,y,z)         = (%13.6g,%13.6g,%13.6g)\n",name,cpt.x,cpt.y,cpt.z);

/* printspt(spoint spt,char *name) */
#define printspt(spt,name) \
  printf("spt %s(rad,pol,azi) = (%13.6g,%13.6g,%13.6g)\n",name,spt.rad,spt.pol,spt.azi);

/* Functions */

float distanceto(cpoint cpt1, cpoint cpt2) ;
cpoint newcpt( float x, float y, float z ) ;
spoint newspt( float rad, float pol, float azi ) ;
cpoint spt2cpt( spoint spt ) ;
spoint cpt2spt( cpoint cpt ) ;
void cart2polr(	float *rad, float *pol, float *azi, 
		float *cx,  float *cy,  float *cz, 
		int nn ) ;

Omega *gen_omega( int *nazim, int *npolr, int setn, int rev );
void *print_Omega( Omega *omega ) ;

#endif
