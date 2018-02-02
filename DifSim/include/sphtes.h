/*<sphtes.h>***********************************************/

#ifndef _SPHTES_HEADER_
#define _SPHTES_HEADER_

/* comment out these three includes for epic */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "coords.h" 

/*------------------------------- Octahedron -------------------------------*/

/* Vertices of a unit octahedron */
/* Six equidistant points lying on the unit sphere */

#define XPLUS {  1,  0,  0 }	/*  X */
#define XMIN  { -1,  0,  0 }	/* -X */
#define YPLUS {  0,  1,  0 }	/*  Y */
#define YMIN  {  0, -1,  0 }	/* -Y */
#define ZPLUS {  0,  0,  1 }	/*  Z */
#define ZMIN  {  0,  0, -1 }	/* -Z */

cpoint octahedron_vertices[] = { XPLUS ,
				 XMIN  ,
				 YPLUS ,
				 YMIN  ,
				 ZPLUS ,
				 ZMIN  };

/* Faces of a unit octahedron */

triangle octahedron_faces[] = {
    { { XPLUS, ZPLUS, YPLUS }, 0.0 },
    { { YPLUS, ZPLUS, XMIN  }, 0.0 },
    { { XMIN , ZPLUS, YMIN  }, 0.0 },
    { { YMIN , ZPLUS, XPLUS }, 0.0 },
    { { XPLUS, YPLUS, ZMIN  }, 0.0 },
    { { YPLUS, XMIN , ZMIN  }, 0.0 },
    { { XMIN , YMIN , ZMIN  }, 0.0 },
    { { YMIN , XPLUS, ZMIN  }, 0.0 }
};

int octahedron_ids[] = {
   0, 4, 2,
   2, 4, 1,
   1, 4, 3,
   3, 4, 0,
   0, 2, 5,
   4, 1, 5,
   1, 3, 5,
   3, 0, 5
};

/* A unit octahedron */

shape_object oct = {
    sizeof(octahedron_faces) / sizeof(octahedron_faces[0]),
    &octahedron_faces[0],
    sizeof(octahedron_vertices) / sizeof(octahedron_vertices[0]),
    &octahedron_vertices[0],
    "Octahedron"
};

/*------------------------------- Octadome -------------------------------*/

/* Vertices of a unit octadome */
/* Six equidistant points lying on the unit sphere */

cpoint octadome_vertices[] = { XPLUS ,
			       XMIN  ,
			       YPLUS ,
			       YMIN  ,
			       ZPLUS };

/* Faces of a unit octadome */

triangle octadome_faces[] = {
    { { XPLUS, ZPLUS, YPLUS }, 0.0 },
    { { YPLUS, ZPLUS, XMIN  }, 0.0 },
    { { XMIN , ZPLUS, YMIN  }, 0.0 },
    { { YMIN , ZPLUS, XPLUS }, 0.0 }
};

/* A unit octadome */

shape_object octd = {
    sizeof(octadome_faces) / sizeof(octadome_faces[0]),
    &octadome_faces[0],
    sizeof(octadome_vertices) / sizeof(octadome_vertices[0]),
    &octadome_vertices[0],
    "Octadome"
};

/*------------------------------- Tetrahedron -------------------------------*/

/* Vertices of a tetrahedron */

#define sqrt_3 0.5773502692
#define PPP {  sqrt_3,	sqrt_3,  sqrt_3 }   /* +X, +Y, +Z */
#define MMP { -sqrt_3, -sqrt_3,  sqrt_3 }   /* -X, -Y, +Z */
#define MPM { -sqrt_3,	sqrt_3, -sqrt_3 }   /* -X, +Y, -Z */
#define PMM {  sqrt_3, -sqrt_3, -sqrt_3 }   /* +X, -Y, -Z */

cpoint tetrahedron_vertices[] = { PPP ,
				  MMP ,
				  MPM ,
				  PMM };

/* Structure describing a tetrahedron */

triangle tetrahedron_faces[] = {
    { { PPP, MMP, MPM }, 0.0} ,
    { { PPP, PMM, MMP }, 0.0} ,
    { { MPM, MMP, PMM }, 0.0} ,
    { { PMM, PPP, MPM }, 0.0}
};

/* A unit tetrahedron */

shape_object tet = {
    sizeof(tetrahedron_faces) / sizeof(tetrahedron_faces[0]),
    &tetrahedron_faces[0],
    sizeof(tetrahedron_vertices) / sizeof(tetrahedron_vertices[0]),
    &tetrahedron_vertices[0],
    "Tetrahedron"
};

/*------------------------------- Tetradome -------------------------------*/

cpoint tetradome_vertices[] = { PPP ,
				MMP };

/* Structure describing a tetradome */

triangle tetradome_faces[] = {
    { { PPP, MMP, MPM }, 0.0} ,
    { { PPP, PMM, MMP }, 0.0} 
};

/* A unit tetradome */

shape_object tetd = {
    sizeof(tetradome_faces) / sizeof(tetradome_faces[0]),
    &tetradome_faces[0],
    sizeof(tetradome_vertices) / sizeof(tetradome_vertices[0]),
    &tetradome_vertices[0],
    "Tetradome"
};

/*------------------------------- Icosahedron -------------------------------*/

/* Twelve vertices of icosahedron on unit sphere */

#define tau 0.8506508084      /* t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)  */
#define one 0.5257311121      /* one=1/sqrt(1+t^2) , unit sphere     */
#define ZA {  tau,  one,    0 }
#define ZB { -tau,  one,    0 }
#define ZC { -tau, -one,    0 }
#define ZD {  tau, -one,    0 }
#define YA {  one,   0 ,  tau }
#define YB {  one,   0 , -tau }
#define YC { -one,   0 , -tau }
#define YD { -one,   0 ,  tau }
#define XA {   0 ,  tau,  one }
#define XB {   0 , -tau,  one }
#define XC {   0 , -tau, -one }
#define XD {   0 ,  tau, -one }

cpoint icosahedron_vertices[] = { ZA,ZB,ZC,ZD ,
				  YA,YB,YC,YD ,
				  XA,XB,XC,XD };

/* Structure for unit icosahedron */

triangle icosahedron_faces[] = {
    { { YA, XA, YD }, 0.0 },
    { { YA, YD, XB }, 0.0 },
    { { YB, YC, XD }, 0.0 },
    { { YB, XC, YC }, 0.0 },
    { { ZA, YA, ZD }, 0.0 },
    { { ZA, ZD, YB }, 0.0 },
    { { ZC, YD, ZB }, 0.0 },
    { { ZC, ZB, YC }, 0.0 },
    { { XA, ZA, XD }, 0.0 },
    { { XA, XD, ZB }, 0.0 },
    { { XB, XC, ZD }, 0.0 },
    { { XB, ZC, XC }, 0.0 },
    { { XA, YA, ZA }, 0.0 },
    { { XD, ZA, YB }, 0.0 },
    { { YA, XB, ZD }, 0.0 },
    { { YB, ZD, XC }, 0.0 },
    { { YD, XA, ZB }, 0.0 },
    { { YC, ZB, XD }, 0.0 },
    { { YD, ZC, XB }, 0.0 },
    { { YC, XC, ZC }, 0.0 }
};
int icosahedron_ids[] = {
    4, 8, 7,
    4, 7, 9,
    5, 6, 11,
    5, 10, 6,
    0, 4, 3,
    0, 3, 5,
    2, 7, 1,
    2, 1, 6,
    8, 0, 11,
    8, 11, 1,
    9, 10, 3,
    9, 2, 10,
    8, 4, 0,
    11, 0, 5,
    4, 9, 3,
    5, 3, 10,
    7, 8, 1,
    6, 1, 11,
    7, 2, 9,
    6, 10, 2
};


/* A unit icosahedron */

shape_object ico = {
    sizeof(icosahedron_faces) / sizeof(icosahedron_faces[0]),
    &icosahedron_faces[0],
    sizeof(icosahedron_vertices) / sizeof(icosahedron_vertices[0]),
    &icosahedron_vertices[0],
    "Icosahedron"
};

/*------------------------------- Icosadome -------------------------------*/

/* Eight (?) vertices of icosadome on unit sphere */

cpoint icosadome_vertices[] = { ZA,ZB,ZC,ZD ,
				YA,      YD ,
				XA,XB,     };

/* Structure for unit icosadome */

triangle icosadome_faces[] = {
    { { YA, XA, YD }, 0.0 },
    { { YA, YD, XB }, 0.0 },
    { { ZA, YA, ZD }, 0.0 },
    { { ZC, YD, ZB }, 0.0 },
    { { XA, YA, ZA }, 0.0 },
    { { YA, XB, ZD }, 0.0 },
    { { YD, XA, ZB }, 0.0 },
    { { YD, ZC, XB }, 0.0 }
};

/* A unit icosadome */

shape_object icod = {
    sizeof(icosadome_faces) / sizeof(icosadome_faces[0]),
    &icosadome_faces[0],
    sizeof(icosadome_vertices) / sizeof(icosadome_vertices[0]),
    &icosadome_vertices[0],
    "Icosadome"
};

/*--- # vertices per tess level for each polygon ----------------------------*/

int  octa_nvert[] = { 6,18, 66,258,1026} ;
int tetra_nvert[] = { 4,10, 34,130, 514} ;
int icosa_nvert[] = {12,42,162,642,2562} ;

int *polygon_nvert[] = {  octa_nvert,
			 tetra_nvert,
			 icosa_nvert };
#endif

