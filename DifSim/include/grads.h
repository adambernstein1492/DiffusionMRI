
/* Tessellation polytopes */

//static char *POLY_TYPES[] = {"Octahedral","Tetrahedral","Icosahedral"} ;

#define POLY_OCT 0
#define POLY_TET 1
#define POLY_ICO 2

//static int POLY_TYPE = POLY_ICO; /* default = icosahedral */

/*--- Gradient definitions -------------------------------------------------*/

typedef struct Grad {		/* Gradient structure */

  int   polytope ;		/* polytope to tesselate */
  int   ntess ;			/* tesselation level */

  float *gx ;			/* gx[nvert] - x grad */
  float *gy ;			/* gy[nvert] - y grad */
  float *gz ;			/* gz[nvert] - z grad */
  int   nvert ;			/* number of grad directions */

  int   *tris ;			/* tris[ntris] - triangles for graphics */
  int   ntris ;			/* number of triangles */

} Grad;

/* pp=0 don't print either grads or tris
   pp=1 print grads
   pp=2 print tris
   pp=3 print grads and tris */

#define PRINT_GRAD(grd,pmode) \
do{ int i;                 \
 if((grd) != NULL) {        \
   fprintf(stderr,           \
	   "polytope = %d\n" \
	   "ntess = %d\n"    \
	   "nvert = %d\n"    \
	   "ntris = %d\n",   \
	   (grd)->polytope,(grd)->ntess,(grd)->nvert,(grd)->ntris) ; \
   if((pmode)==1 || (pmode)==3) {       \
     for(i=0; i<(grd)->nvert; i++) {  \
       fprintf(stderr,"grad[%i] = (%f,%f,%f)\n", \
	       i,(grd)->gx[i],(grd)->gy[i],(grd)->gz[i]); \
     } \
   } \
   if((pmode)==2 || (pmode)==3) { \
     for(i=0; i<(grd)->ntris; i++) { \
       fprintf(stderr, "tris[%i] = (%f,%f,%f)\n", \
	       i,(grd)->tris[i]); \
     } \
   } \
 } else { \
   fprintf(stderr,"PRINT_GRAD:  Grad is NULL!\n") ;\
 }\
} while(0);
