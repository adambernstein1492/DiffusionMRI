/*
 * sphtes - SPHere TESsellation
 *  Generate a triangle mesh approximating a sphere by
 *  recursive subdivision. First approximation is a platonic
 *  solid; each level of refinement increases the number of
 *  triangles by a factor of 4.
 *
 *  Input: sphtes [level] [-options]
 *       level is an integer >= 1 setting the tessellation level (default 1).
 *       Options:
 *       -o starts with an octahedron [def] 
 *       -t starts with a tetrahedron 
 *       -i starts with a icosahedron 
 *       -v vertex output only
 *       -f face output only
 *       -m mathematica output
 *       -p prints output
 *       -c  generate triangles with vertices in counterclockwise
 *           order as viewed from the outside in a RHS coordinate system.
 *           The default is clockwise order.
 *
 * 3/24/89       Jon Leech (leech@cs.unc.edu) 
 *    5/93       icosahedral code added by Jim Buddenhagen 
 *               (jb1556@daditz.sbc.com) 
 *
 * 99/9/15  lrf - add vertex structure
 *                add mathematica output
 *                vertex determination (naive brute force)
 * 99/9/16  lrf - faster determination of vertices
 *                using a hashtable (!)
 * 99/9/21  lrf - cp sphere.c sphses.c
 *                ses optimized version of sphere.c
 *                switches from vertex_union to vertex_unhash
 *                at HASHTHRESH < difstr->ntess
 * 99/9/22  lrf - dome mode
 * 99/11/02 lrf - difstr sole input, update difstr->ndif,difstr->nvert
 * 01/12/28 lrf - more normal input for more general use; returns nvert
 *                Gradients no longer include superfluous (0,0,0) at
 *                first points - only actual vertices in g[] arrays!
 *                Also, sphses no longer requires dif_structs.h (good!)
 *                polytope now has allowable values (0,1,2) instead of (1,2,3)
 * 04/04/05 jlr - remove hashtable stuff, vertex_union
 * 04/04/05 lrf - remove DOME
 */
#include "sphtes.h"

/*--------------------------------------------------------------------------*/

int sdebug   = 0;   /* Don't debug */
int doprint  = 0;   /* Don't print out */
int FACEflag = 0;   /* Faces only */
int VERTflag = 0;   /* Vertices only */
int MATHflag = 1;   /* Don't generate Mathematica format output */
int Flatflag = 0;   /* Don't generate per-vertex normals */

/* Forward declarations */
void flip_object( shape_object *obj );
void print_object( shape_object *obj, int level );
void print_triangle(shape_object *obj,int indx) ;
void print_vertex(shape_object *obj,int indx) ;
void print_pt(cpoint *a, char *name) ;

int create_point(cpoint *verts, int *cur_vcnt, int a, int b, int *midptwho,
                 int *midpts)
{
  int i;
  int l = (a>b)?b:a;
  int r = (a>b)?a:b;
  for( i = 0; i < 6; i++ )
    if( midptwho[6*l+i] == r || midptwho[6*l+i] == -1 )
      break;
  if( i == 6 ) assert(false);
  if( midptwho[6*l+i] == r )
  {
    //printf("point between %d and %d exists\n", l, r);
    return midpts[6*l+i];
  }
  //printf("creating point %d between %d and %d\n", *cur_vcnt, l, r);
  midptwho[6*l+i] = r;
  verts[*cur_vcnt].x = (verts[a].x + verts[b].x)*0.5;
  verts[*cur_vcnt].y = (verts[a].y + verts[b].y)*0.5;
  verts[*cur_vcnt].z = (verts[a].z + verts[b].z)*0.5;
  double mag = verts[*cur_vcnt].x * verts[*cur_vcnt].x
	     + verts[*cur_vcnt].y * verts[*cur_vcnt].y
	     + verts[*cur_vcnt].z * verts[*cur_vcnt].z;
  if( mag != 0 )
  {
    mag = 1.0 / sqrt(mag);
    verts[*cur_vcnt].x *= mag;
    verts[*cur_vcnt].y *= mag;
    verts[*cur_vcnt].z *= mag;
  }

  midpts[6*l+i] = *cur_vcnt;
  return (*cur_vcnt)++;
}

int sphtes( float *gx, float *gy, float *gz, int *tris, 
	    int polytope, int ntess )
{
    int    nvert;
    shape_object *old,                /* Default is octahedron */
           *newobj;
    int     ccwflag = 0,        /* Reverse vertex order if true */
            level,              /* Current subdivision level */
            i;

    switch (polytope) {
      case 0:
          old = &oct;
        break;
      case 1:
          old = &tet;
        break;
      case 2:
          old = &ico;
        break;
      default:
        printf("\nsphtes: Invalid polytope type (= %i)!\a ... exiting\n\n",
               polytope);
        exit(0);
    }

    if(doprint) 
      printf("polyvert = %d\n",polygon_nvert[polytope][ntess-1]);

    if (ccwflag) flip_object(old);

    if(doprint) {
      printf("\nPolygon = %s: #faces = %d, #vertices = %d\n",
             old->type,old->nface,old->nvert);
    }

    /* Subdivide each starting triangle (ntess - 1) times */
    int nv = 12;
    for( i = 0; i < ntess; i++ ) nv = (nv-2)*4 + 2;

    cpoint *all_verts = (cpoint*)calloc(sizeof(cpoint), nv);
    for(i=0; i<old->nvert; i++) all_verts[i] = old->vert[i];

    int vcnt     = old->nvert;
    int ontris   = old->nface;
    int *oldtris = (int*)calloc(3*sizeof(int), ontris);
    for(i=0; i<ontris*3; i++) oldtris[i] = icosahedron_ids[i];

    int *newtris;
    for (level = 1; level < ntess; level++) {

        /* Allocate 4* number of face points in the current approximation */

        newtris = (int*)calloc(3* ontris * 4, sizeof(int));
        if (newtris == NULL) {
          printf("sphere: Out of memory on subdivision level %d\n",
                  level);
          exit(-1);
        }

        int ovc = vcnt;

        int *midptwho = (int*)malloc(6*vcnt*sizeof(int));
        int *midpts = (int*)malloc(6*vcnt*sizeof(int));
        memset(midptwho,-1,6*vcnt*sizeof(int));
        memset(midpts,-1,6*vcnt*sizeof(int));

        /* Subdivide each triangle in the old approximation and normalize
         *  the new points thus generated to lie on the surface of the unit
         *  sphere.
         * Each input triangle with vertices labelled [0,1,2] as shown
         *  below will be turned into four new triangles:
         *
         *                      Make new points
         *                          a = (0+2)/2
         *                          b = (0+1)/2
         *                          c = (1+2)/2
         *       1
         *       /\         Normalize a, b, c
         *      /  \
         *    b/____\ c         Construct new triangles
         *    /\    /\              [0,b,a]
         *   /  \  /  \             [b,1,c]
         *  /____\/____\            [a,b,c]
         * 0      a     2           [a,c,2]
         *
         */

        int nont = 0;
        for (i = 0; i < ontris; i++) {
          int *oldt = &oldtris[3*i],
              *newt = &newtris[3*i*4];
          int a, b, c;

          a = create_point(all_verts, &vcnt, oldt[0], oldt[2], midptwho,midpts);
          b = create_point(all_verts, &vcnt, oldt[0], oldt[1], midptwho,midpts);
          c = create_point(all_verts, &vcnt, oldt[1], oldt[2], midptwho,midpts);
            
          newt[0] = oldt[0];
          newt[1] = b;
          newt[2] = a;
          newt+=3; nont++;

          newt[0] = b;
          newt[1] = oldt[1];
          newt[2] = c;
          newt+=3; nont++;

          newt[0] = a;
          newt[1] = b;
          newt[2] = c;
          newt+=3; nont++;

          newt[0] = a;
          newt[1] = c;
          newt[2] = oldt[2]; nont++;
        }
        ontris = nont;

        free(midptwho);
        free(midpts);
        free(oldtris);

        /* Continue subdividing new triangles */
        oldtris = newtris;
        printf("number of vertices: (level %i) %d\n", level+1,vcnt);
    }

    nvert = vcnt;
    memcpy(tris,oldtris,ontris*3*sizeof(int));
    free(oldtris);

    /* Transfer vertices to img */

    /* sphere components */
    for (i=0;i<nvert;i++) { /* nvert = ndif-1 */
      gx[i] = all_verts[i].x ;
      gy[i] = all_verts[i].y ;
      gz[i] = all_verts[i].z ;
    }

    free(all_verts);
    if(doprint) printf("%d vertices, %d faces\n", nvert, old->nface);

    /* Print out sorted object */
    if(doprint) print_object(old, ntess);

    return(nvert);
}

/* 
   Reverse order of points in each triangle 
*/
void flip_object(shape_object *obj)
{
    int i;
    for (i = 0; i < obj->nface; i++) {
        cpoint tmp;
                       tmp = obj->face[i].pt[0];
        obj->face[i].pt[0] = obj->face[i].pt[2];
        obj->face[i].pt[2] = tmp;
    }
}

void print_pt(cpoint *a,char *name)
{
  printf("%s = (%g %g %g)\n", name,a->x, a->y, a->z);
}

/* 
   Write out all triangles in an object 
*/
void print_object(shape_object *obj,int level)
{
    int i;

    /* Spit out coordinates for each face (triangle) */
    if(FACEflag) {
      if(doprint) printf("\n Faces:\n\n");
      for (i = 0; i < obj->nface; i++)
        print_triangle(obj,i);
    }

    /* Spit out coordinates for each vertex */
    if(VERTflag) {
      if(doprint) printf("\n Vertices:\n\n");
      for (i = 0; i < obj->nvert; i++)
        print_vertex(obj,i);
    }

}

/* 
   Output vertex 
*/
void print_vertex(shape_object *obj,int indx)
{
    cpoint *v = &obj->vert[indx];

    if (MATHflag) {
      if(indx==0) printf("{");
      printf("{%g,%g,%g}", v->x, v->y, v->z);
      if(indx==obj->nvert-1)
        printf("}\n");
      else
        printf(",\n");
    } else {
      printf("cpoint\n");
      printf("\t%g %g %g\n", v->x, v->y, v->z);
    }
}

/* 
   Output a triangle 
*/
void print_triangle(shape_object *obj,int indx)
{
    triangle *t = &obj->face[indx];
    int i;

    if (MATHflag) {
      /* printf("MATH triangle\n"); */
      if(indx==0) printf("{");
        printf("{");
        for (i = 0; i < 3; i++) {
            printf("{%g,%g,%g}", t->pt[i].x, t->pt[i].y, t->pt[i].z);
            if(i<2) printf(",");
        }
        if(indx==obj->nface-1)
          printf("}}\n");
        else
          printf("},\n");
    } else {
        /* Modify this to generate your favorite output format
         * Triangle vertices are in t->pt[0..2].{x,y,z}
         * A generic format is provided.
         */
        printf("triangle\n");
        for (i = 0; i < 3; i++)
            printf("\t%g %g %g\n", t->pt[i].x, t->pt[i].y, t->pt[i].z);
    }
}

