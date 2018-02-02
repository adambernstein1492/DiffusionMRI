/*<alloc.c>***************************************************************

   Title: ALLOCation

*/
#include "alloc.h"

/* ----------------------------- vector allocations --------------------- */

byte *bv_alloc(int n) 
{
  char *program_name = "bv_alloc";

  byte *p = (byte *) calloc((n),sizeof(byte)); 
  if (p == NULL) {
    warnx("\a\n%s: allocation error - %s\n",program_name,strerror(errno));
    fflush(stdout);
    printf("Exiting...\n");
    exit(-1);
  }
  return p;
}

char *cv_alloc(int n) 
{
  char *program_name = "cv_alloc";

  char *p = (char *) calloc((n),sizeof(char)); 
  if (p == NULL) {
    warnx("\a\n%s: allocation error - %s\n",program_name,strerror(errno));
    fflush(stdout);
    printf("Exiting...\n");
    exit(-1);
  }
  return p;
}

short *sv_alloc(int n) 
{
  char *program_name = "sv_alloc";

  short *p = (short *) calloc((n),sizeof(short)); 
  if (p == NULL) {
    warnx("\a\n%s: allocation error - %s\n",program_name,strerror(errno));
    fflush(stdout);
    printf("Exiting...\n");
    exit(-1);
  }
  return p;
}

int *iv_alloc(int n) 
{
  char *program_name = "iv_alloc";

  int *p = (int *) calloc((n),sizeof(int)); 
  if (p == NULL) {
    warnx("\a\n%s: allocation error - %s\n",program_name,strerror(errno));
    fflush(stdout);
    printf("Exiting...\n");
    exit(-1);
  }
  return p;
}

float *fv_alloc(int n) 
{
  char *program_name = "fv_alloc";

  float *p = (float *) calloc((n),sizeof(float)); 
  if (p == NULL) {
    warnx("\a\n%s: allocation error - %s\n",program_name,strerror(errno));
    fflush(stdout);
    printf("Exiting...\n");
    exit(-1);
  }
  return p;
}

double *dv_alloc(int n) 
{
  char *program_name = "dv_alloc";

  double *p = (double *) calloc((n),sizeof(double)); 
  if (p == NULL) {
    warnx("\a\n%s: allocation error - %s\n",program_name,strerror(errno));
    fflush(stdout);
    printf("Exiting...\n");
    exit(-1);
  }
  return p;
}

complex *zv_alloc(int n) 
{
  char *program_name = "zv_alloc";

  complex *p = (complex *) calloc((n),sizeof(complex)); 
  if (p == NULL) {
    warnx("\a\n%s: allocation error - %s\n",program_name,strerror(errno));
    fflush(stdout);
    printf("Exiting...\n");
    exit(-1);
  }
  return p;
}

dcomplex *dzv_alloc(int n) 
{
  char *program_name = "zv_alloc";

  dcomplex *p = (dcomplex *) calloc((n),sizeof(dcomplex)); 
  if (p == NULL) {
    warnx("\a\n%s: allocation error - %s\n",program_name,strerror(errno));
    fflush(stdout);
    printf("Exiting...\n");
    exit(-1);
  }
  return p;
}

/* ----------------------------- matrix allocations --------------------- */

/*<cm_alloc.c>****************************************************************/

char **cm_alloc(int nr, int nc)
{
  char **cm;
  int i;

  cm = (char **)calloc(nr,sizeof(char *));

  for (i=0;i<nr;i++) cm[i] = (char *)calloc(nc,sizeof(char));

  return( cm );
}

/*<sm_alloc.c>****************************************************************/

short **sm_alloc(int nr, int nc)
{
  short **sm;
  int i;

  sm = (short **)calloc(nr,sizeof(short *));

  for (i=0;i<nr;i++) sm[i] = (short *)calloc(nc,sizeof(short));

  return( sm );
}

/*<im_alloc.c>****************************************************************/

int **im_alloc(int nr, int nc)
{
  int **im;
  int i;

  im = (int **)calloc(nr,sizeof(int *));

  for (i=0;i<nr;i++) im[i] = iv_alloc(nc);

  return( im );
}

/*<fm_alloc.c>****************************************************************/

float **fm_alloc(int nr, int nc)
{
  float **fm;
  int i;

  fm = (float **)calloc(nr,sizeof(float *));

  for (i=0;i<nr;i++) fm[i] = fv_alloc(nc);

  return( fm );
}

/*<dm_alloc.c>****************************************************************/

double **dm_alloc(int nr, int nc)
{
  double **dm;
  int i;

  dm = (double **)calloc(nr,sizeof(double *));

  for (i=0;i<nr;i++) dm[i] = dv_alloc(nc);

  return( dm );
}

/*<zm_alloc.c>****************************************************************/

complex **zm_alloc(int nr, int nc)
{
  complex **zm;
  int i;

  zm = (complex **)calloc(nr,sizeof(complex *));

  for (i=0;i<nr;i++) zm[i] = zv_alloc(nc);

  return( zm );
}

/*<dzm_alloc.c>***************************************************************/

dcomplex **dzm_alloc(int nr, int nc)
{
  dcomplex **dzm;
  int i;

  dzm = (dcomplex **)calloc(nr,sizeof(dcomplex *));

  for (i=0;i<nr;i++) dzm[i] = dzv_alloc(nc);

  return( dzm );
}


/*<_mm_free.c>****************************************************************/

void _mm_free( char **mm, int nr, int nc )
{
  int i;

  for (i=0;i<nr;i++) free(mm[i]) ;

  free((char *)mm);

}

/* ----------------------------- tensor allocations --------------------- */

/*<st_alloc.c>****************************************************************/

short ***st_alloc( int ns, int nr, int nc )
{
  short ***st;
  int i;

  st = (short ***)calloc(ns,sizeof(short **));

  for (i=0;i<ns;i++) st[i] = sm_alloc(nr,nc) ;

  return( st );
}

/*<it_alloc.c>****************************************************************/

int ***it_alloc( int ns, int nr, int nc )
{
  int ***it;
  int i;

  it = (int ***)calloc(ns,sizeof(int **));

  for (i=0;i<ns;i++) it[i] = im_alloc(nr,nc) ;

  return( it );
}

/*<f3_alloc.c>****************************************************************/

float ***f3_alloc( int ns, int nr, int nc )
{
  float ***ft;
  int i;

  ft = (float ***)calloc(ns,sizeof(float **));

  for (i=0;i<ns;i++) ft[i] = fm_alloc(nr,nc) ;

  return( ft );
}

/*<d3_alloc.c>****************************************************************/

double ***d3_alloc( int ns, int nr, int nc )
{
  double ***dt;
  int i;

  dt = (double ***)calloc(ns,sizeof(double **));

  for (i=0;i<ns;i++) dt[i] = dm_alloc(nr,nc) ;

  return( dt );
}

/*<zt_alloc.c>****************************************************************/

complex ***zt_alloc( int ns, int nr, int nc )
{
  complex ***zt;
  int i;

  zt = (complex ***)calloc(ns,sizeof(complex **));

  for (i=0;i<ns;i++) zt[i] = zm_alloc(nr,nc) ;

  return( zt );
}

/*<dzt_alloc.c>****************************************************************/

dcomplex ***dzt_alloc( int ns, int nr, int nc )
{
  dcomplex ***dzt;
  int i;

  dzt = (dcomplex ***)calloc(ns,sizeof(dcomplex **));

  for (i=0;i<ns;i++) dzt[i] = dzm_alloc(nr,nc) ;

  return( dzt );
}

/*<d4_alloc.c>****************************************************************/

double ****d4_alloc( int nm, int ns, int nr, int nc )
{
  double ****d4;
  int i;

  d4 = (double ****)calloc(nm,sizeof(double ***));

  for (i=0;i<nm;i++) d4[i] = d3_alloc(ns,nr,nc) ;

  return( d4 );
}

/*<d4_free.c>****************************************************************/

void d4_free( double ****d4, int nm, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<nm;i++) dt_free(d4[i],ns,nr,nc);

  free((char *)d4);

}

/*<z4_alloc.c>****************************************************************/

complex ****z4_alloc( int nm, int ns, int nr, int nc )
{
  complex ****z4;
  int i;

  z4 = (complex ****)calloc(nm,sizeof(complex ***));

  for (i=0;i<nm;i++) z4[i] = zt_alloc(ns,nr,nc) ;

  return( z4 );
}

/*<dz4_alloc.c>****************************************************************/

dcomplex ****dz4_alloc( int nm, int ns, int nr, int nc )
{
  dcomplex ****dz4;
  int i;

  dz4 = (dcomplex ****)calloc(nm,sizeof(dcomplex ***));

  for (i=0;i<nm;i++) dz4[i] = dzt_alloc(ns,nr,nc) ;

  return( dz4 );
}

/*<d5_alloc.c>****************************************************************/

double *****d5_alloc( int nq, int nm, int ns, int nr, int nc )
{
  double *****d5;
  int i;

  d5 = (double *****)calloc(nq,sizeof(double ****));

  for (i=0;i<nq;i++) d5[i] = d4_alloc(nm,ns,nr,nc) ;

  return( d5 );
}

/*<d5_free.c>****************************************************************/

void d5_free( double *****d5, int nq, int nm, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<nq;i++) d4_free(d5[i],nm,ns,nr,nc);

  free((char *)d5);

}

/*<z5_alloc.c>****************************************************************/

complex *****z5_alloc( int nq, int nm, int ns, int nr, int nc )
{
  complex *****z5;
  int i;

  z5 = (complex *****)calloc(nq,sizeof(complex ****));

  for (i=0;i<nq;i++) z5[i] = z4_alloc(nm,ns,nr,nc) ;

  return( z5 );
}

/*<dz5_alloc.c>****************************************************************/

dcomplex *****dz5_alloc( int nq, int nm, int ns, int nr, int nc )
{
  dcomplex *****dz5;
  int i;

  dz5 = (dcomplex *****)calloc(nq,sizeof(dcomplex ****));

  for (i=0;i<nq;i++) dz5[i] = dz4_alloc(nm,ns,nr,nc) ;

  return( dz5 );
}

/*<f3_free.c>****************************************************************/

void f3_free( float ***ft, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<ns;i++) fm_free(ft[i],nr,nc);

  free((char *)ft);

}

/*<d3_free.c>****************************************************************/

void d3_free( double ***dt, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<ns;i++) dm_free(dt[i],nr,nc);

  free((char *)dt);

}

/*<zt_free.c>****************************************************************/

void zt_free( complex ***zt, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<ns;i++) zm_free(zt[i],nr,nc);

  free((char *)zt);

}

/*<dzt_free.c>****************************************************************/

void dzt_free( dcomplex ***dzt, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<ns;i++) dzm_free(dzt[i],nr,nc);

  free((char *)dzt);

}

/*<z4_free.c>****************************************************************/

void z4_free( complex ****z4, int nm, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<nm;i++) zt_free(z4[i],ns,nr,nc);

  free((char *)z4);

}

/*<dz4_free.c>****************************************************************/

void dz4_free( dcomplex ****dz4, int nm, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<nm;i++) dzt_free(dz4[i],ns,nr,nc);

  free((char *)dz4);

}

/*<z5_free.c>****************************************************************/

void z5_free( complex *****z5, int nq, int nm, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<nq;i++) z4_free(z5[i],nm,ns,nr,nc);

  free((char *)z5);

}

/*<dz5_free.c>****************************************************************/

void dz5_free( dcomplex *****dz5, int nq, int nm, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<nq;i++) dz4_free(dz5[i],nm,ns,nr,nc);

  free((char *)dz5);

}

/*<_tt_free.c>****************************************************************/

void _tt_free( char ***tt, int ns, int nr, int nc )
{
  int i;

  for (i=0;i<ns;i++) mm_free(tt[i],nr,nc) ;

  free((char *)tt);
}

/* ----------------------------- ntuple allocations --------------------- */

/*<ftupleN_alloc.c>**********************************************************

  Title: Float nTUPLE in N dimensions

  x[ns][nr][nc][ncomp]

*/
float ****ftupleN_alloc( int ns, int nr, int nc, int ncomp )
{
  float ****ft;
  int i;

  ft = (float ****)calloc(ns,sizeof(float ***));

  for (i=0;i<ns;i++) ft[i] = f3_alloc(nr,nc,ncomp) ;

  return( ft );
}

/*<dtupleN_alloc.c>**********************************************************

  Title: Double nTUPLE in N dimensions

  x[ns][nr][nc][ncomp]

*/
double ****dtupleN_alloc( int ns, int nr, int nc, int ncomp )
{
  double ****dt;
  int i;

  dt = (double ****)calloc(ns,sizeof(double ***));

  for (i=0;i<ns;i++) dt[i] = d3_alloc(nr,nc,ncomp) ;

  return( dt );
}

/*<ztupleN_alloc.c>**********************************************************

  Title: Complex nTUPLE in N dimensions

  x[ns][nr][nc][ncomp]

*/
complex ****ztupleN_alloc( int ns, int nr, int nc, int ncomp )
{
  complex ****zt;
  int i;

  zt = (complex ****)calloc(ns,sizeof(complex ***));

  for (i=0;i<ns;i++) zt[i] = zt_alloc(nr,nc,ncomp) ;

  return( zt );
}

/*<ftuple3_free.c>**********************************************************

  Title: Float nTUPLE in 3 dimensions

  x[ns][nr][nc][ncomp]

*/
void _tupleN_free( char ****ft, int ns, int nr, int nc, int ncomp )
{
  int i;

  for (i=0;i<ns;i++) tt_free(ft[i],nr,nc,ncomp) ;

  free((char *)ft);
}

/*-------------------------------------------------------------------

  Can the following be made obsolete by using the _mm_free above?

-------------------------------------------------------------------*/

/**<fm_free.c>********************************************************

       Title:  Float Matrix FREE

     Purpose:  Deallocate float matrix d[m][n]

*/
void fm_free( float **fm, int nrow, int ncol )
{
  static int ir ;

  if( fm != NULL ){
    for( ir=0 ; ir < nrow ; ir++ ) if( fm[ir] != NULL ) free(fm[ir]) ;
    free(fm) ;
  }

}

/**<dm_free.c>********************************************************

       Title:  Double Matrix FREE

     Purpose:  Deallocate double matrix d[m][n]

*/
void dm_free( double **dm, int nrow, int ncol )
{
  static int ir ;

  if( dm != NULL ){
    for( ir=0 ; ir < nrow ; ir++ ) if( dm[ir] != NULL ) free(dm[ir]) ;
    free(dm) ;
  }

}

/**<im_free.c>********************************************************

       Title:  Int Matrix FREE

     Purpose:  Deallocate int matrix d[m][n]

*/
void im_free( int **im, int nrow, int ncol )
{
  static int ir ;

  if( im != NULL ){
    for( ir=0 ; ir < nrow ; ir++ ) if( im[ir] != NULL ) free(im[ir]) ;
    free(im) ;
  }

}

/**<sm_free.c>********************************************************

       Title:  Short Matrix FREE

     Purpose:  Deallocate short matrix d[m][n]

*/
void sm_free( short **sm, int nrow, int ncol )
{
  static int ir ;

  if( sm != NULL ){
    for( ir=0 ; ir < nrow ; ir++ ) if( sm[ir] != NULL ) free(sm[ir]) ;
    free(sm) ;
  }

}

/**<cm_free.c>********************************************************

       Title:  Char Matrix FREE

     Purpose:  Deallocate char matrix d[m][n]

*/
void cm_free( char **cm, int nrow, int ncol )
{
  static int ir ;

  if( cm != NULL ){
    for( ir=0 ; ir < nrow ; ir++ ) if( cm[ir] != NULL ) free(cm[ir]) ;
    free(cm) ;
  }

}

/**<zm_free.c>********************************************************

       Title:  Complex Matrix FREE

     Purpose:  Deallocate complex matrix d[m][n]

*/
void zm_free( complex **zm, int nrow, int ncol )
{
  static int ir ;

  if( zm != NULL ){
    for( ir=0 ; ir < nrow ; ir++ ) if( zm[ir] != NULL ) free(zm[ir]) ;
    free(zm) ;
  }

}

/**<dzm_free.c>********************************************************

       Title:  Double Complex Matrix FREE

     Purpose:  Deallocate dcomplex matrix d[m][n]

*/
void dzm_free( dcomplex **dzm, int nrow, int ncol )
{
  static int ir ;

  if( dzm != NULL ){
    for( ir=0 ; ir < nrow ; ir++ ) if( dzm[ir] != NULL ) free(dzm[ir]) ;
    free((char *)dzm);
  }

}
