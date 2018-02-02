/*<vector.h>***************************************************************/

#ifndef _VECTOR_HEADER_
#define _VECTOR_HEADER_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdbool.h>

#ifndef TYPEDEF_complex
#define TYPEDEF_complex
typedef struct { float r, i ; } complex ;
#endif

#ifndef TYPEDEF_dcomplex
#define TYPEDEF_dcomplex
typedef struct { double dr, di ; } dcomplex ;
#endif

#ifndef TYPEDEF_byte
#define TYPEDEF_byte
typedef unsigned char byte ;
#endif

#ifndef TWOPI
#define TWOPI  (2.*M_PI)
#endif

#define DEG2RAD (M_PI/180.)
#define RAD2DEG (180./M_PI)

#ifndef MAX
#define MAX(x, y) ( (x) > (y) ? (x) : (y) )
#endif

#ifndef MIN
#define MIN(x, y) ( (x) > (y) ? (y) : (x) )
#endif

#ifndef ABS
#define ABS(x) (((x)<0) ? -(x) : (x))
#endif

#ifndef SWAP
#define SWAP(x,y) \
do{ double temp; temp=(x);(x)=(y);(y)=temp;} while(0) 
#endif

#define BETWEEN(x,a,b)  (((x) >= (a)) && ((x)<=(b)))
#define EVEN(x) (((x) % 2) ? (0) : (1))
#define ODD(x)  (((x) % 2) ? (1) : (0))

#define SIGN2(f,g)  (((g)>0)?(f):(-(f)))
#define SIGN(a)     ((a) >= 0 ? 1 : -1)
#define FMOD(x,y)   ((x)-(y)*floor((x)/(y)))
#define SINC(x) ( ((x)==0.) ? (1) : (sin((x))/(x)) )
#define COSC(x) ( ((x)==0.) ? (1) : (cos((x))/(x)) )
#define POW(x,y) pow((double)(x),(double)(y))

#ifndef SET_LIMS
#define SET_LIMS(lims,xmin,xmax)  ((lims)[0]=(xmin),(lims)[1]=(xmax))
#endif

/* should probably be in timeutils */
long int clock_seed( void );

/******************** special routines for 3-vectors ****************************/

#define DOT_VEC3(x,y) ((x)[0]*(y)[0] + (x)[1]*(y)[1] + (x)[2]*(y)[2])
#define NORM_VEC3(x) sqrt(DOT_VEC3(x,x)) 
#define LOAD_VEC3(v,a,b,c) ((v)[0] = (a), (v)[1] = (b), (v)[2] = (c))
#define UNLOAD_VEC3(a,b,c,v) ((a) = (v)[0], (b)=(v)[1], (c)=(v)[2])
#define MIN_VEC3(v) (MIN(MIN((v)[0],(v)[1]),(v)[2]))
#define MAX_VEC3(v) (MAX(MAX((v)[0],(v)[1]),(v)[2]))

/* Subtraction of two 3-vectors: x-y */
#define SUB_VEC3(z,x,y) ((z)[0]=(x)[0]-(y)[0] , (z)[1]=(x)[1]-(y)[1] , (z)[2]=(x)[2]-(y)[2])

/* Scale vector y and put in x: x = a*y */
#define SCAL_VEC3(x,a,y) \
                 ((x)[0]=(a)*(y)[0],\
		  (x)[1]=(a)*(y)[1],\
		  (x)[2]=(a)*(y)[2])

#define REFLECT_VEC3(x,y) SCAL_VEC3((x),-1.0,(y))

/* Copy vec3 "fromv" into "tov" */
#define COPY_VEC3( tov, fromv ) \
        ((tov)[0] = (fromv)[0], (tov)[1] = (fromv)[1], (tov)[2] = (fromv)[2])

#define VCLIP(v,n,vmin,vmax) \
do{ int i; \
    for(i=0; i<(n); i++) { \
        if(!BETWEEN((v)[i],(vmin),(vmax))) (v)[i]=0.; \
    }\
} while(0)

/************************* some complex math *********************************/

#define SQR(x)   ((x)*(x))
#define CUBE(x)  ((x)*(x)*(x))
#define ZSQR(z)  (SQR(z.r)+SQR(z.i))
#define ZABS(z)  sqrt(ZSQR(z))
#define ZARG(z)  ( ((z).r!=0.0 || (z).i!=0.0) ? atan2((z).i,(z).r) : 0.0 )

#define DZSQR(z)  (SQR(z.dr)+SQR(z.di))
#define DZABS(z)  sqrt(DZSQR(z))
#define DZARG(z)  ( ((z).dr!=0.0 || (z).di!=0.0) ? atan2((z).di,(z).dr) : 0.0 )

/* complex z = 0 */
#define ZCLR(z) ((z).r=0.0,(z).i=0.0)
#define DZCLR(z) ((z).dr=0.0,(z).di=0.0)

/* complex z = y* */
#define ZCNJG(z) ((z).i=-(z).i)
#define DZCNJG(z) ((z).di=-(z).di)

/* complex z = y */
#define ZCPY(z,y) ((z).r=(y).r,(z).i=(y).i)
#define DZCPY(z,y) ((z).dr=(y).dr,(z).di=(y).di)

/* dcomplex dz = fftw_complex fwz */
#define FWZtoDZ(dz,fwz) ((dz).dr=(fwz)[0],(dz).di=(fwz)[1])
#define DZtoFWZ(fwz,dz) ((fwz)[0]=(dz).dr,(fwz)[1]=(dz).di)

/* complex z = fftw_complex fwz */
#define FWZtoZ(z,fwz) ((z).r=(float)(fwz)[0],(z).i=(float)(fwz)[1])
#define ZtoFWZ(fwz,z) ((float)(fwz)[0]=(z).r,(float)(fwz)[1]=(z).i)

/* complex z -> dcomplex dz */
#define ZtoDZ(dz,z) ((dz).dr=(double)(z).r,(dz).di=(double)(z).i)
/* dcomplex dz -> complex z */
#define DZtoZ(z,dz) ((z).r=(float)(dz).dr,(z).i=(float)(dz).di)

/* complex u += v */
#define ZADDTO(u,v) ((u).r += (v).r, (u).i += (v).i)
#define DZADDTO(u,v) ((u).dr += (v).dr, (u).di += (v).di)

/* complex u *= a */
#define ZSCAL(a,z) ((z).r *= (a), (z).i *= (a))
#define DZSCAL(a,z) ((z).dr *= (a), (z).di *= (a))

/* complex u += a*v */
#define ZADDTOSCL(u,a,v) ((u).r += (a)*(v).r, (u).i += (a)*(v).i)
#define DZADDTOSCL(u,a,v) ((u).dr += (a)*(v).dr, (u).di += (a)*(v).di)

/************************* reporting *********************************/

#define PRINT_DZV(v,n,label) \
do{ int i; \
    for(i=0; i<(n); i++) { \
        fprintf(stderr,"%s[%i] = (%13.6g,%13.6g)\n",(label),i,(v)[i].dr,(v)[i].di); \
    }\
} while(0) 

#define PRINT_ZV(v,n,label) \
do{ int i; \
    for(i=0; i<(n); i++) { \
        fprintf(stderr,"%s[%i] = (%13.6g,%13.6g)\n",(label),i,(v)[i].r,(v)[i].i); \
    }\
} while(0) 

#define PRINT_FV(v,n,label) \
do{ int i; \
    for(i=0; i<(n); i++) { \
        fprintf(stderr,"%s[%i] = %13.6g\n",(label),i,(v)[i]); \
    }\
} while(0) 

#define PRINT_DV PRINT_FV

#define PRINT_FM(v,n,m,label) \
do{ int i; \
    for(i=0; i<(n); i++) { \
      for(j=0; j<(m); j++) { \
        fprintf(stderr,"%s[%i][%i] = %13.6g\n",(label),i,j,(v)[i][j]);	\
      }\
    }\
} while(0) 

#define PRINT_IV(v,n,label) \
do{ int i; \
    for(i=0; i<(n); i++) { \
        fprintf(stderr,"%s[%i] = %d\n",(label),i,(v)[i]); \
    }\
} while(0) 

#define PRINT_Z(z,label) \
  do{fprintf(stderr,"%s.r = %f, %s.i = %f\n",(label),(z).r,(label),(z).i);fflush(stdout);} while(0)

#define PRINT_DZ(z,label) \
  do{fprintf(stderr,"%s.dr = %f, %s.di = %f\n",(label),(z).dr,(label),(z).di);} while(0)

#define PRINT_DVAL(str,val) \
   fprintf(stderr,"%s: %13.6g\n",(str),(val))

#define PRINT_VEC3(str,fv) \
   fprintf(stderr,"%s: %13.6g %13.6g %13.6g\n",(str),(fv)[0],(fv)[1],(fv)[2])

/********************** Prototypes *********************************************/

/* Integer routines */

long int factorial(int N) ;

/* Byte vector routines */

int bvfind( int **idx, byte *x, int n );

/* Integer vector routines */

short *svlinsp( short *x, short f1, short f2, int n ) ;
int ivnonzero( int *x, int n );
int ivfind( int **idx, int *x, int n );
int ivmatch( int ival, int *x, int n ) ;
bool ivany( int *x, int ival, int n ) ;
void ivran(int *a, int *amin, int *lamin, int *amax, int *lamax, int n) ;
int ivmnmx(int *x, int n, int *nmin, int *nmax) ;

/* Float vector routines */

float *vabs(float *z, float *x, int n) ;
float *vdif( float *dx, float *x, int n);
float *vgrad( float *gx, float *x, int n);
float *unwrap( float *p, int n, float *mask );
float *vtrn2( float *z, int nx, int ny );
float *vadd( float *z,  float *x,  float a,  float *y,  int N) ;
float *vacc(float *z, float a, float *y, int N) ;
float *vrnd( float *x, int n, char *dist ) ;
float *vsin(float *z, float *x, int n) ;
float *vcos(float *z, float *x, int n) ;
float *vscal(float *z, float a, float *x, float b, int N);
float *vmove( float *z, float *x, int N);
float *vrev( float *z , float *x , int N ) ;
float *vlinsp(float *x,float d1,float d2,int n) ;
float vmsq( float *x, int n ) ;
float vcor(float *x, float *y, int n);
float vstd( float *x, int n ) ;
float vvar( float *x, int n ) ;
float vsum( float *x, int n ) ;
float vmean( float *x, int n ) ;
float vmse(float *x, float *y, int n);
float vmnmx(float *x, int n, int *nmin, int *nmax) ;
void  vran(float *a, float *amin, int *lamin, float *amax, int *lamax, int n) ;
float *vclr(  float *x,  int n );
float *vfill(  float *x, float a,  int n ) ;
float *vexp( float *z,  float *x,  int N) ;
int   vclip(  float *z,  float *x,  int n,
	      float xmin,  float xmax) ;
void vrotxy(float *x, float *y, float phi, int n) ;
float vmedian( float *x, int n );
void vindx(float *x,long int *indx,int n) ;

void print_vran(float *x, int nsamp, char *text);
void print_ivran(int *x, int nsamp, char *text);
void print_dvran(double *x, int nsamp, char *text);

/* Double vector routines */

int dvnan( double *x, int n ) ;
int dvmatch( double ival, double *x, int n ) ;
bool dvany( double *x, double dval, int n );

double *dvsubmat(double *b, double *a, int ma, int na, int rstart, int rend, int cstart, int cend);
double *dvxcor1(double *z, double *x, double *y, int nx, int ny, char *shape);
double *dvcor2(double *c, double *a, double *b, int ma, int na, int mb, int nb, char *shape);

double dvmsq( double *x, int n ) ;
double dvcor( double *x,  double *y,  int n) ;
double dvsum( double *x, int n ) ;
double dvmean( double *x, int n ) ;
double dvmse(double *x, double *y, int n) ;
double dvmedian( double *x, int n );
double dvvar( double *x, int n ) ;
double dvstd( double *x, int n ) ;
double dvmnmx(double *x, int n, int *nmin, int *nmax) ;
void dvran(double *a, double *amin, int *lamin, double *amax, int *lamax, int n) ;

double *dvabs(double *z, double *x, int n) ;
double *dvmul(double *z,double *x,double *y,int n);
double *dvmove( double *z,  double *x, int N) ;
double *dvdif( double *dx, double *x, int n) ;
double *dvlinsp( double *x, double f1, double f2, int n ) ;
double *dvclr(  double *x,  int n );
double *dvfill(  double *x, double a,  int n ) ;
double *dvgrad( double *gx, double *x, int n) ;
double *dvzap2( double *z, int nx, int ny );
double *dvzap3( double *z, int nx, int ny, int nz ) ;
double *dvneg( double *z,  double *x,  int N) ;
double *dvsqrt(double *z, double *x, int N) ;
float  *dvtov(float *z,double *x, int N) ;

/* Complex routines */

complex zrect(complex z) ;
complex zpolr(complex z) ;
complex zcnjg(complex a) ;
complex zmuli(complex a) ;
double  zabs(complex a) ;
complex zexp(complex z) ;
complex zdiv(complex a,complex b) ;
complex zscal(double d, complex z) ;
complex zmul(complex a, complex b) ;
complex zmulc(complex a,complex b) ;
complex zadd(complex a, complex b) ;
complex zsub(complex a, complex b) ;
complex zcplx(double r, double i) ;
complex zadd(complex a, complex b) ;
complex zmean(complex a, complex b) ;

/* Complex vector routines */

float   *zvabs(  float *z, complex *x,  int N) ;
complex *zvadd(complex *z, complex *x, complex a,  complex *y,  int N);
complex *zvcplx( complex *z, float *real, float *imag,  int N );
complex *zvzap( complex *z , int n ) ;
complex *zvzap3( complex *z, int nr, int nc, int nz ) ;
complex *zvzap2( complex *z, int nr, int nc ) ;
complex *zvzap2rc( complex *z, int nr, int nc, int rowcol ) ;
complex *zvpolr(complex *z, complex *x,int N) ;
complex *zvrect(complex *z,complex *x,int N) ;
complex *zvexp(complex *z, complex *x, int n) ;
complex *zvclr(complex *x ,  int N ) ;
complex *zvmove(complex *z, complex *x, int N) ;
complex *zvrev( complex *z, complex *x, int N) ;
dcomplex *zvtodzv( dcomplex *z, complex *x, int N);
complex *zvtori( float *real,  float *imag,complex *z,  int N);
complex *zvcnjg(complex *z, complex *x,  int N) ;
complex  zvcor(complex *x, complex *y,  int n) ;
complex  zvconv(complex *x,complex *y,  int n) ;
complex *zvtrn2( complex *z, int nx, int ny ) ;
complex *zvmul( complex *z,  complex *x,  complex *y,  int n );
complex *zvmulc( complex *z,  complex *x,  complex *y,  int n );
complex *zvdiv( complex *z,  complex *x,  complex *y,  int N) ;
double   zvmnmx(complex *x, int n, int *nmin, int *nmax) ;
void     zvran(complex *a, float *amin, int *lamin, float *amax, int *lamax, int n) ;
complex  zvsum(complex *x,  int n) ;
complex  zvmean( complex *x, int n) ;
double   zvstd( complex *x, int n ) ;
complex *zvscal(complex *z, complex a,  complex *x,  int N) ;
complex zrot(complex z, float theta) ;
complex zexpI(float theta) ;
complex zexp(complex dz) ;
float  zarg(complex z) ;

dcomplex dzvcor( dcomplex *x, dcomplex *y, int n ) ;
dcomplex dzvconv( dcomplex *x, dcomplex *y, int n ) ;

void print_zvran(complex *z,int nsamp, char *text);

/* double complex routines */
double dzarg(dcomplex z) ;
double *dzvarg(double *arg, dcomplex *z, int N) ;
dcomplex dzrect(dcomplex z) ;
dcomplex dzadd(dcomplex a, dcomplex b);
dcomplex dzdiv(dcomplex a, dcomplex b);
double   dzabs(dcomplex a);
dcomplex dzpolr(dcomplex z);
dcomplex dzcnjg(dcomplex a) ;
dcomplex dzmul(dcomplex a, dcomplex b);
dcomplex dzcor(dcomplex a, dcomplex b) ;
dcomplex dzrot(dcomplex dz, double theta) ;
dcomplex dzscal(double doub , dcomplex z ) ;
dcomplex dzexp(dcomplex z) ;
dcomplex dzsub(dcomplex a,dcomplex b) ;
dcomplex dzcplx(double r,double i);
dcomplex dzexpI(double theta) ;
dcomplex *dzvclr(dcomplex *x ,  int N );
dcomplex *dzvmove(dcomplex *z,dcomplex *x,  int N) ;
dcomplex *dzvpolr(dcomplex *z, dcomplex *x, int N) ;
dcomplex *dzvrect(dcomplex *z, dcomplex *x,int N) ;
dcomplex *dzvcplx(dcomplex *z, double *real, double *imag,  int N);
dcomplex *dzvzap( dcomplex *z , int n ) ;
dcomplex *dzvzap2( dcomplex *z, int nx, int ny );
dcomplex *dzvzap3( dcomplex *z, int nr, int nc, int nz ) ;
double *dvadd( double *z,  double *x,  double a,  double *y,  int N) ;
double *dvacc(double *z, double a, double *y, int N) ;
double *dvscal(double *z, double a, double *x, double b, int N) ;
dcomplex *dzvgrad(dcomplex *gx, dcomplex *x, int n) ;
dcomplex *dzvdiv(dcomplex *z,  dcomplex *x,  dcomplex *y,  int N);
dcomplex *dzvscal(dcomplex *z, dcomplex a,  dcomplex *x,  int N) ;
double   *dzvabs( double *z, dcomplex *x, int N) ;
dcomplex *dzvadd(dcomplex *z, dcomplex *x, dcomplex a, dcomplex *y, int N);
dcomplex *dzvrot(dcomplex *zrot, dcomplex *z, double *theta, int sign, int N) ;
dcomplex *dzvtori(double *real, double *imag, dcomplex *z, int N) ;
dcomplex *dzvmul( dcomplex *z,  dcomplex *x,  dcomplex *y,  int n );
dcomplex *dzvmulc( dcomplex *z,  dcomplex *x,  dcomplex *y,  int n ) ;
double *dvadd(double *z, double *x, double a, double *y, int N) ;
dcomplex *dzvsqrt( dcomplex *z, dcomplex *x, int N ) ;
double dzvmnmx(dcomplex *z,int N,int *nmin,int *nmax);
void print_dzvran(dcomplex *z,int nsamp, char *text) ;

/* NULL when numerical value needed */
#define  NULLVAL -666		   
/* NULL when complex numerical value needed */
#define ZNULLVAL zcplx((double)(NULLVAL),(double)(NULLVAL))

/* double routines */
double  *ADD_DVEC3( double *z, double *x , double a , double *y ) ;
double *SCAL_DVEC3( double *z, double a , double *x ) ;
double  *INC_DVEC3( double *z, double a , double *x ) ;

#endif /* _VECTOR_HEADER_ */
