/*<vsechgp.h>***************************************************************/

#ifndef _VSECHGP_HEADER_
#define _VSECHGP_HEADER_

typedef struct {

  float gslew ;			/* readout gradient max slew (mT/m/ms) */
  float vsechbeta;
  float vsechmu;
  float pwvsech;
  int   pw_vsech;
  int   res_vsech;
  int   novrg;
  int   tagx;			/* Slice select on x rather than z */
  float vsechtran  ;		/* Slice transition of VRG sech pulse (cm) */
  float vsechtrunc ;		/* Truncation level of VRG sech pulse (%) */

  float vsechdt ;		/* Time resolution of VRG sech pulse (us) */
  float gmaxvrg ;		/* VRG gradient max amp (G/cm) */
  int   vrgdel ;		/* Offset time of VRG pulses. */

  int   pwsech ;		/* Width of sech pulse (us) */
  float ampsech ;		/* Amplitude of sech pulse */
  float bwsech ;		/* Bandwidth of sech pulse (Hz) */
  float tagthick ;		/* Size of inversion slab (mm) */

  float startgzvsech ;		/* Start amplitude of Gz of VRG sech pulse (G/cm) */
  float endgzvsech ;		/* End amplitude of Gz of VRG sech pulse (G/cm) */
  float startgxvsech ;		/* Start amplitude of Gx of VRG sech pulse (G/cm) */
  float endgxvsech ;		/* End amplitude of Gx of VRG sech pulse (G/cm) */

  short int *Vamp, *Vphase, *Vgz ;

} VsechStr;

#define VRESMAX     15000	/* number of points gp max */
#define MAX_PG_WAMP 32766
#define VSECH_PI    MAX_PG_WAMP	/* PI for vsech pulse */
#define GAMMA       4257.0	/* proton gyromagnetic ratio in Hz/gauss */

#define LIMIT(x, xmin, xmax)   ( (x<xmax)? ((x>xmin)? x:xmin):xmax )

VsechStr *vsechgp( int vrgon, int opgp, float gph, char *froot ) ;
void init_vsech(VsechStr *vstr) ;
int gengpwave( VsechStr *vstr, int gpon, float gph );
int genvrgsech(VsechStr *vstr,float slthick, float gmaxgcm, float slewms, float dtus) ;

#endif /* _VSECHGP_HEADER_ */
