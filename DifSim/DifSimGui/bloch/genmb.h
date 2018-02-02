/*<genmsinc.h>***************************************************************/

#ifndef _GENMSINC_HEADER_
#define _GENMSINC_HEADER_

typedef struct {

  int nsl ;			/* number of slices */
  int res ;			/* number of points */

  float seprat ;		/* ration of separation to width */
  float cyc ;			/* sinc cycles */
  float msscale ;		/* scale relative to single sinc needed for the same flip angle */
  float tshift ;		/* time shift */

  short int *msrho, *mstheta ;

} MbStr;

/* constants */
#define GRESMAX 32768     /* number of points grad max */

/* functions */
#define EVENIZE(a) ((int) 2 * ((int)(a/2)))

/* routines */

MbStr *genmb( char *froot ) ;
void init_genmb(MbStr *mstr);
int genmsinc(MbStr *mstr) ;

#endif /* _GENMSINC_HEADER_ */
