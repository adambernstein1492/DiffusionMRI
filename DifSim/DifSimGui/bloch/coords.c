/**<coords.c>*****************************************************

   Functions for coords

*/
#include "coords.h" 

/* 
   Distance between two Cartesian pts
*/
double cdistance(cpoint cpt1, cpoint cpt2)
{
  double xdiff = cpt1.x - cpt2.x;
  double ydiff = cpt1.y - cpt2.y;
  double zdiff = cpt1.z - cpt2.z;

  return sqrt((xdiff*xdiff)+(ydiff*ydiff)+(zdiff*zdiff));
}

/* Create new Cartesian pt */

cpoint newcpt( float x, float y, float z )
{
  cpoint cpt ;

  cpt.x = x ;
  cpt.y = y ;
  cpt.z = z ;

  return(cpt);
}

/* Create new Spherical pt */

spoint newspt( float rad, float pol, float azi )
{
  spoint spt ;

  spt.rad = rad ;
  spt.pol = pol ;
  spt.azi = azi ;

  return(spt);
}

/* 
   Convert Spherical pt to Cartesian pt
*/
cpoint spt2cpt( spoint spt )
{
  cpoint cpt ;

  cpt.x = spt.rad * sin(spt.pol) * cos(spt.azi);
  cpt.y = spt.rad * sin(spt.pol) * sin(spt.azi);
  cpt.z = spt.rad * cos(spt.pol);

  return(cpt);
}

/* 
   Convert Cartesian pt to Spherical pt
*/
spoint cpt2spt( cpoint cpt )
{
  static spoint spt ;
  cpoint origin ;

  origin  = newcpt(0.0,0.0,0.0);
  spt.rad = cdistance(cpt,origin);
  spt.pol = (spt.rad == 0.0)? 0.0 : acos(cpt.z/spt.rad);
  spt.azi = atan2(cpt.y,cpt.x);

  return(spt);
}


/*
  Title: CARTesian to POLaR conversion
*/
void cart2polr(	float *rad, float *pol, float *azi, 
		float *cx,  float *cy,  float *cz, 
		int nn )
{
  cpoint pcrt ;
  spoint psph ;
  int i;
  int debug = 0;

  for(i=0; i<nn; i++) {
    pcrt = newcpt(cx[i],cy[i],cz[i]);
    psph = cpt2spt( pcrt );
    rad[i] = psph.rad ;	/* better be 1, always! */
    pol[i] = psph.pol ;
    azi[i] = psph.azi ;
    if(debug) {
      char  buf[256];
      sprintf(buf,"Pt %i: ",i);
      printf("\n");
      printcpt(pcrt,buf);
      printspt(psph,buf);
    }
  }

}
