/**<genmb.c>********************************************************
 
     Title: GENerate MultiBand pulses
 
   Purpose: Multiband excitation and refocussing pulses

      Call: MbStr *genmb( int mbtype, char *froot )

     Input: mbtype = type of pulse
            froot = output filename for magnitude and phase files

    Output: froot.rho = short int amplitude file
            froot.theta = short int phase file
            msscale - scale relative to single sinc needed for the same flip angle

 Algorithm: mysterious

    Author: Eric Wong, UCSD 2011.02.15
            LRF added time shift 2013.03.20
            
*/
#include "bloch.h"

void init_genmb(MbStr *mstr)
{
  mstr->nsl = 1;	 /* number of simultaneously excited slices */
  mstr->res = 3200;	 /* number of points in pulse */

  mstr->cyc = 2.0;	 /* number of sinc cycles */
  mstr->msscale = 1.0;   /* scale relative to single sinc needed for the same flip angle */
  mstr->tshift = 0.0;	 /* time shift */

  mstr->msrho = (short int *) calloc(GRESMAX,sizeof(short int));
  mstr->mstheta   = (short int *) calloc(GRESMAX,sizeof(short int));

}
/**<genmsinc.c>********************************************************
 
     Title: GENerate Multiband SINC pulse
 
   Purpose: Multislice hamming windowed sinc pulse

      Call: int vsechgp( int vrgon, int opgp, float gph, char *froot )

     Input: nsl = number of slices
            seprat = ratio of separation to width
	    cyc = sinc cycles
	    res = number of points


    Output: msrho.rho = short int amplitude file
            mstheta.theta = short int phase file
            msscale - scale relative to single sinc needed for the same flip angle

 Algorithm: 

    Author: Eric Wong, UCSD 2011.02.15
            LRF added time shift 2013.03.20
            
*/
int genmsinc(MbStr *mstr)
{
  int i,j, rhomax;
  double mid, hamscale, sincscale, rfi, rfq, ham, tsinc, toff, sinc, p;
  double tsinc_min, tsinc_max;

  double scale[16] = {1.0, 2.0, 2.236013, 2.472041, 2.736891, 3.042766, 3.102312, 3.224903,
		      3.337204, 3.685643, 3.788609, 4.007450, 4.182226, 4.478803, 4.608991, 4.814810};

  double phase[16][16] = { 
    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 3.141593, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 7.013606, 4.602434, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 3.874710, 5.940306, 6.196789, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 3.778016, 5.334583, 0.871946, 6.754247, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 2.005490, 1.674044, 5.011975, 5.736084, 4.123246, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 3.001719, 5.997996, 5.909316, 2.623899, 2.528148, 2.440364, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 1.036132, 3.414200, 3.777500, 3.215258, 1.755637, 4.554525, 2.466791, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 1.250064, 1.782967, 3.557958, 0.739262, 3.318675, 7.578945, 0.521146, 5.331898, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 4.418133, 2.359748, 6.960038, 2.252913, 3.471862, 3.039503, 3.973559, 1.191839, 2.510316, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 5.041031, 4.284611, 3.001074, 5.764831, 4.295455, 6.338741, 4.213457, 6.039541, 1.077508, 2.759304, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 2.754564, 5.490923, 4.446980, 0.231401, 2.499357, 3.538964, 2.931499, 2.759175, 5.375936, 4.553648, 3.479277, 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.603323, 0.008577, 4.179383, 4.360830, 4.837402, 0.816330, 5.995430, 4.150257, 6.699861, 7.803329, 4.517049, 1.729496, 0.0, 0.0, 0.0 },
    { 0.0, 3.996713, 0.830202, 5.712114, 3.837526, 0.084402, 1.685373, 5.328242, 6.519758, 6.788787, 1.355796, 4.024798, 4.482871, 4.083972, 0.0, 0.0 },
    { 0.0, 4.126237, 2.265897, 0.957052, 4.602863, 0.815001, 3.474948, 7.260076, 1.448788, 7.475644, 0.148134, 0.939164, 2.530914, 3.611766, 4.800774, 0.0 },
    { 0.0, 4.359012, 3.509978, 4.409574, 1.749566, 3.356647, 2.061054, 5.948342, 3.000272, 2.822407, 0.626690, 2.768343, 3.874516, 4.172815, 4.224429, 5.941258 }
  };

  mstr->msscale = 1.01 * scale[mstr->nsl-1]; /* 1.01 is to ensure no overranging ECW */

  /* Calculate waveforms */

  mid = 0.5*(mstr->res-1);
  hamscale = M_PI / (mstr->res/2);
  sincscale = 2 * M_PI * mstr->cyc / (mstr->res/2);

  rhomax = 0;
  tsinc_min = 0.;
  tsinc_max = 0.;
  for (i=0;i<mstr->res;i++) {
    rfi = 0.0;
    rfq = 0.0;
    ham = 0.54 + 0.46 * cos((i-mid) * hamscale);

    for (j=0;j<mstr->nsl;j++) {
      toff = j*mstr->tshift*mstr->res;
      tsinc = (i-mid+toff) * sincscale;
      sinc = sin(tsinc) / tsinc;
      p = 2 * mstr->seprat * (j - 0.5*(mstr->nsl-1)) * tsinc + phase[mstr->nsl-1][j];
      rfi += cos(p);
      rfq += sin(p);
    }

    mstr->msrho[i]   = EVENIZE(MAX_PG_WAMP * sinc * ham * sqrt(rfi*rfi+rfq*rfq) / mstr->msscale);
    mstr->mstheta[i] = EVENIZE(MAX_PG_WAMP * atan2(rfq,rfi) / M_PI);

    //fprintf(stderr,"mstr->msrho[%i] = %i, mstr->mstheta[%i] = %i\n",i,mstr->msrho[i],i,mstr->mstheta[i]);

    rhomax = MAX(rhomax,mstr->msrho[i]);
  }

  //fprintf(stderr,"\n\tgenmsinc: tsinc_min=%f, tsinc_max=%f, rhomax = %i\n",tsinc_min,tsinc_max,rhomax);

  return 1;

}

MbStr *genmb( char *froot )
{
	long   int ier,isize;
	char   fout[80],fext[20];
	MbStr *mstr = (MbStr *) malloc(sizeof(MbStr));

	bool debug = 0 ;
	
	/* Initialize mband struct */

	fprintf(stderr,"Initializing genmb ... ");
	init_genmb(mstr);
	fprintf(stderr,"done\n");

	/* Generate pulse */

	fprintf(stderr,"Generating msinc ... ");
	if(!genmsinc(mstr)) {
	  perror("Oops! genmsinc failed");
	  return(0) ;
	} else {
	  fprintf(stderr,"done\n");
	}

	if(froot!=NULL) {	/* Write out files */

	  isize = mstr->res * sizeof(short int);

	  if(debug) {
	    PRINT_IV(mstr->msrho,mstr->res,"mstr->msrho");
	    PRINT_IV(mstr->mstheta,mstr->res,"mstr->mstheta");
	  }

	  sprintf(fout,"%s.rho",froot);
	  fprintf(stderr,"Writing file %s of size %i\n",fout,mstr->res);
	  write_iqm(&ier,fout,&isize,(char *)mstr->msrho);
			
	  sprintf(fout,"%s.theta",froot);
	  fprintf(stderr,"Writing file %s of size %i\n",fout,mstr->res);
	  write_iqm(&ier,fout,&isize,(char *)mstr->mstheta);

	}

	return (mstr) ;
}	
