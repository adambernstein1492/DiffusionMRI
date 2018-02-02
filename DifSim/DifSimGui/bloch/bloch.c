/**<bloch>**************************************************

   Title: BLOCH equation simulation

 Purpose: Simulate Bloch equation 

  Author: LR Frank

 History: June 2003
          July 4, 2003 - pwrf is now always the normal sech pulse
	                 if gp is invoked, pwrf is internally doubled.

          Apr 2013 - VTK garbage eliminated in favor of python plotting
	             routine plot_bloch.py

   Notes: This is an upgraded version of the bloch routines that
          were in ~C/bloch.  The major changes are:

	  1) rf pulse design is now incorporated into the main program. 
	     Only vsechgp.c is used because it is the most general and
	     allows standard or variable rate gradient sech generation,
	     and also gp.


*/
#include "bloch.h"

#define tmp_exit(a) exit(a)

int wrtmag( char *froot, double **mxyz, int n) ;

void instructions(), send_help();

int WRITE_MAG = 1;
int WRITE_SIGNA = 0;
int RF_FROM_FILE = 0;
int CL_INPUT = 0;

int main(int argc,char *argv[])
{
  /* I/O related params */

  int   need_help = 0;
  int   nopt = 0;
  int   required_opt = 1;
  
  int   debug = 0;

  double *mx,*my,*mz,mtot;
  double *amp,*phase,*fmphase;
  double *Mpm,*Mfm,*M0,**Mxyz,**R,**Rfm,**E;
  double dt,theta,et1,et2;
  long  int ier,isize;
  short *rfin;
  int    i,j,nrf,nrf2;
  char   *frame = "PM";
  char   fname[80],fnum[80],*fin,*frf=NULL;
  char junk[80];
  VsechStr *vstr ;
  MbStr *mstr ;
    
  int vrgon, opgp ;
  double gph = 0.0;
    
  int rftype ;

  if (argc < required_opt ){
    instructions();
    exit(0);
  }

  /* Parse arguments */

  for (i = 1; i < argc; i++) { /* ------- Options ------- */

    if (!strncmp(argv[i], "-f", 2)) { /* frame */
      if (++i >= argc) { instructions(); exit(0); }
      frame = argv[i];
      nopt++; nopt++;
    }
    if (!strncmp(argv[i], "-rf", 2)) { /* Read in rf file */
      if (++i >= argc) { instructions(); exit(0); }
      RF_FROM_FILE = 1;
      fin = argv[i];
      nopt++; nopt++;
    }
    if (!strncmp(argv[i], "-signa", 2)) { /* write out signa format */
      if (++i >= argc) { instructions(); exit(0); }
      WRITE_SIGNA = 1;
      frf = argv[i];
      nopt++; nopt++;
    }
    if (!strncmp(argv[i], "-nw", 2)) { /* don't write out magnetization */
      if (i >= argc) { instructions(); exit(0); }
      WRITE_MAG = 0;
      nopt++;
    }
    if (!strncmp(argv[i], "-debug", 2)) { /* print out arrays */
      if (i >= argc) { instructions(); exit(0); }
      debug = 1;
      nopt++;
    }
    if (!strncmp(argv[i], "-help", 2)) { /* help */
      if (i >= argc) { instructions(); exit(0); }
      need_help = 1;
      nopt++;
    }
    if (!strncmp(argv[i], "-c", 2)) { /* command line input of parameters */
      if (i >= argc) { instructions(); exit(0); }
      CL_INPUT = 1;
      nopt++;
    }
      
  } /* end of options */

  if (need_help){
    send_help();
    exit(0);
  }
	
  if (argc < required_opt + nopt){
    instructions();
    exit(0);
  }

  /*---------------- begin input section ----------------*/

  /*for (i=0; i<FRAMES->nlabels; i++) fprintf(stderr,"%s[%i] = %i\n",
    FRAMES->name,i,FRAMES->labels[i]) ; */

  if (CL_INPUT) {
    frame = argv[2];
  }

  if( !MATCH_STR(frame,"PM") && !MATCH_STR(frame,"FM") ) {
    fprintf(stderr,"Unknown frame type %s entered:\n",frame);
    char *slims[] = {"PM","FM"};
    frame = sinput("frame","PM",slims,2);
  }
  fprintf(stderr,"frame = %s\n",frame);
	
  double dlims[2] ;
  float maxamp = 0.0 ;
  float pwrf, pwrf_def,beta, mu, b1, freq, t1, t2 ;

  if (CL_INPUT) {
    t1 = 0.001 * strtod(argv[3], NULL);
    t2 = 0.001 * strtod(argv[4], NULL);
    b1 = 0.1468 * strtod(argv[5], NULL); 
    pwrf = 0.001 * strtod(argv[6], NULL);
    freq = (2.*M_PI*1000.0) * strtod(argv[7], NULL);
  } else {

    /* Relaxation parameters */
    SET_LIMS(dlims,0.,2000.) ;
    t1 =  .001 * pinput("T1 (ms)",0.,dlims);
    SET_LIMS(dlims,0.,500.) ;
    t2 =  .001 * pinput("T2 (ms)",t1/10.0/.001,dlims);
    if(t1==0.) t1 = BIGVAL;
    if(t2==0.) t2 = BIGVAL;

    /* imaging parameters */
    SET_LIMS(dlims,1.,5.) ;
    b1 = .1468 * pinput("b1 (Gauss) x .1468",1.0,dlims);
    SET_LIMS(dlims,1.,100.) ;
    if(RF_FROM_FILE) {
      pwrf_def = 3.2;
    } else {
      pwrf_def = 30.;
    }
    pwrf = 0.001 * pinput("pulse duration (ms)",pwrf_def,dlims);
    SET_LIMS(dlims,0.,10.) ;
    freq = (2.*M_PI*1000.0) * pinput("offset frequency (kHz)",0.0,dlims);

  }

  /*---------------- RF Section ----------------*/

  if(RF_FROM_FILE) {    /* Read in rf train from _old_ versions of sech etc. 
			   This will be obsolete perhaps soon ... */
    /* Read in RF.rho */

    rfin = (short *) calloc((unsigned) BIGRF, sizeof(short));
    isize = 0;
    sprintf(fname,"%s.rho",fin);
    read_iqm(&ier,fname,&isize,(char *)rfin);
    if(ier!=0) {
      fprintf(stderr,"\n**** Problem reading %s! ****\n\a",fname);
      instructions();
      exit(0);
    } else {
      fprintf(stderr,"**** read file %s\n",fname);
    }

    /* From input, determine # rf points */

    nrf = isize/sizeof(short);	
    fprintf(stderr,"... determined nrf=%i\n",nrf);

    /* ... then allocate */

    amp     = dv_alloc(nrf);
    phase   = dv_alloc(nrf);
    fmphase = dv_alloc(nrf);

    for (i=0;i<nrf;i++) {
      if (rfin[i]==0) {
	amp[i] = 1.;
      } else {
	amp[i] = (double) rfin[i] * b1 * GAM/MAXINT;
      }
      if(debug) fprintf(stderr,"rho_in[%i]=%i, amp[%i]=%lf\n",i,rfin[i],i,amp[i]);
    }

    /* Read in RF.theta */

    sprintf(fname,"%s.theta",fin);
    if (file_exists(fname)) {
      read_iqm(&ier,fname,&isize,(char *)rfin);
      if(ier!=0) {
	fprintf(stderr,"\n**** Problem reading %s! ****\n\a",fname);
	instructions();
	exit(0);
      } else {
	fprintf(stderr,"**** read file %s\n",fname);
      }

      for (i=0;i<nrf;i++) {
	phase[i] = rfin[i] * INT_TO_RAD ;
	if(debug) fprintf(stderr,"theta_in[%i]=%i, phase[%i]=%lf\n",i,rfin[i],i,phase[i]);
      }

    } else {
      fprintf(stderr,"*** File %s doesn't exist! ... Setting phases to 0.\n",fname);
      for (i=0;i<nrf;i++) phase[i] = 0. ;
    }

  } else {	     /*---------------- Generate RF ----------------*/


    if (CL_INPUT) {
      rftype = (int) strtol(argv[8], NULL,0);
    } else {
      //char *slims[] = {"i","e","r"};
      //char *rftype = sinput("rf type: [i]nversion, (e)xcitation, (r)efocusing","i",slims,3);
      SET_LIMS(dlims,RF_INVERT,RF_REFOCUS) ;
      int rf_def = RF_INVERT ;
      rftype = (int) pinput("rf type: (0) Inversion, (1) excitation, (2) refocusing",rf_def,dlims);
    }

    switch (rftype) {
    case RF_INVERT:

      fprintf(stderr,"Generating inversion pulse ...\n");

      if (CL_INPUT) {
        beta = strtod(argv[9], NULL);
        mu = strtod(argv[10], NULL); 
        vrgon = (int) strtol(argv[11], NULL, 0);
        opgp = (int) strtol(argv[12], NULL, 0);
        if (vrgon) {
          pwrf *= 2; 
          gph = strtol(argv[13], NULL, 0);
        }
      } else {
        /* sech rf parameters */
        SET_LIMS(dlims,1.,1000.) ;
        beta = pinput("beta (s^-1)",400.,dlims);
        SET_LIMS(dlims,1.,20.) ;
        mu = pinput("mu",10.,dlims);

        fprintf(stderr,"\n");
        vrgon = (int) logqyn("Turn on vrg","no");
        opgp = (int) logqyn("Turn on gp","no");

        if(opgp) {
          pwrf *= 2;		/* July 4, 2003! */
          SET_LIMS(dlims,0.,1.) ;
          gph = pinput("gph (fraction of PI)",.5,dlims);
        }

        /* Generate sech rf pulse */
      }

      if( (!WRITE_SIGNA) && (frf != NULL) ) frf = NULL ;
      vstr = vsechgp( vrgon, opgp, gph, frf ) ;

      /* .. done */

      nrf = vstr->res_vsech ;
      fprintf(stderr,"nrf = %i\n",nrf);
      break;

    case RF_EXCITE:
      fprintf(stderr,"Generating excitation pulse ...\n");
      break;
    case RF_REFOCUS:
      fprintf(stderr,"Generating refocussing pulse ...\n");

      /* Generate multiband rf pulse */

      if( (!WRITE_SIGNA) && (frf != NULL) ) frf = NULL ;
      mstr = genmb( frf ) ;

      /* .. done */

      nrf = mstr->res ;
      fprintf(stderr,"nrf = %i\n",nrf);

      break;
    default:
      fprintf(stderr,"\a\n Invalid RF %i requested! ... Exiting.\n",rftype);
      return(-1);
    }

    /* allocate for amp/phase channels */
    amp     = dv_alloc(nrf);
    phase   = dv_alloc(nrf);
    fmphase = dv_alloc(nrf);

    /* scale channels from SIGNA to actual unit for Bloch sim */

    float AmpScale = b1 * GAM/MAXINT ;
    switch (rftype) {
    case RF_INVERT:
      for (i=0; i<nrf; i++) {
	amp[i] = vstr->Vamp[i] * AmpScale ;
	phase[i] = vstr->Vphase[i] * INT_TO_RAD ;
	if(amp[i]>maxamp) maxamp = amp[i] ;
      }
      break;
    case RF_EXCITE:
      fprintf(stderr,"Generating excitation pulse ...\n");
      break;
    case RF_REFOCUS:
      for (i=0; i<nrf; i++) {
	amp[i]   = mstr->msrho[i]   * AmpScale ;
	phase[i] = mstr->mstheta[i] * INT_TO_RAD ;
	phase[i] += -PI2;
	if(amp[i]>maxamp) maxamp = amp[i] ;
      }
      break;
    default:
      fprintf(stderr,"\a\n Invalid RF %i requested! ... Exiting.\n",rftype);
      exit(-1);
    }

  }

  //PRINT_DV(amp,nrf,"amp");

  printRange(amp,nrf,"amp");
  printRange(phase,nrf,"phase");

  /* Allocate for bloch sim */
  E    = dm_alloc(3,3);	
  R    = dm_alloc(3,3);
  Rfm  = dm_alloc(3,3);
  Mpm  = dv_alloc(3); 
  Mfm  = dv_alloc(3); 
  M0   = dv_alloc(3); 
  Mxyz = dm_alloc(nrf,3);

  /* initialize magnetization */
  switch (rftype) {
  case RF_INVERT:
    M0[0] = 0.;
    M0[1] = 0.;
    M0[2] = 1.;
    break;
  case RF_EXCITE:
    M0[0] = 0.;
    M0[1] = 0.;
    M0[2] = 1.;
    break;
  case RF_REFOCUS:
    M0[0] = 0.;
    M0[1] = 1.;
    M0[2] = 0.;
    break;
  }

  fprintf(stderr,"Initial magnetization {mx,my,mz} = {%lf,%lf,%lf}\n",
	  M0[0],M0[1],M0[2]);
  dvmove(Mpm,M0,3);

  if(MATCH_STR(frame,"FM")) {
    if(opgp) {
      for (i=0;i<nrf/2;i++) { 
	fmphase[i]       =  phase[i];
	fmphase[nrf-i-1] = -phase[i];
      }
    } else {
      for (i=0;i<nrf;i++) fmphase[i] = phase[i];
    }
  }

  if (WRITE_MAG) {
    long int iisize = nrf*sizeof(double);
    char fff[80];
    sprintf(fff,"bloch.rho");
    fprintf(stderr,"Writing file %s of size %i\n",fff,nrf);
    write_iqm(&ier,fff,&iisize,(char *)amp);
    sprintf(fff,"bloch.theta");
    fprintf(stderr,"Writing file %s of size %i\n",fff,nrf);
    write_iqm(&ier,fff,&iisize,(char *)phase);
  }

  /* Bloch simulation */

  dt = pwrf/nrf;		/* time step */

  fprintf(stderr,"\tdt = %lf\n",dt);
  fprintf(stderr,"\tfreq = %lf\n",freq);
  fprintf(stderr,"\tb1 = %lf\n",b1);
  fprintf(stderr,"\tt1 = %lf, t2 = %lf\n",t1,t2);

  for (i=0;i<nrf;i++) { 

    Mpm = propagate_rot(Mpm,R,E,amp[i],phase[i],freq,t1,t2,dt);
    if (dvnan(Mpm,3)) {
      perror("Mpm is nan!");
      exit(-1);
    }

    //PRINT_VEC3("Mpm",Mpm);
    mtot = sqrt(dvcor(Mpm,Mpm,3));

    if(!strcmp(frame,"PM")) {
      
      dvmove(Mxyz[i],Mpm,3);
      
    } else if(!strncmp(frame,"FM",2)) {
      
      /* fm rotation matrix */
      
      rotation_matrix_fm(Rfm,fmphase[i]);  
      dvmmul(Mfm,Rfm,Mpm,3,3);
      dvmove(Mxyz[i],Mfm,3);
      
    }

    if(debug) fprintf(stderr,"Mxyz[%i]=(%lf,%lf,%lf)\n",i,Mxyz[i][0],Mxyz[i][1],Mxyz[i][2]);

  }

  /* Print results */
  fprintf(stderr,"  Final magnetization {mx,my,mz} = {%lf,%lf,%lf}\n",
	 Mxyz[nrf-1][0],Mxyz[nrf-1][1],Mxyz[nrf-1][2]);
  fprintf(stderr,"  Total magnetization = %lf\n",mtot);

  /* Write out results */
  char *outname = "Mxyz";
  if(WRITE_MAG) wrtmag(outname,Mxyz,nrf);

  /* Clean up */
  free(amp); 
  free(phase); 
  free(fmphase);

}

int wrtmag( char *froot, double **mxyz, int n)
{
  int   i;
  char  fout[80];
  FILE *ofp;

  sprintf(fout,"%s.dat",froot);
  fprintf(stderr,"Writing file %s of size %i\n",fout,n);

  ofp = fopen(fout, "w");

  if (ofp == NULL) {
    fprintf(stderr, "Can't open output file %s!\n",fout);
    exit(1);
  }

  for(i=1; i<n; i++) {
    fprintf(ofp,"%lf %lf %lf\n",mxyz[i][0],mxyz[i][1],mxyz[i][2]);
  }

  fclose(ofp);
  
  return(1);
}

void instructions()
{
  fprintf(stderr,"\nBloch equation simulation of response to rf pulse 'rf_root'\n"
	 " (rotation operator formulation).\n\n"
	 "Usage: bloch\n"
	 "       RF file root\n"
	 "       Inversion pulse duration (ms) [30] \n"
	 "       beta (s^-1) [400] \n"
	 "       mu [10] \n"
	 "       b1 (Gauss) x .1468 \n"
	 "       T1 (ms) [0 -> inf]\n"
	 "       T2 (ms) [0 -> inf]\n"
	 "       offres freq (khz) \n"
	 "       outfile_root\n\n"
         "Options:\n"
         "\t-f frame [PM/FM]\n"
         "\t         PM = phase modulated\n"
         "\t         FM = frequency modulated\n"
         "\t-signa   Write out signa readable .rho and .theta files\n\n");
  return;
}

void send_help()
{
  instructions() ;
}
