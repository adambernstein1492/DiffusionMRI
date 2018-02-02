/**<vsechgp.c>********************************************************
 
     Title: Variable rate gradient SECH pulse for GP
 
   Purpose: Generation hyperbolic secant adiabatic rf pulse
            that returns magnetization back to orginal value

      Call: int vsechgp( int vrgon, int opgp, float gph, char *froot )

     Input: vrgon = turn vrg on [0] 
            opgp  = composite pulse on [0] 
	    gp    = gp (fraction of PI) [0] 
	    outfile root

	    Output: writes short int files with amplitude and phase
	            with extensions .rho and .theta

 Algorithm: This program generates a pair of
	    hyperbolic secant pulses grafted together to
	    generate a completely transparent pulse.
	    Parameters as defined in Silver et al Phys Rev A 31(4 2753 (1985)
            
*/
#include "bloch.h"

int signa_phase( float phs ) 
{
  complex ztmp;
  int     iphase ;

  ztmp.r = 0.;
  ztmp.i = phs ;
  ztmp = zpolr(zexp(ztmp));

  iphase = (int) (ztmp.i * MAX_PG_WAMP / M_PI);

  return iphase ;
}

void init_vsech(VsechStr *vstr)
{
  vstr->Vphase = (short int *) calloc(VRESMAX,sizeof(short int));
  vstr->Vamp   = (short int *) calloc(VRESMAX,sizeof(short int));
  vstr->Vgz    = (short int *) calloc(VRESMAX,sizeof(short int));
	  
  vstr->gslew      = 120.0;	/* readout gradient max slew, mT/m/ms */
  vstr->pwsech     = 15000;	/* Width of sech pulse (us) */
  vstr->ampsech    = 1.0;	/* Amplitude of sech pulse */
  vstr->bwsech     = 2546.48;	/* Bandwidth of sech pulse (Hz) */
  vstr->tagthick   = 250.0;	/* Size of inversion slab (mm) */
  vstr->vsechtrunc = 0.495747;	/* Truncation level of VRG sech pulse (%) */
  vstr->vsechtran  = 1.87;	/* Slice transition of VRG sech pulse (cm) */
  vstr->vsechbeta  = 400.0;	/* Beta of VRG sech pulse (1/s) */
  vstr->vsechmu    = 10.0;	/* Mu of VRG sech pulse */
  vstr->pwvsech    = 30.0;	/* Width of pre-squeezed VRG sech pulse (ms) */
  vstr->vsechdt    = 4.0;	/* Time resolution of VRG sech pulse (us) */
  vstr->gmaxvrg    = 2.2;	/* VRG gradient max amp (G/cm) */
  vstr->vrgdel     = 162;	/* Offset time of VRG pulses. */ 

  vstr->startgzvsech = 0.0;	/* Start amp of Gz of VRG sech pulse (G/cm) */
  vstr->endgzvsech   = 0.0;	/* End   amp of Gz of VRG sech pulse (G/cm) */
  vstr->startgxvsech = 0.0;	/* Start amp of Gx of VRG sech pulse (G/cm) */
  vstr->endgxvsech   = 0.0;	/* End   amp of Gx of VRG sech pulse (G/cm) */

}

/*   VRG SECH pulse
*    inputs (args)
*        beta = parameter for sech (1/s)
*        mu = parameter for sech
*        dur = sech duration (s) PRESQUEEZE
*        slthick = slice thickness (2*beta*mu/grad) (cm)
*        gmaxgcm = grad max(G/cm)
*        slewms = grad slew(G/cm/ms)
*	 dt = time resolution (us)	
*    inputs (CVs)
*    outputs
*        Vamp, Vphase, Vgz
*        res_vsech, pw_vsech
*	 startgzvsech, endgzvsech, midgzvsech
*
*	99.09.21   WML   Converted from Eric's vrgsech.c
*	99.10.03   LRF   accomodate tagging on x
*/
int genvrgsech(VsechStr *vstr, float slthick, float gmaxgcm, float slewms, float dtus)
{
  int i,j,k,ii,n,nmax,count,even,nopt,mode,nrflim,nslewlim,ngradlim;
  long int ier,isize;
  double beta,mu,dur,dt,t,sech;
  double gap,gmax,slew,ablip;
  double coef,twist,pshift,cs,ss,ct,st,dn2,tmp,sl,max;
  double g,a,g0,bt;
  double *amp,*phase,*gz;
  short int itmp;
  int   debug = 0 ;

  gmax = GAMMA*gmaxgcm;            /* convert to s^-1/cm */
  slew = 1000.*GAMMA*slewms;       /* convert to s^-1/cm/s */
  dt   = 1.e-6 * dtus; 		   /* convert to seconds */

  printf("Genvrgsech: slthick = %gcm, gmax = %g\n",slthick, gmax);
  printf("	      slew = %g, dt = %gs\n",slew, dt);

  /* Some initial calculations */
  mode = 0;
  if (vstr->novrg) {
    mu = 10.;
  } else {
    mu = 1.87 * slthick/vstr->vsechtran;  /* fixed side-to-width ratio */
  }

  beta = 400.0 * pow(10.0/mu,0.51); /* keep constant threshold */
  g0   = (2.*beta*mu)/(slthick);    /* Hz/cm */
  if (vstr->novrg) gmax = g0;
  dur  = 2.*acosh(100.*gmax/(vstr->vsechtrunc*g0))/beta;  /* sec */

  if (g0 >= gmax) {		    /* go straight to grad limited */
    nmax = (int) ( (g0/gmax) * dur/dt + 20 );
    mode = 3;
  } else {
    nmax = (int) (dur/dt + 2.);
  }

  /* Temporary float arrays for right half of waves */
  amp   = (double *)calloc((unsigned) nmax/2, sizeof(double));
  phase = (double *)calloc((unsigned) nmax/2, sizeof(double));
  gz    = (double *)calloc((unsigned) nmax/2, sizeof(double));

  /* Calculate waveforms */
  count    = 0;
  nrflim   = 0;
  nslewlim = 0;
  ngradlim = 0;
  bt       = 0.;
  if (mode == 0) {
	  for (a=-0.5*beta*dt; bt<beta*dur/2; ) {  /* RF limited */
		  a   += beta*dt;
		  bt   = 2.*atanh(tan(0.5*a));
		  sech = 1./cosh(bt);
		  g = g0 / sech;
		  if (g > gmax) { mode=2; break; }
		  if (count) if (g-gz[count-1] > slew*dt) { mode=1; break; }
		  amp[count]   = 1.;
		  phase[count] = -1.0 * mu * log(sech);
		  gz[count++]  = g;
	  }
	  nrflim = count;
	  a -= beta*dt;
	  bt = 2.*atanh(tan(0.5*a));
  }
  if (mode == 1) {
	  for (; bt<beta*dur/2; ) { /* Slew limited */
		  g = gz[count-1] + slew*dt;
		  if (g > gmax) { mode=2; break; }
		  bt += beta*dt*g/g0;
		  sech = 1./cosh(bt);
		  amp[count]   = sech * g/g0;
		  phase[count] = -1.0 * mu * log(sech);
		  gz[count++]  = g;
	  }
	  nslewlim = count - nrflim;
  }
  if (mode > 1) {
	  if (mode == 3) bt = -0.5*beta*dt*gmax/g0;
	  for (; bt<beta*dur/2-.51*(beta*dt*gmax/g0); ) { /* Gradient limited */
		  bt += beta*dt*gmax/g0;
		  sech = 1./cosh(bt);
		  amp[count]   = sech * gmax/g0;
		  phase[count] = -1.0 * mu * log(sech);
		  gz[count++]  = gmax;
	  }
	  ngradlim = count - nrflim - nslewlim;
  }
  n = count * 2;
  if(n > VRESMAX)  {
      perror("VRG waveform too long");
      printf("VRG waveform too long\n");
      return(0);
    }

  /* Generate sampled sech and gradient waveforms */

  for (i=0;i<count;i++) {

	  itmp = (int) (MAX_PG_WAMP * amp[i]);
	  vstr->Vamp[i+count]   = 2*(itmp/2);
	  vstr->Vamp[count-i-1] = 2*(itmp/2);

	  itmp = signa_phase( phase[i] );

	  itmp = LIMIT(itmp,-MAX_PG_WAMP,MAX_PG_WAMP);
	  vstr->Vphase[i+count]   = 2*(itmp/2);
	  vstr->Vphase[count-i-1] = 2*(itmp/2);

	  itmp = (int) ( MAX_PG_WAMP * (gz[i]/gmax) );
	  vstr->Vgz[i+count]   = 2*(itmp/2);
	  vstr->Vgz[count-i-1] = 2*(itmp/2);

  }

  vstr->vsechbeta = beta;
  vstr->vsechmu   = mu;
  vstr->pwvsech   = dur*1000.;
  vstr->res_vsech = n;
  vstr->pw_vsech  = (int) (n*dtus);

  if(vstr->tagx) {

    vstr->startgxvsech = (float)vstr->Vgz[0]  *gmax/(MAX_PG_WAMP*GAMMA);
    vstr->endgxvsech   = (float)vstr->Vgz[n-1]*gmax/(MAX_PG_WAMP*GAMMA);
    if(debug) printf("start grad = %f end = %f\n",
		     vstr->startgxvsech,vstr->endgxvsech);

  } else {

    vstr->startgzvsech = (float)vstr->Vgz[0]  *gmax/(MAX_PG_WAMP*GAMMA);
    vstr->endgzvsech   = (float)vstr->Vgz[n-1]*gmax/(MAX_PG_WAMP*GAMMA);
    if(debug) printf("start grad = %f end = %f\n",
		     vstr->startgzvsech,vstr->endgzvsech);

  }

  printf("  vsechgp:  pulse params\n"
	 "\tBaseline grad = %f G/cm\n"
	 "\tAmplification  = %f\n"
	 "\tPulse duration = %f (ms)\n"
	 "\tTruncation = %f %%\n",g0/GAMMA,gmax/g0,n*dt*1000,amp[count-1]*100);

  printf("\tTotal points   = %i:\n"
	 "\t RF limited    = %i\n"
	 "\t Slew limited  = %i\n"
	 "\t Grad limited  = %i\n",n,2*nrflim,2*nslewlim,2*ngradlim);

  /* clean up */
  free(amp) ;
  free(phase) ;
  free(gz) ;

  return(1);

}
/* 
   Title: GENerate GP vsech excitation pulse WAVE

   Purpose: 

   Routines: This calls genvsech()

   Inputs (cv's) opgp = 1: generates composite gp pulse
                      = 0: generates standard vrgsech pulse

   Inputs (args): slthick = slice thickness (2*beta*mu/grad) (cm)
                  gmaxgcm = grad max(G/cm)
		   slewms = grad slew(G/cm/ms)
		     dtus = time resolution (us)	

   Algorithm: A constant phase addition does it

   Note: scaling of phase channel in genvsech.e requires
         definition of pi as VSECH_PI

 */
int gengpwave( VsechStr *vstr, int opgp, float gphase )
{
  int i,j,k ;
  float gpinc ;
  short int itmp,jtmp;
  int   debug = 0 ;

  if(!(genvrgsech(vstr, 
		  vstr->tagthick*0.1, 
		  vstr->gmaxvrg, 
		  vstr->gslew*.01, 
		  vstr->vsechdt))) 
    {
      perror("Oops! genvrgsech failed");
      return(0) ;
    }

  if(debug) printf("gengpwave: opgp = %i, gp = %lf\n",opgp,gphase);

  if(opgp) {	/* Replicate rf train for GP */

    for (i=0; i<vstr->res_vsech; i++) {
      itmp = (int) (gphase*VSECH_PI); 
      jtmp = 2*(itmp/2) ;
      vstr->Vamp[  i+vstr->res_vsech] =  vstr->Vamp[i] ;
      vstr->Vphase[i+vstr->res_vsech] = -vstr->Vphase[i] + jtmp ;
      vstr->Vgz[   i+vstr->res_vsech] =  vstr->Vgz[i] ;
      if(debug) printf("vstr->Vphase[%i] = %i, itmp=%i, jtmp=%i\n",
		       i+vstr->res_vsech,vstr->Vphase[i+vstr->res_vsech],
		       itmp,jtmp);
    }

    if(debug*0) {
      printf("gphase = %f\n",gphase);
      printf("vstr->res_vsech = %i\n",vstr->res_vsech);
      PRINT_IV(vstr->Vphase,vstr->res_vsech,"Vphase") ;
    }

    /* Double pulse length for composite pulse */

    vstr->res_vsech *= 2 ;
    vstr->pw_vsech   = (int) (vstr->res_vsech*vstr->vsechdt);
    vstr->pwvsech   *= 2 ;
    
  }

  return(1);
}

VsechStr *vsechgp( int vrgon, int opgp, float gph, char *froot )
{
	long   int ier,isize;
	char   fout[80],fext[20];
	VsechStr *vstr = (VsechStr *) malloc(sizeof(VsechStr));
	int   debug = 0 ;
	
	/* Initialize vsech struct */

	init_vsech(vstr);

	vstr->novrg = !vrgon;

	if(!vrgon && !opgp) printf("\n\tGenerating standard sech pulse\n\n");
	if(!vrgon &&  opgp) printf("\n\tGenerating composite gp pulse\n\n");
	if( vrgon && !opgp) printf("\n\tGenerating vrg sech pulse\n\n");
	if( vrgon &&  opgp) printf("\n\tGenerating composite vrg gp pulse\n\n");
	if( opgp && debug ) printf("vsechgp: gp = %lf\n",gph);

	/* Generate pulse */

	if(!(gengpwave(vstr,opgp,gph))) {
	  perror("Oops! gengpwave failed");
	  return(0) ;
	}

	printf("GP Pulse duration = %i ms\n",vstr->pw_vsech/1000);

	if(froot!=NULL) {

	  /* Write out files */

	  isize = vstr->res_vsech * sizeof(short int);

	  sprintf(fout,"%s.rho",froot);
	  printf("Writing file %s of size %i\n",fout,vstr->res_vsech);
	  write_iqm(&ier,fout,&isize,(char *)vstr->Vamp);
			
	  if(debug) PRINT_IV(vstr->Vphase,vstr->res_vsech,"vstr->Vphase");

	  sprintf(fout,"%s.theta",froot);
	  printf("Writing file %s of size %i\n",fout,vstr->res_vsech);
	  write_iqm(&ier,fout,&isize,(char *)vstr->Vphase);

	  /* sprintf(fout,"%s.gz",froot);
	     printf("Writing file %s of size %i\n",fout,vstr->res_vsech);
	     write_iqm(&ier,fout,&isize,(char *)vstr->Vgz); */

	}

	return (vstr) ;
}	
