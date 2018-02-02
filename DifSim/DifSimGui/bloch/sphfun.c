/* sphfun.c */

#include "sphfun.h"

#define USE_GSL

void *print_Omega( Omega *omega ) 
{
  printf(" Omega: nomega = %i\n",omega->nomega);

  for(int io=0; io<omega->nomega; io++) {
    printf("omega->spts[io].azi = %f, omega->spts[io].pol = %f\n",
	   omega->spts[io].azi,omega->spts[io].pol) ;
    fflush(stdout) ;
  }
  
}

Omega *alloc_omega( int nazim, int npolr ) 
{
  char *program_name = "alloc_omega";
  int debug = 0;
  Omega *omega = (Omega *)malloc(sizeof(Omega)) ;
  int nomega = nazim * npolr ; 

  spoint *spts = sph_alloc(nomega);

  omega->nazim  = nazim ;
  omega->npolr  = npolr ;
  omega->nomega = nomega ;
  omega->spts   = spts ;

  if(debug) {
    printf("%s: omega->nazim = %i, omega->npolr = %i\n",
	   program_name,omega->nazim,omega->npolr);
  }

  return(omega) ;
}

/**<gen_omega>*******************************************

   Title: GENerate OMEGA

   Purpose: Generate angles

    Input: *nazim = # azimuthal angles
           *nelev = # elevation angles
           setn = set nazim and nelev [0=no,1=yes]

   Output: (omega[nomega].azi,omega[nomega].pol) in radians

   History: Taken from tstfib.cxx, eliminate mask_type 
            04.04.06 - introduce Omega struct (coords.h)   
            04.04.21 - setn
*/
Omega *gen_omega( int *nazim, int *npolr, int setn, int rev )
{
  char *program_name = "gen_omega";
  double dlims[2] ;

  if(setn) {
    SET_LIMS(dlims,1.,100.) ;
    *nazim = (int) pinput("number of azimuthal angles",1.,dlims);
    *npolr = (int) pinput("number of     polar angles",1.,dlims);
  }

  Omega *om = alloc_omega(*nazim,*npolr);

  float *azim = fv_alloc(*nazim);
  float *polr = fv_alloc(*npolr);

  /* Azimuthal angles */ 

  float azmin_def = 0.0  ;
  float azmax_def = 180. ;
  if(*nazim==1) {
    float az_def = 0.0 ; /* so that 1 angle is really no change */
    if (setn) {
      az_def = 90.0;
    }
    SET_LIMS(dlims,0.,180.) ;
    float az = pinput("azimuthal angle",az_def,dlims);
    azim[0] = az;
  } else {
    SET_LIMS(dlims,0.,180.) ;
    float azmin = pinput("minimum azimuthal angle",azmin_def,dlims);
    SET_LIMS(dlims,azmin,dlims[1]) ;
    float azmax = pinput("maximum azimuthal angle",azmax_def,dlims);
    if (rev > 0) {
      vlinsp(azim,azmax,azmin,*nazim);
    }
    else {
      vlinsp(azim,azmin,azmax,*nazim);
    }
  }

  /* Polar angles */ 

  float pomin_def = 0.0  ;
  float pomax_def = 180. ;
  if(*npolr==1) {
    float po_def = 90. ; /* so that 1 angle is really no change */
    SET_LIMS(dlims,0.,180.) ;
    float po = pinput("polar angle",po_def,dlims);
    polr[0] = po;
  } else {
    SET_LIMS(dlims,0.,180.) ;
    float pomin = pinput("minimum polar angle",pomin_def,dlims);
    SET_LIMS(dlims,pomin,dlims[1]) ;
    float pomax = pinput("maximum polar angle",pomax_def,dlims);
    if (rev > 0) {
      vlinsp(polr,pomax,pomin,*npolr);
    }
    else {
      vlinsp(polr,pomin,pomax,*npolr);
    }
  }

  /* Convert separate azim,polr into single Omega (allows linear indexing) */ 
  
  for(int iaz=0; iaz<*nazim; iaz++) {
    for(int ipo=0; ipo<*npolr; ipo++) {
      int iomega = ipo + iaz*(*npolr) ;
      
      om->spts[iomega].rad = 1.0 ;
      om->spts[iomega].azi = azim[iaz] * DEG2RAD ;
      om->spts[iomega].pol = polr[ipo] * DEG2RAD ;
    }
  }

  free(azim) ;
  free(polr) ;
  return (om) ;
}

/**<get_angles_from_omega>********************************************/

void get_angles_from_omega( Omega *omega , 
			    float *azim, float *polr,
			    int   nazim, int   npolr) 
{

  for(int iaz=0; iaz<nazim; iaz++) {
    for(int ipo=0; ipo<npolr; ipo++) {
      int iomega = ipo + iaz*npolr ;
      //rad[] = omega[iomega].rad ;
      azim[iaz] = omega->spts[iomega].azi ;
      polr[ipo] = omega->spts[iomega].pol ;
    }
  }

}

/**<sphrnd>*************************************************

    Title: SPHerical RaNDom sequence

  Purpose: Generate random points uniformly 
           distributed on unit sphere

     Call: float *sphrnd( float *x, float *y, float *z, int n )

    Input: x[n],y[n],z[n] = malloc-ed float array
           n  = number of points on sphere
*/
int sphrnd( float *x, float *y, float *z, int n )
{
  int i;
  float *xold = x;
  float *yold = y;
  float *zold = z;
  double vec[3] ;
  double norm,norminv ;

  srand48(clock_seed());

  for(i=0; i<n; i++) {

    do {
      vec[0] = (1.0-2.0*drand48());
      vec[1] = (1.0-2.0*drand48());
      vec[2] = (1.0-2.0*drand48());
      norm = NORM_VEC3(vec) ;
    } while (norm > 1.0 || norm == 0.0);

    norminv = 1. / norm ;

    xold[i] = (float) norminv * vec[0];
    yold[i] = (float) norminv * vec[1];
    zold[i] = (float) norminv * vec[2];

  }

  return(n) ;			/* for consistency with sphtes, sphlin */
}

/**<sphlin>*************************************************

    Title: SPHerical LINear sequence

  Purpose: Generate linear spaced points mapped onto unit sphere

     Call: float *sphlin( float *x, float *y, float *z, int n )

    Input: x[n],y[n],z[n] = malloc-ed float array
           n  = number of points on sphere
*/
int sphlin( float *x, float *y, float *z, int npolr, int nazim )
{
  int i,j;
  float minpolr = 0.0;
  float minazim = 0.0;
  float maxpolr = PI;
  float maxazim = TWOPI;

  float *polr = fv_alloc(npolr);
  float *azim = fv_alloc(nazim);
  int   ii ;

  int   debug = 0;

  spoint gsph;
  cpoint gcrt;

  vlinsp(polr,minpolr,maxpolr,npolr);
  vlinsp(azim,minazim,maxazim,nazim);

  for(i=0; i<npolr; i++) {
    ii = i*nazim ;
    for(j=0; j<nazim; j++) {

      gsph = newspt(1.,polr[i],azim[j]);
      gcrt = spt2cpt(gsph) ;
      x[ii+j] = gcrt.x ;
      y[ii+j] = gcrt.y ;
      z[ii+j] = gcrt.z ;
      if(debug) {
	printcpt(gcrt,"buf");
	printf("x,y,z=%f,%f,%f\n",x[i],y[i],z[i]);
      }

    }
  }

  free(polr);
  free(azim);

  return(npolr*nazim);
}

/*--------------------<Spherical Harmonic Decomposition>--------------------*/

#ifdef USE_GSL
/**<SphericalHarmonicY.c>****************************************************

     Title: SPHERICAL HARMONIC Y_{l}^{m}

   Purpose: Compute the Spherical Harmonic of order l and degree m
            where l = any integer
	         -l <= m <= l 
	    for a given set of angles (theta,phi) where:

	          0 <= theta <= pi    (polar angle)
		  0 <=  phi  <= 2*pi  (azimuthal angle)

      Call: complex **SphericalHarmonicY(complex *Ylm, 
			                 int l, int m, 
					 float *polr, float *azim, 
					 int nomega)

     Input:    l = order
               m = degree
	       polr[nomega]
	       azim[nomega]

    Output: y[n] = Y_{l}^{m}(theta,phi)

 Algorithm: 

*/
complex *SphericalHarmonicY(complex *Ylm, int l, int m, 
			    float *polr, float *azim, 
			    int nomega)
{
  int     iomega;
  int     am = ABS(m);
  double  Plm,saz,caz ;

  /* Identically equal, since we have array of pairs 
     Omega[i] = {theta[i],phi[i]}
     I just use these for clarity.  Could have just used npts */

  /* Check input */

  if (l < 0) {
    printf("Error! SphericalHarmonicY: \a");
    printf("l=%i must be >= 0 \n",l);
    return NULL;
  }
  if (am > l) {
    printf("Error! SphericalHarmonicY: \a");
    printf("m=%i not in range [0,l] = [%i,%i]\n",m,0,l);
    return NULL;
  }

  /* Over all angles */

  for (iomega=0; iomega<nomega; iomega++) {

      saz = sin(m*azim[iomega]) ;
      caz = cos(m*azim[iomega]) ;

      Plm = gsl_sf_legendre_sphPlm( l, am, cos(polr[iomega]) );

      Ylm[iomega].r = Plm * caz ;
      Ylm[iomega].i = Plm * saz ;

  }


  return Ylm;
}

/**<SphericalHarmonicsYlm.c>**************************************************

     Title: SPHERICAL HARMONICS Y_{l}^{m}

   Purpose: Compute the Spherical Harmonic of order l and degree m
            where l = any integer
                 -l <= m <= l 
	    for a given set of angles (theta,phi)

      Call: complex ***SphericalHarmonicsYlm( complex ***Ylm, int Lmin, int Lmax,
				  float *polr, float *azim, int npts )

     Input:

    Output: y[n] = Y_{l}^{m}(theta,phi)

 Algorithm: 

*/
complex ***SphericalHarmonicsYlm( complex ***Ylm, int Lmin, int Lmax,
				  float *polr, float *azim, int npts )
{
  int ll,mm ;
  int debug = 0;

  for( ll=Lmin ; ll<=Lmax ; ll++ ) {

    for( mm=-ll ; mm<=ll ; mm++ ) {

      if(debug) 
	printf("SphericalHarmonicsYlm: (ll,mm)=(%i,%i), npts = %i\n",
	       ll,mm,npts);

      SphericalHarmonicY( Ylm[ll][mm], ll, mm, polr, azim, npts ); 

    }
  }

  return Ylm ;
}
#endif

/**<SpharmTransformV1.c>****************************************************

     Title: SPHerical HARMonic TRANSFORM, Version 1

   Purpose: Compute the Spherical Harmonic Y_{l}^{m} Transform 
            from min order Lmin to max order Lmax 

      Call: complex **SpharmTransform( complex **Alm, float *data , 
			   int Lmin, int Lmax, 
			   float *pol, float *azi, int nvert)

     Input:    data[npolr=nazim] = data to transform
               Lmax = maximum order
               Mmax = maximum degree
	       npolr = polar angle steps
	       nazim = azimuthal angle steps

    Output: Alm[l][m] = Sum Y_{l}^{m}(theta,phi)
            *coeff is the vector returned to MAKER

   Storage: Complex output stored in matrix dif->shd->calc->Alm[L][M] 
            but also returns magnitude in vector coeff, since this is
	    the way the calling routine DIFF_aniso_meth likes the
	    return to be.  
	    Maybe will change that later ...

 Algorithm: Brute force calculation of coefficients

   History: 18 Sep 2000 conversion from my matlab routine
            06 Jan 2002 now incorporates correct spherical weights
            19 Jan 2002 New indexing of Alm and Ylm by actual
	                (l,m) values (see alm_alloc.c, ylm_alloc.c)
            14 Nov 2002 Conversion from plug_diff version.  
	                No dif_structs and also no coeff vector needed.
*/
complex **SpharmTransformV1( complex **Alm, float *data , 
			      int Lmin, int Lmax, 
			      float *pol, float *azi, int nvert)
{
  char *program_name = "SpharmTransformV1";
  bool debug_local = 0;

  int     ii,ll,mm;
  complex ylm, dys;
  double  scale ;

  complex ***Ylm = ylm_alloc(Lmax,nvert) ;
  float      *wt =  fv_alloc(nvert);

#define SHT(a,b,c,d) (zadd((a),zscal((b),zscal((c),zcnjg((d))))));

  if(debug_local) fprintf(stderr,"%s: Constructing weights ... ",program_name); 
  sphere_voronoi_angles( nvert, pol, azi, &wt ) ;
  if(debug_local) fprintf(stderr,"done\n"); 

  if(debug_local) fprintf(stderr,"%s: Constructing spherical harmonics ... ",program_name); 
  SphericalHarmonicsYlm( Ylm, Lmin, Lmax, pol, azi, nvert );
  if(debug_local) fprintf(stderr,"done\n"); 
  if(debug_local>1) printYlm(Ylm,Lmax,nvert);

  /* Compute data .* conjg(Ylm) .* measure */

  if(debug_local) printf("Inside SpharmTransform...\n");

  /* Do the transform */

  for( ll=Lmin ; ll<=Lmax ; ll++ ) {
    for( mm=-ll ; mm<=ll ; mm++ ) {

      ZCLR(dys);

      for( ii=0 ; ii<nvert ; ii++ ) {

	ylm = Ylm[ll][mm][ii] ;	             /* spherical harmonics */
	dys = SHT(dys,data[ii],wt[ii],ylm) ; /* Spherical Harmonic Transform */

      }

      Alm[ll][mm] = dys;	/* store complex coeff in dif */

      if(debug_local>1) printf("Alm[%i][%i] = (%lf,%lf)\n",ll,mm,dys.r,dys.i);

    }
    if(debug_local) printf("\n");
  }

  if(debug_local>1) TMP_EXIT;

  //ylm_free(Ylm,Lmax) ;
  free(wt);

  return Alm;

}
/**<SpharmTransform.c>****************************************************

     Title: SPHerical HARMonic TRANSFORM

   Purpose: Compute the Spherical Harmonic Y_{l}^{m} Transform 
            from min order Lmin to max order Lmax given data
	    and also the gradient sampling vectors

      Call: complex **SpharmTransform( complex **Alm, float *data, float **gvecs, SHD *shd)

     Input:    data[npolr=nazim] = data to transform
               gvecs[npolr=nazim][3] = gradient sampling directions
	       shd = SHD struct

    Output: Alm[l][m] = Sum Y_{l}^{m}(theta,phi)
            *coeff is the vector returned to MAKER

   Storage: Complex output stored in matrix dif->shd->calc->Alm[L][M] 
            but also returns magnitude in vector coeff, since this is
	    the way the calling routine DIFF_aniso_meth likes the
	    return to be.  
	    Maybe will change that later ...

 Algorithm: Brute force calculation of coefficients

   History: 5Nov08 Conversion from SpharmTransformV1 and SpharmTransformOpt

*/
complex **SpharmTransform( complex **Alm, float *data, float **gvecs, SHD *shd)
{
  char *program_name = "SpharmTransform";
  bool debug_local = 0;

  int     ii,ll,mm;
  complex ylm, dys;
  double  scale ;
  float *gx = fv_alloc(shd->nvert);
  float *gy = fv_alloc(shd->nvert);
  float *gz = fv_alloc(shd->nvert);

#define SHT(a,b,c,d) (zadd((a),zscal((b),zscal((c),zcnjg((d))))));

  for( ii=0 ; ii<shd->nvert ; ii++ ) {
    gx[ii] = gvecs[ii][0];
    gy[ii] = gvecs[ii][1];
    gz[ii] = gvecs[ii][2];
  }

  if(debug_local) fprintf(stderr,"%s: Constructing weights ... ",program_name); 
  cart2polr(shd->rad, shd->pol, shd->azi, gx, gy, gz, shd->nvert );
  sphere_voronoi_angles( shd->nvert, shd->pol, shd->azi, &shd->wt ) ;
  if(debug_local) fprintf(stderr,"done\n"); 

  if(debug_local) fprintf(stderr,"%s: Constructing spherical harmonics ... ",program_name); 
  SphericalHarmonicsYlm( shd->Ylm, shd->Lmin, shd->Lmax, shd->pol, shd->azi, shd->nvert );
  if(debug_local) fprintf(stderr,"done\n"); 
  if(debug_local>1) printYlm(shd->Ylm,shd->Lmax,shd->nvert);

  /* Do the transform:  data .* conjg(Ylm) .* measure */

  for( ll=shd->Lmin ; ll<=shd->Lmax ; ll++ ) {
    for( mm=-ll ; mm<=ll ; mm++ ) {

      ZCLR(dys);
      for( ii=0 ; ii<shd->nvert ; ii++ ) {

	ylm = shd->Ylm[ll][mm][ii] ;	             /* spherical harmonics */
	dys = SHT(dys,data[ii],shd->wt[ii],ylm) ; /* Spherical Harmonic Transform */

      }
      Alm[ll][mm] = dys;	/* store complex coeff in dif */

    }
    if(debug_local) printf("\n");
  }

  free(gx);
  free(gy);
  free(gz);

  return Alm;

}

/**<SpharmTransformOpt.c>****************************************************

     Title: SPHerical HARMonic TRANSFORM - OPTimized

   Purpose: Compute the Spherical Harmonic Y_{l}^{m} Transform 
            from min order Lmin to max order Lmax with 
	    pre-computed spherical harmonic components Ylm and
	    pre-computed weights wt

      Call: complex **SpharmTransformOpt( complex **Alm, float *data, SHD *shd) 

     Input:    data[npolr=nazim] = data to transform
               Lmax = maximum order
               Mmax = maximum degree
	       npolr = polar angle steps
	       nazim = azimuthal angle steps

    Output: Alm[l][m] = Sum Y_{l}^{m}(theta,phi)
            *coeff is the vector returned to MAKER

   Storage: Complex output stored in matrix dif->shd->calc->Alm[L][M] 
            but also returns magnitude in vector coeff, since this is
	    the way the calling routine DIFF_aniso_meth likes the
	    return to be.  
	    Maybe will change that later ...

 Algorithm: Brute force calculation of coefficients

   History: 18 Sep 2000 conversion from my matlab routine
            06 Jan 2002 now incorporates correct spherical weights
            19 Jan 2002 New indexing of Alm and Ylm by actual
	                (l,m) values (see alm_alloc.c, ylm_alloc.c)
            14 Nov 2002 Conversion from plug_diff version.  
	                No dif_structs and also no coeff vector needed.
*/
complex **SpharmTransformOpt( complex **Alm, float *data, SHD *shd) 
{
  int     ii,ll,mm;
  complex ylm, dys;
  double  scale ;

  int     debug = 0;

  /* Compute data .* conjg(Ylm) .* measure */

#define SHT(a,b,c,d) (zadd((a),zscal((b),zscal((c),zcnjg((d))))));

  if(debug) printf("Inside SpharmTransform...\n");

  /* Do the transform */

  for( ll=shd->Lmin ; ll<=shd->Lmax ; ll++ ) {
    for( mm=-ll ; mm<=ll ; mm++ ) {

      ZCLR(dys);

      for( ii=0 ; ii<shd->nvert ; ii++ ) {

	ylm = shd->Ylm[ll][mm][ii] ;              /* spherical harmonics */
	dys = SHT(dys,data[ii],shd->wt[ii],ylm) ; /* Spherical Harmonic Transform */

      }

      Alm[ll][mm] = dys;	/* store complex coeff in dif */

      if(debug) printf("Alm[%i][%i] = (%lf,%lf)\n",ll,mm,dys.r,dys.i);

    }
    if(debug) printf("\n");
  }

  return Alm;
}

/*-----------------------<Spherical Wave Decomposition>----------------------*/

#ifdef USE_GSL
/**<sbesselj>*************************************************

      Title: Spherical BESSEL function (J)

    Purpose: Compute spherical bessel function of 
             order n and argument x[n]:  j(l,x) = j (x)
                                                   l
      Call: float *sbesselj(float *sbj, int l, float *x, int n )

   Returns: double *j[n]

     Notes: This routine only works up through l=5 

*/
float *sbesselj(float *sbj, int l, float *x, int n )
{
  int i;

  for(i=0; i<n; i++) sbj[i] = gsl_sf_bessel_jl(l,(double)x[i]) ;
  return(sbj);

}
#endif

/**<sbesseln>*************************************************

      Title: Spherical BESSEL function (N)

    Purpose: Compute spherical bessel function of 
             order n and argument x[n]:  n(l,x) = n (x)
                                                   l

      Call: double *sbesseln(double *s, int l, float *x, int n )

   Returns: double *s[n]

     Notes: This routine only works up through l=5 

*/
float *sbesseln(float *SphBesN, int l, float *x, int n )
{
  int i,k,Lmax;
  double Nl, NlMin1, NlMin2, sx, dx;

#define Lmax_local 5		/* HACK!! */

  Lmax = Lmax_local ;		/* HACK!! */

  printf("\a\n"
	 "Error in sbesseln: \n"
	 "THIS ROUTINE DOES NOT WORK YET!\n"
	 " ... exiting\n");
  exit(0);

  /* This routine only works up through l=5 */

  if(l>Lmax) {
    printf("\a\n"
	   "Error in sbesseln: \n"
	   "Requested order %i > Max order %i allowed!\n"
	   " ... exiting\n",l,Lmax);
    exit(0);
  }

  if(l==0) {

    for(i=0; i<n; i++) SphBesN[i] = -COSC(x[i]) ;

  } else if(l==1) {

    for(i=0; i<n; i++) {
      if(x[i]==0.) {dx=0.;} else {dx=1./x[i];}
      SphBesN[i] = -dx * (COSC(x[i])-cos(x[i])) ;
    }

  } else {

    for(i=0; i<n; i++) {

      if(x[i]==0.) {

	Nl = 0.0;

      } else {
	
	dx=1./x[i];

	NlMin1 = -COSC(x[i]);
	Nl = dx * (NlMin1-cos(x[i])) ;

	for(k=2; k<=l; k++) {
	  NlMin2 = NlMin1 ;
	  NlMin1 = Nl;
	  Nl = (2*k-1)/x[i]*NlMin1 - NlMin2;
	}
      }

      SphBesN[i] = Nl;
    }

  }

  return(SphBesN) ;
}

#ifdef USE_GSL
/**<SphericalBesseljl.c>**************************************************

     Title: SPHERICAL HARMONICS Y_{l}^{m}

   Purpose: Compute the Spherical Harmonic of order l and degree m
            where l = any integer
                 -l <= m <= l 
	    for a given set of angles (theta,phi)

      Call: float **SphericalBesseljl( float **jl, int Lmin, int Lmax,
			   float *jrad, int nrad )

     Input:

    Output: y[n] = Y_{l}^{m}(theta,phi)

 Algorithm: 

*/
float **SphericalBesseljl( float **jl, int Lmin, int Lmax,
			   float *jrad, int nrad )
{
  int ll ;
  int debug_local = 0;

  for( ll=Lmin ; ll<=Lmax ; ll++ ) sbesselj(jl[ll],ll,jrad,nrad);

  if(debug_local) {
    for( ll=Lmin ; ll<=Lmax ; ll++ ) {
      printf("\nSphericalBesseljl: L = %i\n",ll);
      PRINT_FV(jl[ll],nrad,"jl");
    }
  }

  return jl ;
}
#endif
#ifdef USE_GSL
/**<jl_scale.c>*******************************************************

     Title: Jl SCALE coefficients for SWT

   Purpose: 

      Call: float ***jl_scale( float ***jlscal, int Lmax, 
                               float *jrad, int nrad ) 

     Input: 
    Output: 
 Algorithm: 

***************************************************************************/

float ***jl_scale( float ***jlscal, int Lmax, float *jrad, int nrad ) 
{
  int   ll,mm,rr,am;
  float fm,fp,A,B,b2,j2 ;
  int   debug_local = 0;

  for( ll=0 ; ll<=Lmax ; ll++ ) {

    for( mm=-ll ; mm<=ll ; mm++ ) {

      am = (int) ABS(mm);
      fm = (float)factorial((ll)-(am));
      fp = (float)factorial((ll)+(am));
      if(debug_local) printf("fm=%13.6g, fp=%13.6g\n",fm,fp);

      A = 2.*(2.*ll+1.)*fm/(3.*fp) ;

      for( rr=0 ; rr<nrad ; rr++ ) {

	b2  = jrad[rr]*jrad[rr];
	j2  = gsl_sf_bessel_jl(ll,(double)jrad[rr]) ;
	j2 *= j2 ;
	B   = b2/(b2 - ll*(ll+1))/j2 ;

	if(debug_local) 
	  printf("A=%13.6g, b2=%13.6g, j2=%13.6g, B=%13.6g\n",A,b2,j2,B);

	jlscal[rr][ll][mm] = sqrt(A*B) ;
	jlscal[rr][ll][mm] = jrad[rr]*jrad[rr]/M_PI ;

	if(debug_local) printf("jlscal[%i][%i][%i] = %13.6g\n",
			 rr,ll,mm,jlscal[rr][ll][mm]);
      }
    }
  }

  return(jlscal);
}
#endif

#ifdef USE_GSL
/**<SphericalWavesWlmr.c>**************************************************

     Title: SPHERICAL HARMONICS Y_{l}^{m}

   Purpose: Compute the Spherical Harmonic of order l and degree m
            where l = any integer
                 -l <= m <= l 
	    for a given set of angles (theta,phi)

      Call: complex ****SphericalWavesWlmr( complex ****Wlmr, 
				float *jrad, int nrad,
				int    Lmin, int Lmax,
				float *polr, float *azim, int nvert )

     Input:

    Output: y[n] = Y_{l}^{m}(theta,phi)

 Algorithm: 

*/
complex ****SphericalWavesWlmr( complex ****Wlmr, 
				float *jrad, int nrad,
				int    Lmin, int Lmax,
				float *polr, float *azim, int nvert )
{
  int     ll,mm,rr ;
  complex ***Ylm    =  ylm_alloc(Lmax,nvert) ;
  float   ***jlscal = rlmr_alloc(Lmax,nrad) ;
  float   **jl      =   fm_alloc(Lmax+1,nrad) ;
  float   j2 ;

  int   doscale = 1;
  int   debug_local = 0;

  SphericalHarmonicsYlm( Ylm, Lmin, Lmax, polr, azim, nvert );
  SphericalBesseljl(      jl, Lmin, Lmax,        jrad, nrad  );

  /* Scaling for SWT */

  if(doscale) jl_scale( jlscal, Lmax, jrad, nrad ) ;

  for( rr=0 ; rr<nrad ; rr++ ) {
	
    if(doscale) {
      j2 = jrad[rr]*jrad[rr]/M_PI;
    } else {
      j2 = 1.0;
    }
    for( ll=Lmin ; ll<=Lmax ; ll++ ) {

      for( mm=-ll ; mm<=ll ; mm++ ) {

      if(debug_local) printf("SphericalWavesWlmr: (ll,mm,rr)=(%i,%i,%i), "
		       "nvert = %i, nrad = %i\n",
		       ll,mm,rr,nvert,nrad);

      /* Normalization */

      /* W = jl * Ylm */
      zvscal( Wlmr[rr][ll][mm], zcplx(j2*jl[ll][rr],0.), Ylm[ll][mm], nvert );

      }
    }
  }

  /* The following routines not written yet.
     ylm_free(Ylm,Lmax,nvert) 
     fm_free(jl,nrat) */

    return Wlmr ;
}
#endif

/**<SphWaveTransform.c>****************************************************

     Title: SPHerical HARMonic TRANSFORM

   Purpose: Compute the Spherical Wave Harmonic Transform 
            from min order Lmin to max order Lmax 
	    and between radii Rmin and Rmax

      Call: complex **SphWaveTransform( double *data , float *coeff, DIF_state *dif )

     Input:    data[npolr=nazim] = data to transform
               Lmax = maximum order
               Mmax = maximum degree
	       npolr = polar angle steps
	       nazim = azimuthal angle steps

    Output: Alm[l][m] = Sum Y_{l}^{m}(theta,phi)
            *coeff is the vector returned to MAKER

   Storage: Complex output stored in matrix dif->shd->calc->Alm[L][M] 
            but also returns magnitude in vector coeff, since this is
	    the way the calling routine DIFF_aniso_meth likes the
	    return to be.  
	    Maybe will change that later ...

 Algorithm: Brute force calculation of coefficients

   History: 17 Nov 2002 first version
*/
complex ***SphWaveTransform( complex ***Almr, float **data , 
			     int Lmin, int Lmax, int nvert,
			     float *rad, int nrad,
			     complex ****Wlmr, float *wt )
{
  int     ii,ll,mm,rr;
  complex wlmr, dys;
  double  scale ;

  int     debug_local = 0;

  /* Compute data .* conjg(Wlmr) .* measure */

#define SHT(a,b,c,d) (zadd((a),zscal((b),zscal((c),zcnjg((d))))));

  if(debug_local) {printf("Inside SphWaveTransform...\n");fflush(stdout);}

  /* Do the transform */

  for( rr=0 ; rr<nrad ; rr++ ) {

    for( ll=Lmin ; ll<=Lmax ; ll++ ) {

      for( mm=-ll ; mm<=ll ; mm++ ) {

	ZCLR(dys);

	for( ii=0 ; ii<nvert ; ii++ ) {

	  if(debug_local) {
	    printf("Wlmr[%i][%i][%i][%i] = %13.6g,%13.6g\n",
		   rr,ll,mm,ii,Wlmr[rr][ll][mm][ii].r,Wlmr[rr][ll][mm][ii].i) ;
	    fflush(stdout);}

	  /* Spherical Wave */

	  wlmr = Wlmr[rr][ll][mm][ii] ;

	  /* Spherical Harmonic Transform */

	  dys = SHT(dys,data[rr][ii],wt[ii],wlmr) ;

	}

	Almr[rr][ll][mm] = dys;	/* store complex coeff in dif */

	if(debug_local) printf("Almr[%i][%i][%i] = (%13.6g,%13.6g)\n",
			 rr,ll,mm,dys.r,dys.i);

      }

    }
    if(debug_local) printf("\n");
  }

  return Almr;

}
