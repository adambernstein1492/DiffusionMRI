/**<alm.c>*******************************************************

     Title: spherical harmonic transform coefficients ALM

   Purpose: Routines for manipulation of the spherical 
            harmonic transform coefficients Alm

    Author: Lawrence R. Frank

   History: 19-Jan-02 initial version (conversion from matlab)
   
     Usage: Routines called "alm2something" take the complex
            array "alm" as input and convert to "something"

    Notes: Index into Alm is now just by actual (l,m) value 
           (see comments in alm_alloc)

  History: This is alm.c from ~/Apps/afni/diff/diffuse, with
           anything using the plug_diff.c structs taken out.
	   ylm_scale put here (in shd.c in plug_diff)
*/
#include "alm.h"

/**<numAlm.c>*******************************************************

     Title: NUMber of ALM coefficients

   Purpose: quickly determine nl*nm

      Call: int numAlm( int Lmax )

     Input: int Lmax = maximum order

    Output: int Nalm = number of alm coeffs

     Notes: There must be an analytical expression for this ...

***********************************************************************/

int numAlm( int Lmax )
{
  int i,j;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;
  int nn = 0;

  for (i=0; i<nl; i++) {
    for (j=-i; j<=i; j++) nn++ ;
  }

  return nn ;
}

/**<ylm_scale.c>*******************************************************

     Title: Ylm SCALE coefficients for SHT

   Purpose: 

      Call: int ylm_scale( int Lmin, int Lmax, float **Ylm_scale ) 

     Input: 
    Output: 

 Algorithm: 

  Return value is 0 if an error transpired, or 1 if all went OK.

***************************************************************************/

float **ylm_scale( float **Ylm_scale, int Lmin, int Lmax) 
{
  int ll,mm,am;
  float ffm,ffp ;
  int debug = 0;

  for( ll=Lmin ; ll<=Lmax ; ll++ ) {
    for( mm=-ll ; mm<=ll ; mm++ ) {
      am = (int) ABS(mm);
      ffm = (float)factorial((ll)-(am));
      ffp = (float)factorial((ll)+(am));
      /*Ylm_scale[ll][mm] = YLM_COEF(ll,am) ;*/
      Ylm_scale[ll][mm] = ((2*(ll)+1)/(4.*PI))*(ffm/ffp);
      if(debug) printf("Ylm_scale[%i][%i] = %g\n",
		       ll,mm,Ylm_scale[ll][mm]);
    }
  }

  return(Ylm_scale);
}

/**<almr_alloc.c>*******************************************************

     Title: ALMR ALLOCation

      Call: int almr_alloc( int Lmax, int npts )

   Purpose: Allocate space for Almr

     Input: int Lmax = maximum L value

    Output: Almr[0,...,nr][0,..,nl][-nl,..,nl]

     Notes: This routine allocated almr so that indexing into
            Almr can be done in the most obvious manner, ie
	    Y_{l}{m}[k] = Almr[l][m][k], eg, Y_{2}{-2}[4] = Almr[2][-2][4]
	    (negative indexing - that's right!).

	    Note also that the allocation only requires the
	    maximum L value, so you don't have to worry about
	    those pesky +1's and -1's.

	    Printing out the values of Almr can then be done
	    using the utility (defined in diffuse.h) "printAlmr(Lmax)"

    THIS ROUTINE IS NOW OBSOLETE!  Use almN_alloc() instead!

***********************************************************************/

complex ***almr_alloc(int Lmax, int npts)
{
  complex ***almr,***zzz;
  int i,j;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;

   zzz = (complex ***)calloc(npts,sizeof(complex **));
  almr = (complex ***)calloc(npts,sizeof(complex **));

  for (j=0;j<npts;j++) {
    zzz[j] = (complex **)calloc(nl,sizeof(complex *));
    almr[j] = (complex **)calloc(nl,sizeof(complex *));
    for (i=0;i<nl;i++) {
      zzz[j][i] = (complex *)calloc(nm,sizeof(complex));
      almr[j][i] = (complex *)calloc(nm,sizeof(complex));
      almr[j][i] = zzz[j][i]+nl-1;
    }
  }

  free(zzz);
  return( almr );
}

/**<rlmr_alloc.c>*******************************************************

     Title: RLMR ALLOCation

      Call: int rlmr_alloc( int Lmax, int npts )

   Purpose: Allocate space for Rlmr

     Input: int Lmax = maximum L value

    Output: Rlmr[0,...,nr][0,..,nl][-nl,..,nl]

     Notes: This routine allocated rlmr so that indexing into
            Rlmr can be done in the most obvious manner, ie
	    Y_{l}{m}[k] = Rlmr[l][m][k], eg, Y_{2}{-2}[4] = Rlmr[2][-2][4]
	    (negative indexing - that's right!).

	    Note also that the allocation only requires the
	    maximum L value, so you don't have to worry about
	    those pesky +1's and -1's.

	    Printing out the values of Rlmr can then be done
	    using the utility (defined in diffuse.h) "printRlmr(Lmax)"

***********************************************************************/

float ***rlmr_alloc(int Lmax, int npts)
{
  float ***rlmr,***zzzz;
  int i,j;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;

  zzzz = (float ***)calloc(npts,sizeof(float **));
  rlmr = (float ***)calloc(npts,sizeof(float **));

  for (j=0;j<npts;j++) {
    zzzz[j] = (float **)calloc(nl,sizeof(float *));
    rlmr[j] = (float **)calloc(nl,sizeof(float *));
    for (i=0;i<nl;i++) {
      zzzz[j][i] = (float *)calloc(nm,sizeof(float));
      rlmr[j][i] = (float *)calloc(nm,sizeof(float));
      rlmr[j][i] = zzzz[j][i]+nl-1;
    }
  }

  free(zzzz);
  return( rlmr );
}

/**<alm_alloc.c>*******************************************************

     Title: ALM ALLOCation

      Call: complex **alm_alloc( int Lmax )

   Purpose: Allocate space for Alm

     Input: int Lmax = maximum L value

    Output: Alm[0,..,nl][-nl,..,nl]

     Notes: This routine allocated alm so that indexing into
            Alm can be done in the most obvious mannter, ie
	    A_{l}{m} = Alm[l][m], eg, A_{2}{-2} = Alm[2][-2]
	    (negative indexing - that's right!).

	    Note also that the allocation only requires the
	    maximum L value, so you don't have to worry about
	    those pesky +1's and -1's.

	    Printing out the values of Alm can then be done
	    using the utility (defined in diffuse.h) "printAlm(Lmax)"

	    The indexing loop for Alm is either:

	    for(i=0; i<nl; i++) {  
	      for(j=-i; j<=i; j++) {}
	    }

	    where nl = Lmax + 1, or
	    
	    for(i=0; i<=Lmax; i++) {  
	      for(j=-i; j<=i; j++) {}
	    }

***********************************************************************/

complex **alm_alloc(int Lmax)
{
  complex **alm,**zzz;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;
  int i;

  zzz = (complex **)calloc(nl,sizeof(complex *));
  alm = (complex **)calloc(nl,sizeof(complex *));

  for (i=0;i<nl;i++) {
    zzz[i] = (complex *)calloc(nm,sizeof(complex));
    alm[i] = zzz[i]+nl-1;
  }

  free(zzz);
  return( alm );
}

/**<almN_alloc.c>*******************************************************

     Title: ALMN ALLOCation

      Call: complex ***almN_alloc(int Lmax, int npts)

   Purpose: Allocate space for AlmN

     Input: int Lmax = maximum L value

    Output: AlmN[0,...,N][0,..,nl][-nl,..,nl]

     Notes: This routine allocated almN so that indexing into
            AlmN can be done in the most obvious manner, ie
	    Y_{l}{m}[k] = AlmN[l][m][k], eg, Y_{2}{-2}[4] = AlmN[2][-2][4]
	    (negative indexing - that's right!).

	    Note also that the allocation only requires the
	    maximum L value, so you don't have to worry about
	    those pesky +1's and -1's.

	    Printing out the values of AlmN can then be done
	    using the utility (defined in diffuse.h) "printAlmN(Lmax)"

***********************************************************************/

complex ***almN_alloc(int Lmax, int npts)
{
  complex ***almN;
  int i;

  almN = (complex ***)calloc(npts,sizeof(complex **));

  for (i=0;i<npts;i++) almN[i] = alm_alloc(Lmax) ;

  return( almN );
}

/**<almNM_alloc.c>*******************************************************

     Title: ALMN ALLOCation

      Call: complex ****almNM_alloc(int Lmax, int npts, int mpts)

   Purpose: Allocate space for AlmN

     Input: int Lmax = maximum L value

    Output: AlmN[0,...,M][0,...,N][0,..,nl][-nl,..,nl]

     Notes: This routine allocated almNM so that indexing into
            AlmN can be done in the most obvious manner, ie
	    Y_{l}{m}[k] = AlmN[l][m][k], eg, Y_{2}{-2}[4] = AlmN[2][-2][4]
	    (negative indexing - that's right!).

	    Note also that the allocation only requires the
	    maximum L value, so you don't have to worry about
	    those pesky +1's and -1's.

	    Printing out the values of AlmN can then be done
	    using the utility (defined in diffuse.h) "printAlmN(Lmax)"

***********************************************************************/

complex ****almNM_alloc(int Lmax, int npts, int mpts)
{
  complex ****almNM;
  int i;

  almNM = (complex ****)calloc(mpts,sizeof(complex ***));

  for (i=0;i<mpts;i++) almNM[i] = almN_alloc(Lmax,npts) ;

  return( almNM );
}

/**<alm_free.c>*******************************************************

     Title: ALM FREE

      Call: complex **alm_free( complex **alm, int Lmax )

   Purpose: Allocate space for Alm

     Input: int Lmax = maximum L value

    Output: Alm[0,..,nl][-nl,..,nl]

***********************************************************************/

void alm_free(complex **alm, int Lmax)
{
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;
  int i;

  for (i=0;i<nl;i++) free(alm[i]) ;

  free((complex *)alm);
  alm=NULL ;
}

/**<rlm_alloc.c>*******************************************************

     Title: RLM ALLOCation

      Call: float *rlm_alloc( int Lmax )

   Purpose: Allocate space for Rlm

     Input: int Lmax = maximum L value

    Output: Rlm[0,..,nl][-nl,..,nl]

     Notes: This routine allocated rlm so that indexing into
            Rlm can be done in the most obvious mannter, ie
	    R_{l}{m} = Rlm[l][m], eg, A_{2}{-2} = Rlm[2][-2]
	    (negative indexing - that's right!).

	    Note also that the allocation only requires the
	    maximum L value, so you don't have to worry about
	    those pesky +1's and -1's.

	    Printing out the values of Rlm can then be done
	    using the utility (defined in diffuse.h) "printRlm(Lmax)"

	    The indexing loop for Rlm is either:

	    for(i=0; i<nl; i++) {  
	      for(j=-i; j<=i; j++) {}
	    }

	    where nl = Lmax + 1, or
	    
	    for(i=0; i<=Lmax; i++) {  
	      for(j=-i; j<=i; j++) {}
	    }

***********************************************************************/

float **rlm_alloc(int Lmax)
{
  float **rlm,**zzz;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;
  int i;

  zzz = (float **)calloc(nl,sizeof(float *));
  rlm = (float **)calloc(nl,sizeof(float *));

  for (i=0;i<nl;i++) {
    zzz[i] = (float *)calloc(nm,sizeof(float));
    rlm[i] = zzz[i]+nl-1;
  }

  free(zzz);
  return( rlm );
}

/**<wlmr_alloc.c>*******************************************************

     Title: WLMR ALLOCation

      Call: complex ****wlmr_alloc(int Lmax, int nrad, int nvert)

   Purpose: Allocate space for Wlmr

     Input: int Lmax = maximum L value

    Output: Wlmr[0,..,nrad][0,..,nl][-nl,..,nl][0,..,nvert]

     Notes: This routine allocated wlmr so that indexing into
            Wlmr can be done in the most obvious manner, ie
	    Y_{l}{m}[k] = Wlmr[l][m][k], eg, Y_{2}{-2}[4] = Wlmr[2][-2][4]
	    (negative indexing - that's right!).

	    Note also that the allocation only requires the
	    maximum L value, so you don't have to worry about
	    those pesky +1's and -1's.

	    Printing out the values of Wlmr can then be done
	    using the utility (defined in diffuse.h) "printWlmr(Lmax)"

***********************************************************************/

complex ****wlmr_alloc(int Lmax, int nrad, int nvert)
{
  complex ****wlmr,****zzzz;
  int i,j,k;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;

  zzzz = (complex ****)calloc(nrad,sizeof(complex ***));
  wlmr = (complex ****)calloc(nrad,sizeof(complex ***));

  for (k=0;k<nrad;k++) {
    zzzz[k] = (complex ***)calloc(nl,sizeof(complex **));
    wlmr[k] = (complex ***)calloc(nl,sizeof(complex **));
    for (j=0;j<nl;j++) {
      zzzz[k][j] = (complex **)calloc(nm,sizeof(complex *));
      wlmr[k][j] = (complex **)calloc(nm,sizeof(complex *));
      wlmr[k][j] = zzzz[k][j]+nl-1;
      for (i=0;i<nm;i++) {
	zzzz[k][j][i] = (complex *)calloc(nvert,sizeof(complex));
      }
    }
  }

  free(zzzz);
  return( wlmr );
}

/**<ylm_alloc.c>*******************************************************

     Title: YLM ALLOCation

      Call: complex ***ylm_alloc(int Lmax, int npts)

   Purpose: Allocate space for Ylm

     Input: int Lmax = maximum L value

    Output: Ylm[0,..,nl][-nl,..,nl][0,...,npts]

     Notes: This routine allocated ylm so that indexing into
            Ylm can be done in the most obvious manner, ie
	    Y_{l}{m}[k] = Ylm[l][m][k], eg, Y_{2}{-2}[4] = Ylm[2][-2][4]
	    (negative indexing - that's right!).

	    Note also that the allocation only requires the
	    maximum L value, so you don't have to worry about
	    those pesky +1's and -1's.

	    Printing out the values of Ylm can then be done
	    using the utility (defined in diffuse.h) "printYlm(Lmax)"

***********************************************************************/

complex ***ylm_alloc(int Lmax, int npts)
{
  complex ***ylm,***zzz;
  int i,j;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;

  zzz = (complex ***)calloc(nl,sizeof(complex **));
  ylm = (complex ***)calloc(nl,sizeof(complex **));

  for (j=0;j<nl;j++) {
    zzz[j] = (complex **)calloc(nm,sizeof(complex *));
    ylm[j] = (complex **)calloc(nm,sizeof(complex *));
    ylm[j] = zzz[j]+nl-1;

    for (i=0;i<nm;i++) zzz[j][i] = (complex *)calloc(npts,sizeof(complex));

  }

  free(zzz);
  return( ylm );
}

/**<ylm_free.c>*******************************************************

     Title: YLM FREE

      Call: complex **ylm_free( complex **ylm, int Lmax )

   Purpose: Allocate space for Ylm

     Input: int Lmax = maximum L value

    Output: Ylm[0,..,nl][-nl,..,nl]

***********************************************************************/

void ylm_free(complex ***ylm, int Lmax)
{
  int ll,mm;

  for( ll=0 ; ll<=Lmax ; ll++ ) {
    for( mm=-ll ; mm<=ll ; mm++ ) {
      free(ylm[ll][mm]) ;
    }
    //free(ylm[ll]) ;
  }

  free((complex *)ylm);
  ylm=NULL ;
}

/**<pickalm.c>*******************************************************

     Title: PICK out an element of Alm

   Purpose: Extract element (Lpick,Mpick) out of Alm

      Call: complex pickalm( complex **Alm, 
		             int Lpick, int Mpick, 
			     int Lmax)

     Input: complex **Alm[][] - SHT coefficients

    Output: returns complex coeff Alm[l][m] if successful
            or ZNULLVAL if something went wrong

     Notes: Since alm_alloc now makes indexing easy,
            this routine just checks integrity of requested
	    L and M values before returning value.

***********************************************************************/

complex pickalm( complex **Alm, int Lpick, int Mpick, int Lmax )
{
  complex almpick ;

  if(!CHECK_L(Lpick,Lmax)) NEAR_FATAL_ZERROR("Invalid Lpick requested");
  if(!CHECK_M(Mpick,Lmax)) NEAR_FATAL_ZERROR("Invalid Mpick requested");
  
  almpick = Alm[Lpick][Mpick] ;

  return almpick;
}

/**<zvtoAlm.c>*******************************************************

     Title: Complex Vector TO ALM

   Purpose: Load Alm stored in a standard complex vector
            into the Alm format

      Call: complex *zvtoAlm( complex **Alm, complex *zv, int Lmax )

     Input: complex *zv[nlm] - Alm coeff's in linear vector

    Output: complex **Alm[][] 

***********************************************************************/

complex **zvtoAlm( complex **Alm, complex *zv, int Lmax )
{
  int i,j;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;
  int nn = 0;

  for (i=0; i<nl; i++) {
    for (j=-i; j<=i; j++) Alm[i][j] = zv[nn++] ;
  }

  return Alm ;
}

/**<Almtozv.c>*******************************************************

     Title: ALM TO Complex Vector 

   Purpose: Load Alm stored into a standard complex vector

      Call: complex *Almtozv( complex *zv, complex **Alm, int Lmax )

     Input: complex **Alm[][] 

    Output: complex *zv[nlm] - Alm coeff's in linear vector

***********************************************************************/

complex *Almtozv( complex *zv, complex **Alm, int Lmax )
{
  int i,j;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;
  int nn = 0;

  for (i=0; i<nl; i++) {
    for (j=-i; j<=i; j++) zv[nn++] = Alm[i][j] ;
  }

  return zv ;
}

int fac2(int l) {
  int r = 1;
  while (l > 0) {
    r *= l;
    l -= 2;
  }
  return r;
}

void alm2olm(complex **Alm, complex **Olm, int Lmax) {
  int ll, mm;
  complex zero;
  zero.r = 0.0;
  zero.i = 0.0;
  for( ll=0 ; ll<=Lmax ; ll+=2 ) {
    for( mm=-ll ; mm<=ll ; mm++ ) {
      if (ll%2 != 0) {
	Olm[ll][mm] = zero;
      }
      else if ((ll/2)%2 == 0) {
	Olm[ll][mm] = zscal((-1.0*fac2(ll+1))/(1.0*(ll+1)*fac2(ll)), Alm[ll][mm]);
      }
      else {
	Olm[ll][mm] = zscal((1.0*fac2(ll+1))/(1.0*(ll+1)*fac2(ll)), Alm[ll][mm]);
      }
    }
  }
}

/**<alm2shape.c>*******************************************************

     Title: Alm TO SHAPE

   Purpose: Reconstruct shape from alm

      Call: float **alm2shape( float **shape, complex **Alm, 
		               int Lmin, int Lmax,
			       float *az, float *po, int naz, int npo )

     Input: float **shape[nazim][npolr]
            float *az[nazim*npolr]
            float *po[nazim*npolr]

    Output: **shape[naz][npo]

***********************************************************************/

float **alm2shape( float **shape, complex **Alm, 
		   int Lmin, int Lmax,
		   float *az, float *po, int naz, int npo )
{
  int     i,j,ll,mm;
  int     nl =   Lmax+1;
  int     nm = 2*Lmax+1;
  int     iaz,ipo,iomega ;
  int     nomega = naz * npo ;
  float   AZmin = 0.0,   POmin = 0.0;
  float   AZmax = TWOPI, POmax = M_PI;
  float   daz = (AZmax-AZmin)/(naz-1.);
  float   dpo = (POmax-POmin)/(npo-1.);
  complex zsum ;

  complex ***Ylm ;

  Ylm = ylm_alloc( Lmax, nomega ) ;

  /* Create angles */

  for (ipo=0; ipo<npo; ipo++) {
    for (iaz=0; iaz<naz; iaz++) {
      iomega = iaz + ipo*naz ;
      az[iomega] = AZmin + iaz*daz ;
      po[iomega] = POmin + ipo*dpo ;
    }
  }

  /* generate spherical harmonics over all angles */

  SphericalHarmonicsYlm( Ylm, Lmin, Lmax, po, az, nomega );
  
  /* scale each (l,m) component of Ylm by Alm 
     (Computation done here in Ylm - notice that they
     are no longer just Ylm, but scaled versions ... */

  for( ll=Lmin ; ll<=Lmax ; ll++ ) {
    for( mm=-ll ; mm<=ll ; mm++ ) {

      zvscal( Ylm[ll][mm], Alm[ll][mm], Ylm[ll][mm] , nomega ) ;

    }
  }

  /* At each Omega, sum components */

  for (ipo=0; ipo<npo; ipo++) {
    for (iaz=0; iaz<naz; iaz++) {

      iomega = iaz + ipo*naz ;

      ZCLR(zsum) ;

      for( ll=Lmin ; ll<=Lmax ; ll++ ) {
	for( mm=-ll ; mm<=ll ; mm++ ) {

	  zsum = zadd(zsum,Ylm[ll][mm][iomega]) ;
	
	}
      }
      shape[iaz][ipo] = ZABS(zsum) ;
    }
  }

  /* The following routines not written yet.
     ylm_free(Ylm,Lmax,nomega) */

  return(shape) ;
}
/**<alm2egy.c>*******************************************************

     Title: Alm TO EnerGY

   Purpose: Compute energy in the l'th channel of Alm for all L

      Call: complex *alm2egy( complex *egyL, complex **Alm, int Lmax )

     Input: complex **Alm[][] - SHT coefficients

    Output: returns float energy in L channels of A[l][m] if successful
            or NULLVAL if something went wrong

***********************************************************************/

float *alm2egy( float *egyL, complex **Alm, int Lmax )
{
  int i,j,nn;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;
  complex zsum, zrow[nm];
  int debug = 0;

  /* Do we need scaling for each L channel?? */

  if(debug) {
    printf("alm2egy:\n");
    printAlm(Alm,Lmax);
  }

  for (i=0; i<nl; i++) {
    nn = 2*i + 1;
    getMrowFromAlm(zrow,Alm,i);
    zsum = zvcor(zrow,zrow,nn); 
    egyL[i] = zabs(zsum) ;
    if(debug) printf("alm2egy: nn = %i, zsum = (%f,%f), egyL[%i] = %f\n",
		     nn,zsum.r,zsum.i,i,egyL[i]);
  }

  return egyL ;
}

/**<alm2gra.c>*******************************************************

     Title: Alm TO Generalized Relative Anisotropy

   Purpose: Compute mean squared anisotropy from Alm

      Call: float alm2gra( complex **Alm, int Lmax )

     Input: complex **Alm[][] = SHT coefficients
            int Lmax = maximum L value computed in SHD

    Output: returns float gra

 Algorithm: This is sqrt[<Da/Di>^2] 
            where Da is the anisotropic components 
	    and   Di is the isotropic component

***********************************************************************/

float alm2gra( complex **Alm, int Lmax )
{
  int debug = 0;

  int l,nl =   Lmax+1;
  int m,nm = 2*Lmax+1;

  float Ctmp,Cl[nl] ;

  float scal, aniso ;

  if(debug) {
    printf("alm2gra:\n");
    printAlm(Alm,Lmax);
  }

  /* rescale coeff's by 1/(mean diffusion) */

  scal = 1./(Alm[0][0].r*4.*M_PI) ; /* This is 1/(a00*Y00) 
				     = 1/(mean diffusion) */
  if(debug) printf("alm2gra: scal = %f\n",scal);
  scaleAlm(Alm,scal,Lmax);

  /* Multipole coefficients 
     Compute all, even though we only use some below. */

  for (l=0; l<nl; l++) {
    Ctmp = 0.0 ;
    for(m=-l; m<=l; m++) Ctmp += ZSQR(Alm[l][m]) ;
    Cl[l] = Ctmp/(2.*l+1) ;
  }

  /* Anisotropy 
     Note: l=2,4,...,Lmax */

  aniso = 0.0 ;
  for (l=2; l<nl; l+=2) aniso += Cl[l] ;

  aniso = sqrt(aniso) ;

  return aniso ;
}

/**<alm2gfa_ozarslan.c>*******************************************************

     Title: Alm TO Generalized Fractional Anisotropy

   Purpose: Compute mean squared anisotropy from Alm

      Call: float alm2gfa( complex **Alm, int Lmax )

     Input: complex **Alm[][] = SHT coefficients
            int Lmax = maximum L value computed in SHD

    Output: returns float gfa

 Algorithm: This is sqrt[<Da/D>^2] 
            where Da is the anisotropic components 
	    and   D  is the full tensor

***********************************************************************/

float alm2gfa_ozarslan( complex **Alm, int Lmax )
{
  int debug = 0;

  int l,nl =   Lmax+1;
  int m,nm = 2*Lmax+1;

  float Ctmp,Cl[nl] ;

  float magni, aniso ;
  float eV;

#define gfa_scal sqrt(4.*M_PI) 

  if(debug) {
    printf("alm2gfa_ozarslan:\n");
    printAlm(Alm,Lmax);
  }

  /* Multipole coefficients 
     Compute all, even though we only use some below. */

  /*  for (l=0; l<nl; l++) { */
  /* Multipole coefficients 
     Compute just even */
  for (l=0; l<nl; l+=2) {
    Ctmp = 0.0 ;
    for(m=-l; m<=l; m++) Ctmp += ZSQR(Alm[l][m]) ;
    /*    Cl[l] = Ctmp/(2.*l+1) ; */
    /* shouldn't need scaling */
    Cl[l] = Ctmp;
  }

  /* Tensor magnitude */

  magni = 0.0 ;
  magni = Cl[0];
  /*  for (l=0; l<nl; l+=2) magni += Cl[l] ; */

  /* Anisotropy 
     Note: l=2,4,...,Lmax */

  aniso = 0.0 ;
  for (l=2; l<nl; l+=2) aniso += Cl[l] ;

  /* this gives V, where V is as defined in eq. 20 of 
     Ozarslan et al. 2006. */

  aniso = (aniso/magni)/9.0 ;

  /* eq. 20 of Ozarslan et al. 2005 */
  eV = 1.0 + (1.0 / (1.0 + 5000.0*aniso));

  /* eq. 19 of Ozarslan et al. 2005 */
  aniso = 1.0 - (1.0 / (1.0 + pow((250.0*aniso), eV)));

  return aniso ;
}

/**<alm2gfa.c>*******************************************************

     Title: Alm TO Generalized Fractional Anisotropy

   Purpose: Compute mean squared anisotropy from Alm

      Call: float alm2gfa( complex **Alm, int Lmax )

     Input: complex **Alm[][] = SHT coefficients
            int Lmax = maximum L value computed in SHD

    Output: returns float gfa

 Algorithm: This is sqrt[<Da/D>^2] 
            where Da is the anisotropic components 
	    and   D  is the full tensor

***********************************************************************/

float alm2gfa( complex **Alm, int Lmax )
{
  int debug = 0;

  int l,nl =   Lmax+1;
  int m,nm = 2*Lmax+1;

  float Ctmp,Cl[nl] ;

  float magni, aniso ;

#define gfa_scal sqrt(4.*M_PI) 

  if(debug) {
    printf("alm2gfa:\n");
    printAlm(Alm,Lmax);
  }

  /* Multipole coefficients 
     Compute all, even though we only use some below. */

  for (l=0; l<nl; l++) {
    Ctmp = 0.0 ;
    for(m=-l; m<=l; m++) Ctmp += ZSQR(Alm[l][m]) ;
    Cl[l] = Ctmp/(2.*l+1) ;
  }

  /* Tensor magnitude */

  magni = 0.0 ;
  for (l=0; l<nl; l+=2) magni += Cl[l] ;

  /* Anisotropy 
     Note: l=2,4,...,Lmax */

  aniso = 0.0 ;
  for (l=2; l<nl; l+=2) aniso += Cl[l] ;

  aniso = sqrt(aniso/magni) ;

  if(isnan(aniso)) aniso = 0.0 ; /* set nan's to zero */

  return aniso ;
}

/**<sumalm.c>*******************************************************

     Title: SUM Alm

   Purpose: Sum Alm over all M values for a given L

      Call: complex *sumalm( complex *sumL, complex **Alm, int Lmax )

     Input: complex **Alm[][] - SHT coefficients

    Output: returns complex sum over l of A[l][m] if successful
            or ZNULLVAL if something went wrong

***********************************************************************/

complex *sumalm( complex *sumL, complex **Alm, int Lmax )
{
  int i,j;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;
  complex zsum ;

  for (i=0; i<nl; i++) {
    ZCLR(zsum) ;
    for (j=-i; j<=i; j++) zsum = zadd(zsum,Alm[i][j]);
    sumL[i] = zsum ;
  }

  return sumL ;
}

/**<alm2tensor.c>*******************************************************

     Title: Alm TO TENSOR conversion

      Call: float **alm2tensor( float **Dap, complex **Alm, int Lmax )

   Purpose: Transformation between the Spherical and the 
            Cartesian tensor representations

     Input: Alm = SHD coefficients
             nl = # L values
             nm = # M values
	   Lmax = maximum L

    Output: Dap = Cartesian (apparent) diffusion tensor

***********************************************************************/

float **alm2tensor( float **Dap, complex **Alm, int Lmax )
{
  int i ;
  int debug = 1 ;

#define sqrt2 sqrt(2.)
#define sqrt3 sqrt(3.)
#define sqrt6 sqrt(6.)

#define isqrt2 (1./sqrt2)
#define isqrt3 (1./sqrt3)
#define isqrt6 (1./sqrt6)

#define a0scal sqrt(3.)  /* guessed at this one! */
#define a1scal sqrt(4.*M_PI/3.) 
#define a2scal sqrt(3./(10.*M_PI)) 

  complex a00 ;
  complex a10,a1p1,a1n1 ;
  complex a20,a2p1,a2p2,a2n1,a2n2 ;

  complex Axx,Axy,Axz ;
  complex Ayx,Ayy,Ayz ;
  complex Azx,Azy,Azz ;

  float   Dxx,Dxy,Dxz ;
  float   Dyx,Dyy,Dyz ;
  float   Dzx,Dzy,Dzz ;

  a00 = zscal( a0scal, pickalm(Alm,0,0,Lmax) ) ;

  a1n1 = zscal( a1scal, pickalm(Alm,1,-1,Lmax) ) ;
  a10  = zscal( a1scal, pickalm(Alm,1, 0,Lmax) ) ;
  a1p1 = zscal( a1scal, pickalm(Alm,1, 1,Lmax) ) ;

  a2n2 = zscal( a2scal, pickalm(Alm,2,-2,Lmax) ) ;
  a2n1 = zscal( a2scal, pickalm(Alm,2,-1,Lmax) ) ;
  a20  = zscal( a2scal, pickalm(Alm,2, 0,Lmax) ) ;
  a2p1 = zscal( a2scal, pickalm(Alm,2, 1,Lmax) ) ;
  a2p2 = zscal( a2scal, pickalm(Alm,2, 2,Lmax) ) ;

  /*
  Axx = (a2p2 + a2n2)/2 - a20/sqrt6 - a00/sqrt3;
  Axy = i*a10/sqrt2 - i*(a2p2 - a2n2)/2 ;
  Axz = ( (a1p1 + a1n1) - (a2p1 - a2n1) )/2 ;
  Ayx = -i*a10/sqrt2 - i*(a2p2-a2n2)/2 ;
  Ayy = -(a2p2 + a2n2)/2 - a20/sqrt6 - a00/sqrt3 ;
  Ayz = i*(-(a1p1-a1n1) + (a2p1+a2n1) )/2;
  Azx =  -( (a1p1+a1n1) + (a2p1-a2n1) )/2;
  Azy = i*( (a1p1-a1n1) + (a2p1+a2n1) )/2;
  Azz = (sqrt2*a20 - a00)/sqrt3;
  */

  { /* Axx */
    complex z1,z2,z3 ;
     z1 = zscal(.5,zadd(a2p2,a2n2)) ;
     z2 = zscal(-isqrt6,a20) ;
     z3 = zscal(-isqrt3,a00) ;
    Axx = zadd(zadd(z1,z2),z3);
  }

  { /* Axy */
    complex z1,z2,z3 ;
     z1 = zmuli(zscal(isqrt2,a10));
     z2 = zsub(a2p2,a2n2) ;
     z3 = zmuli(zscal(-.5,z2));
    Axy = zadd(z3,z1);
  }

  { /* Axz */
    complex z1,z2,z3 ;
     z1 = zadd(a1p1,a1n1);
     z2 = zsub(a2p1,a2n1);
     z3 = zsub(z1,z2) ;
    Axz = zscal(.5,z3);
  }

  { /* Ayx */
    complex z1,z2,z3 ;
    z1 = zmuli(zscal(-isqrt2,a10)); 
    z2 = zsub(a2p2,a2n2);
    z3 = zmuli(zscal(-.5,z2)); 
    Ayx = zadd(z2,z3);
  }

  { /* Ayy */
    complex z1,z2,z3 ;
     z1 = zscal(-.5,zadd(a2p2,a2n2)) ;
     z2 = zscal(-isqrt6,a20) ;
     z3 = zscal(-isqrt3,a00) ;
    Ayy = zadd(zadd(z1,z2),z3);
  }

  { /* Ayz */
    complex z1,z2,z3 ;
    z1 = zsub(a1n1,a1p1);
    z2 = zadd(a2p1,a2n1);
    z3 = zscal(.5,zadd(z1,z2)) ;
    Ayz = zmuli(z3);
  }

  { /* Azx */
    complex z1,z2,z3 ;
     z1 = zadd(a1p1,a1n1);
     z2 = zsub(a2p1,a2n1);
     z3 = zadd(z1,z2) ;
    Azx = zscal(-.5,z3);
  }

  { /* Azy */
    complex z1,z2,z3 ;
    z1 = zsub(a1p1,a1n1);
    z2 = zadd(a2p1,a2n1);
    z3 = zscal(.5,zadd(z1,z2));
    Azy = zmuli(z3);
  }

  { /* Azz */
    complex z1,z2,z3 ;
     z1 = zscal(sqrt2,a20);
     z2 = zscal(-1.,a00);
     z3 = zadd(z1,z2);
    Azz = zscal(isqrt3,z3);
  }

  /* Compute matrix components from A's */

  Dxx = Axx.r ;  Dxy = Axy.r ;  Dxz = Axz.r ;
  Dyx = Ayx.r ;  Dyy = Ayy.r ;  Dyz = Ayz.r ;
  Dzx = Azx.r ;  Dzy = Azy.r ;  Dzz = Azz.r ;

  /* Fill output matrix */

  Dap[0][0] = Dxx;  Dap[0][1] = Dxy;  Dap[0][2] = Dxz;
  Dap[1][0] = Dyx;  Dap[1][1] = Dyy;  Dap[1][2] = Dyz;
  Dap[2][0] = Dzx;  Dap[2][1] = Dzy;  Dap[2][2] = Dzz;

  return(Dap) ;
}

/*--------------------- Fiber routines. 10-Aug-03 ----------------------*/

/**<almcorr.c>*******************************************************

     Title: CORRelate 2 Alm

   Purpose: Correlate two Alm's

      Call: complex almcorr( complex **Alm1, complex **Alm2, int Lmax )

     Input: complex **Alm1[][] = SHT coefficients
            complex **Alm2[][] = SHT coefficients
	    int Lmax = maximum order

    Output: complex corr(Alm1,Alm2)

***********************************************************************/

complex almcorr( complex **Alm1, complex **Alm2, int Lmax )
{
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;
  complex zc;
  int debug = 0;

  int nlm = numAlm(Lmax) ;
  complex *vAlm1 = zv_alloc(nlm) ;
  complex *vAlm2 = zv_alloc(nlm) ;

  if(debug) {
    printf("almcorr:\n");
    printAlm(Alm1,Lmax);
    printAlm(Alm2,Lmax);
  }

  Almtozv(vAlm1,Alm1,Lmax) ;
  Almtozv(vAlm2,Alm2,Lmax) ;

  zc = zvcor(vAlm1,vAlm2,nlm) ;

  free(vAlm1) ;
  free(vAlm2) ;
  return zc ;
}

/**<loadalm.c>*******************************************************

     Title: SUM Alm

   Purpose: load AlmFrom into AlmTo

      Call: *loadalm( complex **AlmTo, complex **AlmFrom, int Lmax )

     Input: complex **AlmFrom[][] - SHT coefficients

    Output: complex **AlmTo[][] - SHT coefficients

***********************************************************************/

complex **loadalm( complex **AlmTo, complex **AlmFrom, int Lmax )
{
  int i,j;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;

  for (i=0; i<nl; i++) {
    for (j=-i; j<=i; j++) AlmTo[i][j] = AlmFrom[i][j];
  }

  return AlmTo ;
}

/**<almr_realloc.c>*******************************************************

     Title: ALMR REALLOCation

      Call: int almr_realloc( complex **almr, int Lmax, int npts )

   Purpose: Allocate space for Almr

     Input: int Lmax = maximum L value

    Output: Almr[0,...,nr][0,..,nl][-nl,..,nl]

     Notes: This routine allocated almr so that indexing into
            Almr can be done in the most obvious manner, ie
	    Y_{l}{m}[k] = Almr[l][m][k], eg, Y_{2}{-2}[4] = Almr[2][-2][4]
	    (negative indexing - that's right!).

	    Note also that the allocation only requires the
	    maximum L value, so you don't have to worry about
	    those pesky +1's and -1's.

	    Printing out the values of Almr can then be done
	    using the utility (defined in diffuse.h) "printAlmr(Lmax)"

***********************************************************************/

complex ***almr_realloc(complex ***almr, int Lmax, int npts)
{
  complex ***zzzz;
  int i,j;
  int nl =   Lmax+1;
  int nm = 2*Lmax+1;
  int debug = 1;
  int nn = 0;

  zzzz = (complex ***)calloc(npts,sizeof(complex **));
  almr = (complex ***)realloc((char *)almr,sizeof(complex **)*npts);

  for (j=0;j<npts;j++) {
    zzzz[j] = (complex **)calloc(nl,sizeof(complex *));
    almr[j] = (complex **)realloc((char *)almr[j],sizeof(complex *)*nl);
    for (i=0;i<nl;i++) {
      zzzz[j][i] = (complex *)calloc(nm,sizeof(complex));
      almr[j][i] = (complex *)realloc((char *)almr[j][i],sizeof(complex)*nm);
      almr[j][i] = zzzz[j][i]+nl-1;
    }
  }

  free(zzzz);
  return( almr );
}
