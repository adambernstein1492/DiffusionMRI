/*********************************************************************************
Bloch Equation simulation of an RF pulse

USAGE: bloch2(rho,phi,g,rast,Min,z,t1,t2,v,off,tkz)
    Input:
        rho      row vector of B1 envelop (G)
        phi      row vector of phi (radians)
        g        row vector of grad envelope (G/cm)
        rast     raster (us)
        Min      1x3 vector of initial magnetization (M/Mo)
        z        row vector of Z positions (cm)
        t1       (ms) enter 0 for no relaxation
        t2       (ms) enter 0 for no relaxation
        v        row vector of velocities (cm/s)
        off      (Hz)  
        tkz      time at Kz=0 (fraction of duration)
        
    Output:
        Mout     (nt*nz*nv*3) matrix of output magnetization
                   looping in output is m[xyz] inside, then t, z, and v outside

************************
This is a Matlab Mex file adapted from ecwong's bloch.c code.
Runtime is approx 200x faster than a native bloch.m script.
040414 Matt Cronin
050513 ECW fixed bug in velocity input
050602 ECW created bloch2 with restructured inputs
*********************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/fcntl.h>
#include "mex.h"

#define GAM 26747.5    /* gamma in 1/(s-G) */

int bloch(double* rho, double* phi, double* g, int cgrad, int nrf, double rast, double* Min, double* z,
      int nz, double t1,double t2, double* v, int nv, double off, double tkz, double *Mout)

{ 
   register int i,j,k,w2,ind,count,size;
   double tx,ty,tz,et1,et2,mx,my,mz;
   double ca,sa,cp,sp,cp2,sp2,frot,cf,sf,a,p;
   double time,dt,*freq,*cflow,*sflow,tmp;
   double ztmp,*gg,*sr,*cr,gc;
   double *r11,*r12,*r13,*r21,*r22,*r23,*r31,*r32,*r33;        
        
/* get standard arguments */
    off = 2*M_PI*off;					/* convert from Hz to rads */
    t1 = 0.001 * t1;
    t2 = 0.001 * t2;
    dt = 0.000001 * rast;				/* us to seconds */
	time = nrf*dt;
	ind = 0;
            
/* Allocate memory */
    r11 = (double *) mxCalloc((unsigned) nrf, sizeof(double));
	r12 = (double *) mxCalloc((unsigned) nrf, sizeof(double));
	r13 = (double *) mxCalloc((unsigned) nrf, sizeof(double));
	r21 = (double *) mxCalloc((unsigned) nrf, sizeof(double));
	r22 = (double *) mxCalloc((unsigned) nrf, sizeof(double));
	r23 = (double *) mxCalloc((unsigned) nrf, sizeof(double));
	r31 = (double *) mxCalloc((unsigned) nrf, sizeof(double));
	r32 = (double *) mxCalloc((unsigned) nrf, sizeof(double));
	r33 = (double *) mxCalloc((unsigned) nrf, sizeof(double));
	gg = (double *) mxCalloc((unsigned) nrf, sizeof(double));
	sr = (double *) mxCalloc((unsigned) nz, sizeof(double));
	cr = (double *) mxCalloc((unsigned) nz, sizeof(double));


/* make rotation matrices */
	for (i=0;i<nrf;i++) {
		a = GAM*dt*rho[i];
		p = phi[i];
		ca = cos(a);
		sa = sin(a);
		cp = cos(p);
		sp = sin(p);
		cp2 = cp*cp;
		sp2 = sp*sp;
		r11[i] = cp2 + ca*sp2;
		r12[i] = (ca-1)*cp*sp;
		r13[i] = sa*sp;
		r21[i] = r12[i];
		r22[i] = ca*cp2 + sp2;
		r23[i] = cp*sa;
		r31[i] = -r13[i];
		r32[i] = -r23[i];
		r33[i] = ca;
	}

/* Calculate relaxation factors */
	if (t1*t1>1.e-10) et1 = exp(-dt/t1);
	else et1 = 1.;
	if (t2*t2>1.e-10) et2 = exp(-dt/t2);
	else et2 = 1.;
        
/* Precalculate flow and frequency table and refocussing phase twist if possible */
	if (cgrad) {

		/* convert grad to radians/s/cm */
		gc = GAM * g[0];
	
		/* Calculate flow factors */	
		sflow = (double *) mxCalloc((unsigned) nv, sizeof(double));
		cflow = (double *) mxCalloc((unsigned) nv, sizeof(double));
		for(i=0;i<nv;i++) {
			tmp = v[i]*gc*dt*dt;
			cflow[i] = cos(tmp);
			sflow[i] = sin(tmp);
		}
	
		/* Construct frequency table */
		freq = (double *) mxCalloc((unsigned) nz, sizeof(double));
		for (i=0;i<nz;i++) {
			freq[i] = z[i]*gc + off;
			tmp = -(z[i]*gc)*(1.-tkz)*time;
			sr[i] = sin(tmp);
			cr[i] = cos(tmp);
		}
	}
	else {
        for (i=0;i<nrf;i++) gg[i] = GAM * g[i];
		if (tkz<0.) tkz = 0.;
		else if (tkz>1.) tkz = 1.;
		for (tmp=0.,i=tkz*nrf+0.5;i<nrf;i++) tmp += gg[i]*dt;
		for (i=0;i<nz;i++) {
			ztmp = -z[i]*tmp;
			sr[i] = sin(ztmp);
			cr[i] = cos(ztmp);
		}
	}

/* GO !! */
	if (cgrad) {
		count=0;
		for(k=0;k<nv;k++) {
			for (i=0;i<nz;i++) {
				frot = freq[i]*dt;
				cf = cos(frot);
				sf = sin(frot);
				mx = Min[0];
				my = Min[1];
				mz = Min[2];
				for (j=0;j<nrf;j++) {
					tx = mx*r11[j] + my*r12[j] + mz*r13[j];
					ty = mx*r21[j] + my*r22[j] + mz*r23[j];
					tz = mx*r31[j] + my*r32[j] + mz*r33[j];
					mx = et2*(tx*cf+ty*sf);
					my = et2*(-tx*sf+ty*cf);
					mz = 1.-(1.-tz)*et1;
					tmp = cf*cflow[k] + sf*sflow[k];
					sf = -cf*sflow[k] + sf*cflow[k];
					cf = tmp;
					ind = 3*count;
					Mout[ind++] = mx;
					Mout[ind++] = my;
					Mout[ind] = mz;
					count++;
				}
				Mout[ind-2] = mx*cr[i]+my*sr[i];
				Mout[ind-1] = -mx*sr[i]+my*cr[i];
			}
		}
	} 
	else {
		count=0;
		for(k=0;k<nv;k++) { 
			for (i=0;i<nz;i++) {
				ztmp = z[i];
				mx = Min[0];
				my = Min[1];
				mz = Min[2];
				for (j=0;j<nrf;j++) {
					frot = (gg[j]*ztmp+off)*dt;
					cf = cos(frot);
					sf = sin(frot);
					tx = mx*r11[j] + my*r12[j] + mz*r13[j];
					ty = mx*r21[j] + my*r22[j] + mz*r23[j];
					tz = mx*r31[j] + my*r32[j] + mz*r33[j];
					mx = et2*(tx*cf+ty*sf);
					my = et2*(-tx*sf+ty*cf);
					mz = 1.-(1.-tz)*et1;
					ztmp += v[k]*dt;
					ind = 3*count;
					Mout[ind++] = mx;
					Mout[ind++] = my;
					Mout[ind] = mz;
					count++;
				}
				Mout[ind-2] = mx*cr[i]+my*sr[i];
				Mout[ind-1] = -mx*sr[i]+my*cr[i];
			}
		}
	}
	mxFree(r11); mxFree(r12); mxFree(r13);
	mxFree(r21); mxFree(r22); mxFree(r23);
	mxFree(r31); mxFree(r32); mxFree(r33);
	mxFree(gg); mxFree(sr); mxFree(cr);
	if (cgrad) {
		mxFree(sflow); mxFree(cflow); 
		mxFree(freq); 
	}
	
	return;
}

/****************************************************************************************
                        START MEX GATEWAY FUNCTION
****************************************************************************************/
                        
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
    int i,nrf,ng,cgrad,nz,nv,w2,size;
    char name[40] = "myvar";
    int *dims;
    double rast,t1,t2,off,tkz;
    double *rho;
    double *phi;
    double *g;
    double *Min;
    double *z;
    double *v;
    double *Mout;
    
    /* Check number of inputs */
    if (nrhs != 11) mexErrMsgTxt("Oops! Incorrect number of inputs to bloch.");
       
    /* Assign Pointers. */
    rho = (double *)mxGetPr(prhs[0]);
    phi = (double *)mxGetPr(prhs[1]);
    g = (double *)mxGetPr(prhs[2]);
    rast = mxGetScalar(prhs[3]);
    Min = (double *)mxGetPr(prhs[4]);
    z = (double *)mxGetPr(prhs[5]);
    t1 = mxGetScalar(prhs[6]);
    t2 = mxGetScalar(prhs[7]);
    v = (double *)mxGetPr(prhs[8]);
    off = mxGetScalar(prhs[9]);
    tkz = mxGetScalar(prhs[10]);
    
    dims = (int *)mxGetDimensions(prhs[0]);
    nrf = dims[0]*dims[1]; /* let it be a row or column vector */
    dims = (int *)mxGetDimensions(prhs[1]);
    if (dims[0]*dims[1] != nrf) mexErrMsgTxt("rho and phi not the same size\n");
    dims = (int *)mxGetDimensions(prhs[2]);
    ng = dims[0]*dims[1];
    if (ng == 1) cgrad = 1;
    else if (ng == nrf) cgrad = 0;
    else mexErrMsgTxt("grad is a funny size\n");
    
    dims = (int *)mxGetDimensions(prhs[5]);
    nz = dims[0]*dims[1];
    dims = (int *)mxGetDimensions(prhs[8]);
    nv = dims[0]*dims[1];
    
/*    mexPrintf("nrf=%i nz=%i nv=%i\n", nrf, nz, nv); */
            
    /* Alocate memory for output array. */
    plhs[0] = mxCreateDoubleMatrix(1,nrf*nv*nz*3, mxREAL);
       
    /* Assign output array */
    Mout = mxGetPr(plhs[0]);

    /* Go C routine */
    bloch(rho,phi,g,cgrad,nrf,rast,Min,z,nz,t1,t2,v,nv,off,tkz,Mout);
    return;
}

