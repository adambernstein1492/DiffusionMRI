% /*********************************************************************************
% Wrapper function for bloch2.c that returns output in a multidimensional array
% 
% USAGE: m = blok(rho,phi,g,rast,Min,z,t1,t2,v,off,tkz)
%     Input:
%         rho      row vector of B1 envelop (G)
%         phi      row vector of phi (radians)
%         g        row vector of grad envelope (G/cm)
%         rast     raster (us)
%         Min      1x3 vector of initial magnetization (M/Mo)
%         z        row vector of Z positions (cm)
%         t1       (ms) enter 0 for no relaxation
%         t2       (ms) enter 0 for no relaxation
%         v        row vector of velocities (cm/s)
%         off      (Hz)  
%         tkz      time at Kz=0 (fraction of duration)
%         
%     Output:
%         m       3 by nt by nz by nv matrix of output magnetization
% 
% ************************
% This is a Matlab Mex file adapted from ecwong's bloch.c code.
% Runtime is approx 200x faster than a native bloch.m script.
% 040414 Matt Cronin
% 050513 ECW fixed bug in velocity input
% 050602 ECW created bloch2 with restructured inputs
% *********************************************************************************/

function m = blok(rho,phi,g,rast,Min,z,t1,t2,v,off,tkz);

nv = length(v);
nz = length(z);
nt = length(rho);

m = bloch2(rho,phi,g,rast,Min,z,t1,t2,v,off,tkz);
m = reshape(m,3,nt,nz,nv);
