#!/usr/bin/env python
import math
q = 45
for s in [0, 1, 2, 4, 8]:
    for phi in range(0, 360, 15):
        parts = q*q*100

    # if q >= 30:
    # else:
    #     parts = 2500000
        
        fname = "dpfg-s-%d-q-%d-phi-%d.txt" % (s, q, phi)
        f = open(fname, 'w')
        G = 10000.0*math.pi*q/26752.0
        f.write("""# q = %f mm^(-1)
periodic = 1
random seed = %d
Delta = 38000  # in microseconds (Delta - delta using standard definitions)
delta = 2000  # in microseconds
mixing time = 0 # in microseconds (only for DPFG pulses)
gradient strength = %f  # in G/cm
ramp = 0  # in microseconds (for both leading and trailing edges)
pulse = 2  # 0 = gradient echo; 1 = spin echo; 2 = DPFG
snr = 100  # signal-to-noise ratio (percent)
voxel size = 80.0  # in microns
x dim = %f # 8 for spacing == 0
y dim = %f # 13.85640646 for spacing = 0
z dim = 7
cylinder array = 1
cylinder diameter = 7.992
cylinder radius = 4.0
nstep = 1  # number of steps (not used)
step size = 25  # time step length in microseconds
ntess = 0  # tesselation number (not used for DPFG pulses)
nparts = %d  # number of diffusing molecules
DPFG angle = %f
signal file = signal-dpfg-s-%d-q-%d-phi-%d.dat  # output file for signal data
complex signal file = complex-dpfg-s-%d-q-%d-phi-%d.dat  # output file for signal data
structure file = cylinder-8.mdl  # mcell file describing geometry
""" % (q, phi*100+1, G, 8.0+s, math.sqrt(3)*(8.0+s), parts, phi, s, q, phi, s, q, phi))
    f.close()

            

