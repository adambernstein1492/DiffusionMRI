#!/usr/bin/env python
import math
for q in range(15,91,15):
    for phi in range(0, 360, 15):
        parts = q*q*1000

    # if q >= 30:
    # else:
    #     parts = 2500000
        
        fname = "dpfg-q-%d-phi-%d.txt" % (q, phi)
        f = open(fname, 'w')
        G = 10000.0*math.pi*q/26752.0
        f.write("""# q = %f mm^(-1)
random seed = %d
Delta = 38000  # in microseconds (Delta - delta using standard definitions)
delta = 2000  # in microseconds
mixing time = 0 # in microseconds (only for DPFG pulses)
gradient strength = %f  # in G/cm
ramp = 0  # in microseconds (for both leading and trailing edges)
pulse = 2  # 0 = gradient echo; 1 = spin echo; 2 = DPFG
snr = 100  # signal-to-noise ratio (percent)
voxel size = 8.0  # in microns
volume size = 8.0  # in microns
cylinder radius = 4.0
cylinder height = 8.8
nstep = 1  # number of steps (not used)
step size = 25  # time step length in microseconds
ntess = 0  # tesselation number (not used for DPFG pulses)
nparts = %d  # number of diffusing molecules
DPFG angle = %f
signal file = signal-dpfg-q-%d-phi-%d.dat  # output file for signal data
complex signal file = complex-dpfg-q-%d-phi-%d.dat  # output file for signal data
structure file = z-cylinder-radius-1.mdl  # mcell file describing geometry
""" % (q, q*100, G, parts, phi, q, phi, q, phi))
    f.close()

            

