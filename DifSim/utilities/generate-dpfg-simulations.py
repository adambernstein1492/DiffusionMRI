#!/usr/bin/env python
import math
for q in range(0,62,2):
    if q >= 30:
        parts = 10000000
    else:
        parts = 2500000
    fname = "dpfg-dt-25-tm-100-q-%d.txt" % (q)
    f = open(fname, 'w')
    G = 10000.0*math.pi*q/26752.0
    f.write("""# q = %f mm^(-1)
random seed = %d
Delta = 148000  # in microseconds (Delta - delta using standard definitions)
delta = 2000  # in microseconds
mixing time = 100000 # in microseconds (only for DPFG pulses)
gradient strength = %f  # in G/cm
ramp = 0  # in microseconds (for both leading and trailing edges)
pulse = 2  # 0 = gradient echo; 1 = spin echo; 2 = DPFG
snr = 100  # signal-to-noise ratio (percent)
voxel size = 30  # in microns
volume size = 30  # in microns
volume height = 28  # in microns
cylinder diameter = 28.639
cylinder radius = 14.32025
cylinder height = 30.0
nstep = 1  # number of steps (not used)
step size = 25  # time step length in microseconds
ntess = 0  # tesselation number (not used for DPFG pulses)
nparts = %d  # number of diffusing molecules
signal file = signal-dpfg-dt-25-tm-100-q-%d.dat  # output file for signal data
complex signal file = complex-dpfg-dt-25-tm-100-q--%d.dat  # output file for signal data
structure file = new-cylinder.mdl  # mcell file describing geometry
""" % (q, q*100, G, parts, q, q))
    f.close()

            

