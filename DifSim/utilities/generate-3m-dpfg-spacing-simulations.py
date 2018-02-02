#!/usr/bin/env python
import math
import os
q = 45
s = 2
for t in [0.1, 0.5, 0.9]:
    for phi in range(0, 360, 15):
        parts = q*q*100

    # if q >= 30:
    # else:
    #     parts = 2500000
        
        dirname = "dpfg-s-%d-t-%0.1f-q-%d-phi-%d" % (s, t, q, phi)
        os.mkdir(dirname)
        fname = "%s/dpfg-s-%d-t-%0.1f-q-%d-phi-%d.txt" % (dirname, s, t, q, phi)
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
x dim = 10.0 # 8 for spacing == 0, but we have spacing = 2.0
y dim = 17.3205080757  # sqrt(3) * 10
z dim = 15
cylinder array = 1
three region = 1  # core-sheath-bath like in Sen-Basser paper
core diffusion = 2.0e-5
sheath concentration = 0.15
core radius = 4.0
sheath radius = %0.1f
cylinder height = 20
nstep = 1  # number of steps (not used)
step size = 25  # time step length in microseconds
ntess = 0  # tesselation number (not used for DPFG pulses)
nparts = %d  # number of diffusing molecules
DPFG angle = %f
signal file = signal-dpfg-s-%d-t-%0.1f-q-%d-phi-%d.dat  # output file for signal data
complex signal file = complex-dpfg-s-%d-t-%0.1f-q-%d-phi-%d.dat  # output file for signal data
structure file = ../z-cylinder-radius-1-c.mdl  # mcell file describing geometry
sheath file = ../z-cylinder-radius-1-s.mdl  # mcell file describing geometry
""" % (q, phi*100+1, G, 4.0+t, parts, phi, s, t, q, phi, s, t, q, phi))
        f.close()

            

