#!/usr/bin/env python

# Code for testing results for diffusion within a single cylinder of
# unbounded extent.

import subprocess
import sys
import os
import os.path
import math

# Get the cylinder building utility.
sys.path.append('../utilities')
import cylinder

executable = sys.argv[1]
full_executable = os.path.abspath(executable)
if (not os.path.exists(full_executable)) or (not os.access(full_executable, os.X_OK)):
    print "The given program does not exist or is not executable."
    sys.exit(1)

structure_file = open('cylinder.mdl','w')
cylinder.write(structure_file, 64, 0.5)
structure_file.close()

g_squared = 26752.0 * 26752.0   # gamma in rad/G/s

standard_values = '''Delta = %d
delta = %d
mixing time = 0  # in microseconds (only for DPFG pulses)
gradient strength = %f  # in G/cm
ramp = 0  # in microseconds (for both leading and trailing edges)
pulse = %d  # 0 = gradient echo; 1 = spin echo; 2 = DPFG
periodic = %d  # NOT periodic boundary conditions
snr = 100  # signal-to-noise ratio (percent)
voxel size = %f  # in microns
cylinder array = 0  # periodic array
three region = 0  # core-sheath-bath like in Sen-Basser paper
permeable = 0  # default is permeable = 1
gradient directions = 2

core radius = %f  # in microns (really the radius)
cylinder height = %f  # in microns
x dim = %f  # in microns
y dim = %f  # in microns
z dim = %f  # in microns

nstep = 1  # number of steps (not used)
step size = 10  # time step length in microseconds
nparts = 1000000  # number of diffusing molecules
structure file = ./cylinder.mdl  # mcell file describing geometry
signal file = %s_signal.dat # output file for signal amplitude
complex signal file = %s_complex_signal.dat # output file for complex signal data
'''

input_filename = 'cylinder_test.txt'

for b in [600, 1800]:
    for Delta_eff in [1800, 3600]:
        for delta_frac in [1, 15, 29]:
            delta = delta_frac*Delta_eff/30
            Delta = Delta_eff + delta/3
            Delta_ds = Delta - delta
            gradient_strength = math.sqrt(1.0e20 * b/(g_squared 
                                          * delta * delta * Delta_eff))
            for pulse in [0, 1]:
                for diam in [4.0, 8.0]:
                    periodic = 0
                    vx_size = 10.0
                    vl_size = 25.0
                    variant_name = 'cyl-b%d-d%d-D%d-%d-%d' % (b, delta, Delta_eff,
                                                          pulse, periodic)
                    if os.path.exists(variant_name):
                        print 'Removing previous %s ...' % (variant_name)
                        if os.path.isdir(variant_name):
                            shutils.rmtree(variant_name)
                        else:
                            os.remove(variant_name)
                    os.mkdir(variant_name)
                    f = open(os.path.join(variant_name,input_filename),'w')
                    f.write(standard_values % (Delta_ds, delta,
                                               gradient_strength,
                                               pulse,
                                               periodic,
                                               vx_size, diam,
                                               2.0*vl_size, vl_size,
                                               vl_size, vl_size,
                                               variant_name, variant_name))
                    f.close()
                    if (subprocess.call([executable, input_filename])):
                        print 'An error occured running %s.' % (executable)
                        sys.exit(1)

                    periodic = 1
                    vx_size = 15.0
                    vl_size = 10.0
                    variant_name = 'cyl-b%d-d%d-D%d-%d-%d' % (b, delta, Delta_eff,
                                                          pulse, periodic)
                    f = open(input_filename,'w')
                    f.write(standard_values % (Delta_ds, delta,
                                               gradient_strength,
                                               pulse,
                                               periodic,
                                               vx_size, diam,
                                               2.0*vl_size, vl_size,
                                               vl_size, vl_size,
                                               variant_name, variant_name))
                    f.close()
                    if (subprocess.call([executable, input_filename])):
                        print 'An error occured running %s.' % (executable)
                        sys.exit(1)


