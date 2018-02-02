#!/usr/bin/env python

import sys
import math
import numpy
import optparse

# Define gamma squared, since it's used most often in that form.
gamma_squared = 26752.0 * 26752.0 # (gamma is in rad/G/s)

# allowable input combinations:
#
# d D b -> g
# d e b -> g
# d dg b -> g
# D e b -> g
# D dg b -> g
# e dg b -> g
# d D g -> b
# d e g -> b
# d dg g -> b
# D e g -> b
# D dg g -> b
# e dg g -> b
# d b g -> dg
# D b g -> dg (& d)
# e b g -> dg (& d)
# dg b g -> d

# Set up option parser.
parser = optparse.OptionParser(description='Generate input files for DiffSim', version='1.0')
parser.add_option('-g', '--gradient', type=float, metavar='G', help='Gradient strength in G/cm')
parser.add_option('-d', '--delta', type=int, help='delta in microseconds')
parser.add_option('-D', '--Delta', type=int, help='Delta in microseconds')
parser.add_option('-e', '--Delta_eff', type=int, help='Delta - delta/3 (or tau) in microseconds')
parser.add_option('--dg', type=int, metavar='Delta_gap', help='Delta - delta (i.e. the time between pulses)')
parser.add_option('-b', '--bval', type=float, help='b value in s/mm^2')
parser.add_option('-t', '--timestep', type=int, help='timestep in microseconds')

# Parse command line options.
(options, args) = parser.parse_args()

delta = None
Delta = None
Delta_eff = None
dg = None
bval = None
gradient = None

# sort out timing parameters: need delta and dg for DiffSim
input_count = 0
if options.delta:
    input_count += 1
    delta = options.delta
if options.Delta:
    input_count += 1
    Delta = options.Delta
    if input_count == 2:
        dg = Delta - delta
if options.Delta_eff:
    if input_count == 2:
        print 'Pulse timing is overspecified: ignoring Delta_eff.'
    else:
        input_count += 1
        Delta_eff = options.Delta_eff
        if delta:
            dg = Delta_eff - 2*delta/3.0
        elif Delta:
            delta = 3*(Delta - Delta_eff)
            dg = Delta - delta
if options.dg:
    if input_count == 2:
        print 'Pulse timing is overspecified: ignoring dg.'
    else:
        input_count += 1
        dg = options.dg
        if Delta:
            delta = Delta - dg
        if Delta_eff:
            delta = 3*(Delta_eff - dg)/2.0

if input_count == 2:
    if options.gradient and options.bval:
        print 'Pulse sequence is overspecified: ignoring b'
    if options.gradient:
        input_count += 1
        gradient = options.gradient
        bval = gamma_squared*gradient*gradient*delta*delta \
            *(dg + 2.0*delta/3.0)/1.0e20
    elif options.bval:
        input_count += 1
        bval = options.bval
        gradient = math.sqrt(1.0e20 * bval
                             /(gamma_squared*delta*delta
                               *(dg + 2.0*delta/3.0)))
elif input_count == 1:
    if options.bval and options.gradient:
        bval = options.bval
        gradient = options.gradient
        input_count += 2
        if delta:
            dg = 1.0e20 * bval / (gamma_squared*gradient*gradient
                                  *delta*delta) - 2*delta/3.0
        elif dg:
            delta = numpy.roots([2.0/3.0, dg, 0.0,
                                       -1e20*bval
                                       /(gamma_squared*gradient*gradient)
                                       ])[2].real
        elif Delta:
            roots = numpy.roots([-1.0/3.0, Delta, 0.0,
                                        -1e20*bval
                                        /(gamma_squared
                                          *gradient*gradient)
                                        ])
            # take the first real root delta value greater than Delta
            for r in roots:
                if numpy.isreal(r) and r > Delta:
                    delta = r
                    break
            dg = Delta - delta
        elif Delta_eff:
            delta = math.sqrt(1.0e20*bval
                              /(Delta_eff*gamma_squared))/gradient
            dg = Delta_eff - 2*delta/3.0

if input_count < 3:
    print '''Pulse sequence is underspecified: please supply one time parameter
and both G AND b, or two time parameters and either G OR b.'''
    sys.exit()

print delta, dg, bval, gradient
sys.exit()



