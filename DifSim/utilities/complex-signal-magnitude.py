#!/usr/bin/env python

import sys
import cmath

if len(sys.argv) == 1:
    print "Usage: complex-signal-magnitude.py <complex-signal-file> [<signal-file>]"

complexSignalFile = sys.argv[1]
signalFile = ""
if len(sys.argv) == 3:
    signalFile = sys.argv[2]

cf = open(complexSignalFile)
cv = []
for l in cf:
    r, i = map(float, l.strip().split(" "))
    cv.append(complex(r, i))
    print abs(complex(r, i))
