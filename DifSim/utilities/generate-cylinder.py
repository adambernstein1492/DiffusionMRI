#!/usr/bin/env python

import math

facets = 64
radius = 0.5
print "structure POLYGON_LIST {\n  VERTEX_LIST {"

# points at one end
for i in range(facets):
    x = radius*math.cos(2.0*math.pi*i/facets)
    y = radius*math.sin(2.0*math.pi*i/facets)
    print "    [ %.5f, %.5f, %.5f ]" % (x, y, -0.5)

print "    [ %.5f, %.5f, %.5f ]" % (0.0, 0.0, -0.5)

# points at the other end
for i in range(facets):
    x = radius*math.cos(2.0*math.pi*i/facets)
    y = radius*math.sin(2.0*math.pi*i/facets)
    print "    [ %.5f, %.5f, %.5f ]" % (x, y, 0.5)

print "    [ %.5f, %.5f, %.5f ]" % (0.0, 0.0, 0.5)

print "  }\n  ELEMENT_CONNECTIONS {"

# triangles on the side of the cylinder
for i in range(facets-1):
    print "    [ %d, %d, %d ]" % (i, i+1, i+facets+1)
    print "    [ %d, %d, %d ]" % (i+1, i+facets+2, i+facets+1)

# completing the side
print "    [ %d, %d, %d ]" % (facets-1, 0, 2*facets)
print "    [ %d, %d, %d ]" % (0, facets+1, 2*facets)

# one end cap
for i in range(facets-1):
    print "    [ %d, %d, %d ]" % (i+1, i, facets)

# completing the end
print "    [ %d, %d, %d ]" % (0, facets-1, facets)

# the other end cap
for i in range(facets-1):
    print "    [ %d, %d, %d ]" % (i+facets+1, i+facets+2, 2*facets+1)

# completing the other end
print "    [ %d, %d, %d ]" % (2*facets, facets+1, 2*facets+1)

print "  }\n}"

