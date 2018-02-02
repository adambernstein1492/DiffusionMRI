#!/usr/bin/env python

import math

def write(output, facets=64, radius=0.5):
    # facets = 64
    #radius = 0.5
    output.write("structure POLYGON_LIST {\n  VERTEX_LIST {\n")

    # points at one end
    for i in range(facets):
        x = radius*math.cos(2.0*math.pi*i/facets)
        y = radius*math.sin(2.0*math.pi*i/facets)
        output.write("    [ %.5f, %.5f, %.5f ]\n" % (x, y, -0.5))

    output.write("    [ %.5f, %.5f, %.5f ]\n" % (0.0, 0.0, -0.5))

    # points at the other end
    for i in range(facets):
        x = radius*math.cos(2.0*math.pi*i/facets)
        y = radius*math.sin(2.0*math.pi*i/facets)
        output.write("    [ %.5f, %.5f, %.5f ]\n" % (x, y, 0.5))

    output.write("    [ %.5f, %.5f, %.5f ]\n" % (0.0, 0.0, 0.5))

    output.write("  }\n  ELEMENT_CONNECTIONS {\n")

    # triangles on the side of the cylinder
    for i in range(facets-1):
        output.write("    [ %d, %d, %d ]\n" % (i, i+1, i+facets+1))
        output.write("    [ %d, %d, %d ]\n" % (i+1, i+facets+2, i+facets+1))

    # completing the side
    output.write("    [ %d, %d, %d ]\n" % (facets-1, 0, 2*facets))
    output.write("    [ %d, %d, %d ]\n" % (0, facets+1, 2*facets))

    # one end cap
    for i in range(facets-1):
        output.write("    [ %d, %d, %d ]\n" % (i+1, i, facets))

    # completing the end
    output.write("    [ %d, %d, %d ]\n" % (0, facets-1, facets))

    # the other end cap
    for i in range(facets-1):
        output.write("    [ %d, %d, %d ]\n" % (i+facets+1, i+facets+2, 2*facets+1))

    # completing the other end
    output.write("    [ %d, %d, %d ]\n" % (2*facets, facets+1, 2*facets+1))

    output.write("  }\n}\n")
