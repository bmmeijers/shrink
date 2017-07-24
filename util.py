"""
Utility functions.
"""
import math


def compute_envelope(points):
    """
    For a set of 2d points compute the bounding box
    (to be able to give a good initial big radius)
    """
    result = [[float('inf'), float('inf')],
              [float('-inf'), float('-inf')]]
    for pt in points:
        for i in range(2):
            result[0][i] = min(pt[i], result[0][i])
            result[1][i] = max(pt[i], result[1][i])
    return result


def distance2(a, b):
    """
    Return squared distance between two 2d points a and b
    """
    dist = 0.0
    for i in range(2):
        dist += pow(b[i] - a[i], 2)
    return dist


# -- debugging
def circle(center, radius):
    ct = 540
    alpha = math.pi * 2 / ct
    circ = []
    for n in range(ct):
        pt = center[0] + math.cos(n * alpha) * radius, center[1] + math.sin(n * alpha) * radius
        circ.append(pt)
    circ.append(circ[0])
    return circ


def format_circle(pts):
    return "POLYGON(({0}))".format(",".join(["{0[0]} {0[1]}".format(pt) for pt in pts]))


def format_point(pt):
    return "POINT({0[0]} {0[1]})".format(pt)


def format_line(start, end):
    return "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})".format(start, end)
