"""
Approximate the 2D medial axis given a set of points sampled from a
shape and the normal vectors to the shape at those points

The original 3D algorithm is described in:

3D medial axis point approximation using nearest neighbors and the normal field
Ma, J., Bae, S.W. & Choi, S.
Vis Comput (2012) 28: 7.
https://doi.org/10.1007/s00371-011-0594-7
"""

import math
from kdtree import KdTree
from normals import compute_normals
from vectorops import mul, add, norm, sub, dot
from util import distance2, compute_envelope

# -- for visual debugging
from util import  circle, format_circle, format_point, format_line
import random
import time
import sys


# -- some important calculations for the algorithm
def cos_angle(p, q):
    """Calculate the cosine of angle between vector p and q
    See http://en.wikipedia.org/wiki/Law_of_cosines#Vector_formulation
    """
    cos_angle = dot(p, q) / mul(norm(p), norm(q))
    # clamp value in range [-1, 1]
    return max(min(1, cos_angle), -1)


def compute_radius(p, n, q):
    """
    Compute *signed* radius of ball that touches points p and q and 
    is centered on along the normal n of p.

    Raises ValueError in case point p === q

    Note: Can return negative value!
    """
    if p == q:
        raise ValueError("For equal points we can *not* compute a radius")
    d = norm(sub(p, q))
    cos_theta = dot(n, sub(p, q)) / d
    if cos_theta == 0:
        raise ValueError('Vector (n) and (p-q) are perpendicular')
    radius = d / (2.0 * cos_theta)
    return radius


# -- structures
class Ball(object):
    """
    A ball with its center near the axis of the shape
    """
    def __init__(self, center, radius, p, q, iterations):
        self.center = center
        self.radius = radius
        self.p = p
        self.q = q
        self.iterations = iterations


# -- the shrinking ball algorithm
def shrink_balls(ring, balls, inner = True):
    """
    Approximate the 2D medial axis by shrinking 'balls'
    """

    # steps
    # -----
    # 1. load data
    #    * transform edges to (point, normal) list:
    #    -> every corner as (point, normal) and every edge also half way
    #    * index (point) list with KdTree
    #    * compute min/max box (needed for determining large radius)
    # 2. run algorithm to construct Balls: 
    #    (center point, radius, 2 touching points, # iterations)
    # 3. output Balls

    # FIXME take proper polygon as input, not a ring 
    # (or take topology structure and face id)
    points, normals = compute_normals(ring)

    # -- VISUAL DEBUG with QGIS, see qgiswatch.py
    debug = True
    if debug:
        with open("/tmp/normals.wkt", 'w') as fh:
            fh.write("wkt\n")
            for p, n in zip(points, normals):
                fh.write(format_line(p, add(p, mul(n, 10))) + "\n")
                del p, n
        del fh
    # -- END

    tree = KdTree(points[:])
    envelope = compute_envelope(points)

    # set big radius based on envelope of points
    bigradius = max([envelope[1][i] - envelope[0][i] for i in range(2)]) * 2

    stop_eps = 1.0
    radius = bigradius
    radius_prev = None
    q_prev = None
    for p, n in zip(points, normals):
        # use previous close candidate to reduce number of steps 
        # for convergence, see sec. 3.2 of (Ma, 2012)
        if radius != bigradius and p != q_prev:
            try:
                radius = compute_radius(p, n, q_prev)
            except ValueError:
                radius = bigradius
        else:
            radius = bigradius

        # in case we have a previous second point (q_prev)
        # and compute_radius gets a negative radius
        # then we start the algorithm with the big radius
        if radius <= 0:
            radius = bigradius

        radius = bigradius

        # find inner balls or outer balls?
        if not inner:
            n = mul(n, -1.0) # swap normal direction 
                             # in case we want to find balls
                             # at other side of shape (exterior)

        # find the smallest ball possible given point p and normal n
        for _ in xrange(len(points)):  # this should be the max iter. needed
            # *center*       ball its center point in this iteration
            # *q*            2nd point that defines a ball (with p and n)
            # *radius*       radius of the ball in this iteration
            # *radius_prev*  radius found in the previous iteration
            # *q_prev*       2nd point from previous iteration

            # compute ball center
            center = sub(p, mul(n, radius))

            # find 2 closest points to the center
            candidates = tree.k_nearest(center, k=2)
            q = candidates[0][1]
            if q == p:
                # if we found the p as closest point we take the other point
                q = candidates[1][1]
                # but we make sure that this point still lies at correct side
                # of halfplane through p
                # otherwise radius would grow
                if (candidates[1][0] - radius) > 1e-5:
                    b = Ball(center, radius, p, q_prev, _)
                    balls.append(b)
                    break
            q_prev = q

            # -- VISUAL DEBUG with QGIS
            if False:
                one = add(p, mul(n, -3 * bigradius))
                other = add(p, mul(n, 3 * bigradius))
                with open("/tmp/line.wkt", 'w') as fh:
                    fh.write("wkt\n")
                    fh.write(format_line(one, other) + "\n")
                del one
                del other
                with open("/tmp/ptp.wkt", 'w') as fh:
                    fh.write("wkt\n")
                    fh.write(format_point(p) + "\n")
                with open("/tmp/ptq.wkt", 'w') as fh:
                    fh.write("wkt\n")
                    fh.write(format_point(q) + "\n")
                with open("/tmp/center.wkt", 'w') as fh:
                    fh.write("wkt\n")
                    fh.write(format_point(center) + "\n")
                with open("/tmp/ball.wkt", 'w') as fh:
                    ballstr = format_circle(circle(center, radius))
                    fh.write("wkt\n")
                    fh.write(ballstr + "\n")
                    del ballstr
                # -- let QGIS refresh, once file has changed...
                with open("/home/martijn/signal", "w") as fh:
                    fh.write("{0}".format(random.randint(0, int(1e6))))
                    # wait a bit to let QGIS refresh
                    time.sleep(0.5)
                # -- ask user to press key
                #raw_input('pause')
                del fh
            # -- END VISUAL DEBUG

            # store old radius and compute new radius
            radius_prev = radius
            radius = compute_radius(p, n, q)
            # terminate when radius of new ball is not significantly different
            if abs(radius - radius_prev) < stop_eps:
                # print >> sys.stderr, "Found ball with", _, "step(s)"
                b = Ball(center, radius, p, q, _)
                balls.append(b)
                break
    return balls

def _test():
    from data import poly as ring
    inner = []
    shrink_balls(ring, inner, True)
    outer = []
    shrink_balls(ring, outer, False)

    if True:
        # FIXME: split in inner and outer set of balls
        with open("/tmp/balls.wkt", 'w') as fh:
            #
            fh.write("wkt;iters;inner\n")
            lines = "\n".join(format_circle(circle(b.center, b.radius)) + ";" + str(b.iterations) + ";True" for b in inner)
            lines += "\n"
            lines += "\n".join(format_circle(circle(b.center, b.radius)) + ";" + str(b.iterations) + ";False" for b in outer)
            fh.write(lines)

        with open("/tmp/centers.wkt", 'w') as fh:
            #
            fh.write("wkt;iters;inner\n")
            lines = "\n".join(format_point(b.center) + ";" + str(b.iterations)+ ";True" for b in inner)
            lines += "\n"
            lines += "\n".join(format_point(b.center) + ";" + str(b.iterations)+ ";False" for b in outer)
            fh.write(lines)

if __name__ == "__main__":
    _test()
