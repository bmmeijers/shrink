from vectorops import make_vector, unit, rotate90cw, add, mul


def segment_normal(end, start):
    """
    Return normal to a segment defined by two points (end, start)
    The normal is a unit vector (normalized by its length)
    """
    return rotate90cw(unit(make_vector(end, start)))


def compute_normals(ring):
    """
    Make normals perpendicular to a ring of vertices
    The normals are created on the end points and half way the vertices
    """
    # output lists
    normals = []
    points = []
    # create normals half way each segment
    ct = len(ring)
    for i in xrange(ct - 1):
        cur, nxt = ring[i], ring[i+1]
        n = segment_normal(nxt, cur)
        center = mul(add(cur, nxt), 0.5)
        normals.append(n)
        points.append(center)
    # create normals on every point, using normals for every segment
    ct = len(normals)
    for i in xrange(ct):
        cur, nxt = normals[i], normals[(i+1) % ct]
        n = unit(add(cur, nxt))
        pt = ring[i+1]
        normals.append(n)
        points.append(pt)
    return points, normals


if __name__ == "__main__":
    # ring = [(0,0), (10,0), (5,10), (0,0)]
    from data import poly as ring
    points, normals = compute_normals(ring)
    print "wkt"
    for p, n in zip(points, normals):
        print "LINESTRING({0[0]} {0[1]}, {1[0]} {1[1]})".format(p, add(p, n))

