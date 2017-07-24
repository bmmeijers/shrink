import operator
import logging
from util import distance2, compute_envelope, format_line
import unittest

class KdNode(object):
    """
    A node in the KdTree

    Attributes:
        items: the items (vertices) stored in this node
        active: list with booleans, indicating whether items are active
        left: the left child node
        right: the right child node
        cutdim: the axis that was used for sorting while creating this node
        cutval: the median value for the split
    """
    __slots__ = ('_items', 'left', 'right', 'parent', 'cutdim', 'cutval',
                 'active')

    def __init__(self, items=None, left=None, right=None):
        """
        Constructor.

        Note, all properties can be set later and by default are set to None.
        """
        self.active = None
        self.items = items
        self.left = left
        self.right = right
        self.cutdim = None  # cut on dimension
        self.cutval = None  # median value for the cut

    @property
    def items(self):
        """The items that belong with this node."""
        return self._items

    @items.setter
    def items(self, value):
        """
        Set the items that belong with this node

        Args:
            value: The items that have to go with this node.

        Returns:
            None.
        """
        if value is None:
            self._items = None
            self.active = None
        else:
            self._items = value
            self.active = [True] * len(self._items)


class KdTree(object):
    """
    A KdTree

    This implementation is based on:
    Jon Louis Bentley, K-d Trees for Semi-dynamic Point Sets (DOI_)
    .. _DOI: http://dx.doi.org/10.1145/98524.98564

    Attributes:
        root: The KdNode instance that forms the root of the tree.

    Note: The list with points L is *not* copied and the algorithm in the
    constructor is sorting the list, so make a copy (L[:]) if you rely 
    on the order of the point list later.
    """
    def __init__(self, L, axes=(0, 1), cutoff=10):
        """
        Constructor
        """
        self.root = KdNode()
        stack = [(L, self.root, 0)]
        while stack:
            items, node, axis = stack.pop()
            if len(items) <= cutoff:
                node.items = items
            else:
                items.sort(key=operator.itemgetter(axis))
                halfway = len(items) // 2
                lft, rgt = items[:halfway], items[halfway:]
                nodelft, nodergt = KdNode(), KdNode()
                node.left = nodelft
                node.right = nodergt
                node.cutdim = axes[axis]
                node.cutval = items[halfway][axes[axis]]
                next_axis = (axis + 1) % len(axes)
                stack.append((lft, nodelft, next_axis))
                stack.append((rgt, nodergt, next_axis))

    def intersection(self, envelope):
        """
        Returns the items that interact with the given axis aligned box.

        Args:
            envelope: Axis-aligned box that forms the query region
                (includes boundaries of the box)
                Has to be given as: [(xmin, ymin), (xmax, ymax)]

        Returns:
            List with items that intersect with the query region
        """
        stack = [self.root]
        result = []
        while stack:
            node = stack.pop()
            if node.items:
                for item, active in zip(node.items, node.active):
                    if not active:
                        continue
                    ok = True
                    for axis in xrange(2):
                        c = item[axis]
                        low = envelope[0][axis]
                        high = envelope[1][axis]
                        if c < low or c > high:
                            ok = False
                            break
                    if ok:
                        result.append(item)
            else:
                axis = node.cutdim
                median = node.cutval
                low = envelope[0][axis]
                high = envelope[1][axis]
                if median <= high:
                    stack.append(node.right)
                if median >= low:
                    stack.append(node.left)
        return result

    def __contains__(self, item):
        """
        Returns bool whether item is *in* the tree.
        """
        _, active, _ = self._find(item)
        return active

    def _find(self, item):
        """
        Traverse the tree to find the given item

        Args:
            item: The object to be found

        Returns:
            A tuple. The tuple has three items:
            0. boolean, whether item is found,
            1. boolean, whether item is active,
            2. the KdNode instance in which the item is eventually stored,
               or None.
        """
        stack = [self.root]
        while stack:
            node = stack.pop()
            if node.items and item in node.items:
                return True, node.active[node.items.index(item)], node
            elif not node.items:
                axis = node.cutdim
                cur = item[axis]
                median = node.cutval
                if median >= cur:
                    stack.append(node.left)
                if median <= cur:
                    stack.append(node.right)
        return False, False, None

    def delete(self, item):
        """
        Deletes the item from the tree.

        Args:
            item: The object to be marked inactive in the tree

        Returns:
            A boolean indicating whether the item was deleted.
        """
        is_found, active, node = self._find(item)
        if is_found and active:
            idx = node.items.index(item)
            node.active[idx] = False
            return True
        else:
            return False

    def undelete(self, item):
        """
        Undeletes the item from the tree.

        Args:
            item: The object to be marked active in the tree

        Returns:
            A boolean indicating whether the item was undeleted.
        """
        is_found, active, node = self._find(item)
        if is_found and not active:
            idx = node.items.index(item)
            node.active[idx] = True
            return True
        else:
            return False

    def k_nearest(self, pt, k):
        """
        Search for the k nearest points to point pt.

        Returns a list with tuples, where a tuple represents: (distance, point).
        """
        if k < 1:
            raise ValueError('k should be at least 1')
        result = []
        visit_ct = k_nearest(self.root, pt, k, result)
        logging.debug('Visited {0} leaf nodes'.format(visit_ct))
        return result  # [item for (d, item) in result]


def visit_k_nearest(node, pt, k, result):
    """
    Visit node to see if there are close points to point pt in this node
    Only k close points are kept in the result list
    """
    # rather brute force but because cut off and k expected to be rather small
    # not further optimized
    # (result could instead of list be a bin heap with at most k items)
    for active, item in zip(node.active, node.items):
        # check active items
        if active:
            d = distance2(pt, item)
            result.append( (d, item) )
        # sort on distance
        result.sort(key=lambda x: x[0])
        # keep max k items
        while len(result) > k:
            result.pop()


def k_nearest(node, pt, k, result):
    """
    Recursive nearest neighbours implementation
    Adds candidate points to the result argument (initially an empty list)

    Returns: count of nodes visited during the search
    """
    if node.items:
        visit_k_nearest(node, pt, k, result)
        return 1
    else:
        dx = pt[node.cutdim] - node.cutval
        if dx <= 0:
            near = node.left
            far = node.right
        else:
            near = node.right
            far = node.left
        ct_near = k_nearest(near, pt, k, result)
        # check if we found results, 
        # if we have sufficient results and the closest of these
        # is closer than the split line, we do not have to search further
        if result and len(result) >= k and pow(dx, 2) >= result[0][0]:
            return ct_near 
        ct_far = k_nearest(far, pt, k, result)
        return ct_near + ct_far


def node2path(node, lowx, lowy, highx, highy, polygons, lines, points):
    """
    Output the KdTree as set of points / split lines / polygons

    Input: node, envelope of the node (lowx, ..., highy)
           polygons, lines, points (initially as empty lists)
    """
    if node.items:
        ll = lowx, lowy
        lr = highx, lowy
        ur = highx, highy
        ul = lowx, highy
        polygons.append((ll, lr, ur, ul))
        for pt in node.items:
            points.append(pt)
        return
    else:
        if (node.cutdim % 2 == 0):
            items.append( ((node.cutval, lowy), (node.cutval, highy)) )
            node2path(node.left, lowx, lowy, node.cutval, highy, items, points)
            node2path(node.right, node.cutval, lowy, highx, highy, items, points)
        else:
            items.append((( lowx, node.cutval),( highx, node.cutval)))
            node2path(node.left, lowx, lowy, highx, node.cutval, items, points)
            node2path(node.right, lowx, node.cutval, highx, highy, items, points)
        return 


def visualize(root, envelope):
    """
    Output a visualization of the KdTree (for display in QGIS)
    """
    polygons = []
    lines = []
    points = []
    node2path(root,
              envelope[0][0], envelope[0][1], envelope[1][0], envelope[1][1], 
              polygons, lines, points)
    with open('/tmp/splitlines.wkt', 'w') as fh:
        fh.write('wkt\n')
        for line in lines:
            fh.write(format_line(line[0], line[1]) + '\n')


class TestCase(unittest.TestCase):
    def test(self):
        """Tests."""
        # -- construct the tree, only first 2 dims are used for the tree construction
        L = range(100)
        L = [(i, i, i, i) for i in L]
        tree = KdTree(L)
        for item in [(i, i, i, i) for i in range(100)]:
            assert item in tree
        for item in [(i, i, i, i) for i in range(101, 150)]:
            assert item not in tree
        for item in [(i, i, i, i) for i in range(-1, -10, -1)]:
            assert item not in tree
        # -- query the tree
        expected = [(i, i, i, i) for i in range(50, 61)]
        assert tree.intersection([(50, 50), (60, 60)]) == expected
        # -- modify the tree and query again
        assert tree.delete((1, 1, 1, 1))
        assert (1, 1, 1, 1) not in tree
        assert not tree.delete((-1, -1, 1, 1))
        assert tree.intersection([(0, 0), (2, 2)]) == [(0, 0, 0, 0), (2, 2, 2, 2)]
        assert tree.undelete((1, 1, 1, 1))
        assert (1, 1, 1, 1) in tree
        assert not tree.undelete((1, 1, 1, 1))

    def test_k_nearest(self):
        """Tests for knn searching"""
        L = range(100)
        L = [(i, i, i, i) for i in L]
        tree = KdTree(L)
        # remove distance, only keep points from the result
        items = lambda items: [x for (d, x) in items] 
        assert items(tree.k_nearest((-1, -1), 1)) == [(0, 0, 0, 0)]
        assert items(tree.k_nearest((100, 100), 1)) == [(99, 99, 99, 99)]
        assert items(tree.k_nearest((50, 50), 1)) == [(50, 50, 50, 50)]
        assert items(tree.k_nearest((-1, -1), 2)) == [(0, 0, 0, 0),
                                                      (1, 1, 1, 1)]


def _test():
    import unittest
    unittest.main()


def _visualize():
    from datapts import pts
    kdt = KdTree(pts)
    envelope = compute_envelope(pts)
    visualize(kdt.root, envelope)


if __name__ == "__main__":
    _test()

