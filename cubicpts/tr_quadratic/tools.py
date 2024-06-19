"""
Auxilliary functions for the surfaces treated by tr_quadratic.
"""
from itertools import product
from math import pi, log

from sympy import isprime, legendre_symbol, primefactors

from cubicpts import tr_quadratic as quad
from cubicpts.tools import reducepoint

def filter_quadratic(pointlist, a=-2, b=3, c=-3, d=5):
    """
    Filter out points on quadratically embedded lines.
    """
    result = []
    for (x, y, z) in pointlist:
        if -z not in (
                a*(c*d*y + c + d),
                b*(c*d*y + c + d),
                c*(a*b*x + a + b),
                d*(a*b*x + a + b)
        ):
            result.append((x, y, z))
    return result


def filter_cubic_ex(pointlist):
    """Filter out points on cubics in xyz=q(x,y) with a,b,c,d=-2,3,-3,5"""
    results = []
    points_c1 = []
    points_c2 = []
    points_c3 = []
    points_c4 = []
    for (x, y, z) in pointlist:
        if 2700*x == z**3 + 5*z**2 - 240*z and 450*y == -z**2 - 5*z + 150:
            points_c1.append((x, y, z))
            continue
        if -z**3 + 3*z**2 + 144*z == 1620*x and z**2 - 3*z - 54 == y:
            points_c2.append((x, y, z))
            continue
        if 270*x == -z**2 - 6*z + 135 and 4050*y == z**3 + 6*z**2 - 225*z:
            points_c3.append((x, y, z))
            continue
        if 180*x == z**2 - 4*z - 60 and 2700*y == -z**3 + 4*z**2 + 150*z:
            points_c4.append((x, y, z))
            continue
        results.append((x, y, z))
    print(points_c1)
    print(points_c2)
    print(points_c3)
    print(points_c4)
    return results


def filter_quartic_ex(pointlist):
    """Filter out points on cubics in xyz=q(x,y) with a,b,c,d=-2,3,-3,5"""
    results = []
    points_c1 = []
    points_c2 = []
    points_c3 = []
    points_c4 = []
    for (x, y, z) in pointlist:
        if 2700*x == z**3 + 5*z**2 - 240*z:
            points_c1.append((x, y, z))
            continue
        if -z**3 + 3*z**2 + 144*z == 1620*x:
            points_c2.append((x, y, z))
            continue
        if 4050*y == z**3 + 6*z**2 - 225*z:
            points_c3.append((x, y, z))
            continue
        if 2700*y == -z**3 + 4*z**2 + 150*z:
            points_c4.append((x, y, z))
            continue
        results.append((x, y, z))
    return results


def test_reduction(p, plist, a=-2, b=3, c=-3, d=5):
    local_points = set(quad.points_mod_p(p, a, b, c, d))
    for point in plist:
        local_points.discard(reducepoint(point, p))
    return local_points



def get_faces():
    """
    Return generators of effective cones for alpha.  Cf. the paper.
    """
    rays = [[1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 1], [1, -1, 0, 0, - 1, 0, 0],
            [1, -1, 0, 0, 0, -1, 0], [1, -1, 0, 0, 0, 0, -1],
            [1, 0, -1, 0, -1, 0, 0], [1, 0, -1, 0, 0, -1, 0],
            [1, 0, -1, 0, 0, 0, -1], [1, 0, 0, -1, -1, 0, 0],
            [1, 0, 0, -1, 0, -1, 0], [1, 0, 0, -1, 0, 0, -1],
            [1, -1, -1, -1, 0, 0, 0], [1, 0, 0, 0, -1, -1, -1],
            [-3, 1, 1, 1, 1, 1, 1]
           ]

    face1, face2, face3 = [], [], []
    for r in rays:
        rn1 = [r[2], r[3], r[5], r[6]]
        face1.append(rn1)

        rn2 = [r[0] + r[3] + r[4],
               r[2] - r[3],
               r[5],
               r[6]
              ]
        face2.append(rn2)

        rn3 = [r[0] + r[2]+ r[4],
               r[3] - r[2],
               r[5] - r[4],
               r[6] - r[4]
              ]
        face3.append(rn3)

    return(face1, face2, face3)


def testtwomods(p1, p2, k1, k2, pointlist):
    """
    Test strong approximation modulo two prime powers.

    Return a list of all pairs (P,Q) of points P modulo p1^k1 and
    Q modulo p2^k2 which are not simultaneously lift to a point in pointlist.
    """
    loc_points1 = set(quad.points_mod(p1**k1))
    loc_points2 = set(quad.points_mod(p2**k2))
    loc_comb = set(product(loc_points1, loc_points2))
    points_red = set()
    for pt in pointlist:
        points_red.add((reducepoint(pt, p1**k1), reducepoint(pt, p2**k2)))
    return loc_comb.difference(points_red)


def complete_euler(a=-2, b=3, c=-3, d=5):
    '''
    Return the Euler product for xyz=q(x,y).
    '''
    ceuler = 36/(pi**4)
    for p in primefactors(a*b*c*d):
        if a*b % p == 0 and c*d % p == 0:
            ceuler *= (1 + 1/(p*p)) / ((1 + 1/p)**2)
        else:
            ceuler *= (1 + 1/p + 1/(p*p)) / ((1 + 1/p)**2)
    return ceuler


def euler_l(l):
    '''
    Return the Euler product for lxyz=q(x,y) with a,b,c,d=-2,3,-3,5.

    Arguments:
        l: 1 or prime not equal to 2,3,5,7,19 (default: 1)

    Raises:
        ValueError: If requirements on l are not satisfied.
    '''
    c2 = 7/9
    c3 = 5/8
    c5 = 31/36
    if l == 1:
        cl = 1
    elif isprime(l) and 2*3*5*7*19 % l != 0:
        cl = (1 - legendre_symbol(-10, l)/l) / (1 + 2/l + 1/(l*l))
    else:
        raise ValueError("l has to be 1 or a prime different from 2,3,5,7,19")
    return c2*c3*c5*cl* 36 / (pi**4)


def testmodh(p, k1, k2):
    """
    Return all points modulo p^k1 that do not lift to p^k2.
    """
    points1 = set(quad.points_mod(p**k1))
    points2red = set()
    for pt in quad.points_mod(p**k2):
        points2red.add(reducepoint(pt, p**k1))
    return points1.difference(points2red)


def testmod(p, k1, k2, pointlist):
    """
    Return all classes mod p^k1 that lift to p^k2 but not to pointlist.
    """
    loc_points = set()
    for pt in quad.points_mod(p**k2):
        loc_points.add(reducepoint(pt, p**k1))

    points_red = set()
    for pt in pointlist:
        points_red.add(reducepoint(pt, p**k1))

    return loc_points.difference(points_red)
