"""
Provides tools to analyze the number of points of bounded height.
"""

from math import log

from numpy import unique


def pointheight(P):
    return max(abs(P[0]), abs(P[1]), abs(P[2]))


def reducepoint(P, mod):
    """
    Reduces the point P modulo mod.
    """
    return (P[0]%mod, P[1]%mod, P[2]%mod)


def import_pointlist(filename):
    """
    Returns a list of points imported from a file.

    Args:
        filename: Name of a file containing points: one point per line,
            coordinates separated by spaces.
    """
    pointlist = []
    for line in open(filename):
        values = line.split()
        pointlist.append((int(values[0]), int(values[1]), int(values[2])))
    return pointlist


def height_statistics(pointlist):
    """
    Count points of bounded height.

    Returns a sorted list of tuples (h, N), where h runs over all heights
    attained by points in pointlist, and N is the number of points of
    height <= h.
    """
    heights = [max([abs(p[0]), abs(p[1]), abs(p[2])]) for p in pointlist]
    t = unique(heights, return_counts=True)
    tstat = list(zip(t[0], t[1]))
    N = 0
    stat = []
    for (h, n) in tstat:
        N += n
        stat.append((h, N))
    return stat


def asymptotic(height_stat, const=1, exp=4):
    """
    Compares a height statistic to an expected formula of the form c (log B)^e.

    Args:
        height_stat: height statistic as returned by height_statistics().
        const: expected constant of the asymptotic formula. (default: 1)
        exp: expected exponent of log B in the asymptotic formula. (default: 4)

    Returns:
        A list of tuples (h, N, C), where (h, N) are as in height_stat and
        C is the quotient of N by the expected value.
    """
    asymp = []
    for (h, N) in height_stat:
        asymp.append((h, N, N/(const*log(h)**exp)))
    return asymp


def write_points(pointlist, filename):
    """
    Write the list of points to a file.
    """
    ofile = open(filename, 'w')
    for (x, y, z) in pointlist:
        print(x, y, z, file=ofile)


def write_stat(pointlist, filename):
    """
    Write the height statistics generated from a list of points to a file.
    """
    ofile = open(filename, 'w')
    stat = height_statistics(pointlist)
    for (a, b) in stat:
        print(log(a), b, file=ofile)


