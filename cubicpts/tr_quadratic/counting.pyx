"""
Finds points on affine cubic surfaces U of Picard rank 2 defined by equations 
of the form
q(x,y)=xyz, where q(x, y) = ab x^2 + (a+b)x + cdy^2 + (c+d)y + 1.
"""

from sympy import divisors
from cython.parallel cimport prange

cimport cython
cimport openmp
from libcpp.list cimport list as cpplist
from cubicpts.nth cimport mod_inverse, crt, gcd


cdef extern from *:
    ctypedef int int128 "__int128_t"


@cython.cdivision(True)
def points(long long bound, long long a=-2, long long b=3, long long c=-3,
           long long d=5):
    """
    Returns all solutions to xyz = q(x,y) of height <= bound.

    Here, q(x,y) = ab x^2 + (a+b)x + cdy^2 + (c+d)y + 1,
    e.g., xyz = - 6x^2 - 15y^2 + x + 2y + 1 (default).
    Exclude points on one of the lines defined over the integers.

    Args:
        bound: The maximal height up to which to search.
        a, b, c, d: The parameters of the surface.

    Returns:
        A list of 3-tuples (x, y, z) of integers that are solutions
        to the equation.
    """
    cdef long long x_0, x_1, x_2, y_0, y_1, y_2, t, k
    cdef int128 x, y, z
    cdef long long bound_k_l, bound_k_u
    cdef int128 q_x, q_xy#, z_l
    ts = divisors(c-d)
    sols = []
    cdef int128 bad_z1 = -(a*d+b*c)
    cdef int128 bad_z2 = -(a*c+b*d)
    cdef int128 cd = <int128> c*d
    cdef int128 c_p_d = <int128> c+d

    # Factor x = t x_1 x_2, (x_1, x_2)=1 to find y with
    # x_1 | (cy+1), x_2 | (dy+1).  This is a necessary condition for solubility,
    # and sufficient if t=1.
    for t in ts:
        for x_1 in prange(-bound//t, bound//t + 1, nogil = True,
                          schedule='static', chunksize = 1):
            if x_1 == 0 or gcd(x_1, c) != 1:
                continue
            for x_2 in range(1, bound//abs(t*x_1) + 1):
                # TODO: This is inefficient; we could generate them in a way
                # that guarantees they are coprime, probably saving a log here.
                # This might reduce the leading constant.
                if gcd(x_2, d) != 1 or gcd(x_1, x_2) != 1:
                    continue
                x_0 = x_1 * x_2
                x = x_0 * t
                y_1 = mod_inverse(-c, x_1)
                y_2 = mod_inverse(-d, x_2)
                y_0 = crt(y_1, y_2, x_1, x_2)
                # Solve the congruences mod t x_1, t x_2.  We know there
                # is a solution.  Since t is small, we just use brute force.
                #for _ in range(0, t):
                #    if (c*y_0+1) % (x_1*t) == 0 and (d*y_0+1) % (x_2*t) == 0:
                #        break
                #    y_0 = y_0 + x_1*x_2
                #assert (c*y_0+1) % (x_1*t) == 0 and (d*y_0+1) % (x_2*t) == 0
                q_x = <int128> (a*x+1)*(b*x+1)
                if x > 0:
                    bound_k_l = (-bound-y_0) // x_0
                    bound_k_u = (bound -y_0) // x_0
                else:
                    bound_k_l = (bound -y_0) // x_0
                    bound_k_u = (-bound-y_0) // x_0
                #for y in range(y_0-bound-x_0, y_0+bound+x_0, x_0):
                for k in range(bound_k_l - 1, bound_k_u + 1):
                    y = <int128> (y_0 + k*x_0)
                    if y == 0:
                        continue
                    if (q_x % y != 0):
                        continue
                    q_xy = q_x + cd*y*y + c_p_d*y
                    if (q_xy % (<int128> x*y) != 0):
                        continue
                    z = (q_xy // (x*y))
                    if (max(
                            abs(<long long> x),
                            abs(<long long> y),
                            abs(<long long> z)
                        )<= bound
                        and z not in (bad_z1, bad_z2)
                        ):
                        with gil:
                            sols.append((x, y, z))
    return set(sols)


def points_mod_p(long long p, long long a=-2, long long b=3, long long c=-3,
                 long long d=5):
    """
    Returns all solutions to xyz = q(x,y) modulo a prime (power).
    """
    if a*b*c*d*(a-b)*(c-d) % p == 0:
        raise NotImplementedError
    results = []
    cdef long long x_1, x_2, y_1, y_2, q_x, z

    # Find all solutions with x, y invertible.
    for x in range(1, p):
        q_x = (a*b*x*x + (a+b)*x) % p
        for y in range(1,p):
            inverse_xy = mod_inverse(x*y, p)
            z = ((q_x + (c*y+1)*(d*y+1)) * inverse_xy) % p
            results.append((x, y, z))
    
    # Find all sol's with y=0: (ax+1)(bx+1)=0, and analogously for x=0.
    x_1 = mod_inverse(-a, p) % p
    x_2 = mod_inverse(-b, p) % p
    y_1 = mod_inverse(-c, p) % p
    y_2 = mod_inverse(-d, p) % p
    for z in range(0, p):
        results.append((0, y_1, z))
        results.append((0, y_2, z))
        results.append((x_1, 0, z))
        results.append((x_2, 0, z))
    
    return results
                

@cython.cdivision(True)
def points_mod(int pk, int a=-2, int b=3, int c=-3, int d=5, int l=1):
    """
    Returns all solutions to l xyz = q(x,y) modulo pk.
    
    Here, q(x,y) = ab x^2 + (a+b)x + cdy^2 + (c+d)y + 1,
    e.g. xyz = - 6x^2 - 15y^2 + x + 2y + 1 (default).

    Arguments:
        pk: The modulus, a positive integer.
        a, b, c, d, l: The parameters of the surface.  l has to be non-zero.
    
    Returns:
        A list of integer tuples (x,y,z) of solutions modulo pk,
        with 0 <= x, y, z <= pk-1.
    """
    cdef int x, y, z
    cdef int q_x = a*b
    cdef int l_x = a+b
    cdef int q_y = c*d
    cdef int l_y = c+d
    cdef int poly_xy
    sols=[]

    if l == 0:
        raise ValueError("l has to be non-zero.")
    if pk <= 0:
        raise ValueError("pk has to be positive.")

    for x in prange(pk, nogil=True, schedule='dynamic'):
        for y in range(pk):
            poly_xy = ((q_x*x + l_x)*x + (q_y*y + l_y)*y + 1) % pk
            for z in range(pk):
                if (poly_xy - l*x*y*z) % pk == 0:
                    with gil:
                        sols.append((x, y, z))

    return sols