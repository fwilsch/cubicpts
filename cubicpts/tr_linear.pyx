"""
Finds points on surfaces defined by an equation of the form
(x^2 -ay^2) = l(x, y), where l is a linear form.
"""

from sympy import divisors
from cython.parallel cimport prange

cimport cython
cimport openmp

from libc.math cimport sqrt, llround
#from libquadmath import sqrtq

from cubicpts.nth cimport is_square, u_is_square

cdef extern from *:
    ctypedef int int128 "__int128_t"
    ctypedef int uint128 "__uint128_t"



#cdef uint128 I_1 = <uint128> 1, I_2 = <uint128> 2
#cdef uint128 I_3 = <uint128> 3, I_4 = <uint128> 4
#cdef uint128 I_5 = <uint128> 5
#cdef uint128 I_7 = <uint128> 7, I_0 = <uint128> 0,

#cpdef int128 perfect_square_root(int128 arg):
#    if arg < 0:
#        return 

#@cython.boundscheck(False)
#cpdef bint is_square(int128 arg) nogil:
#    """
#    Checks whether an integer is a square.
#
#    Args:
#        arg: Integer < 2^100
#    """
#    cdef int128 i_sqr
#    if arg < 0:
#        return False
#    return u_is_square(<uint128> arg)
#
#
#@cython.boundscheck(False)
#cpdef bint u_is_square(uint128 arg) nogil:
#    """
#    Checks whether an unsigned integer is a square.
#
#    Args:
#        arg: Unsigned nteger < 2^100
#    """
#    cdef uint128 i_sqr
#    # Test whether arg is a square mod 8.  Other moduli don't seem to yield a
#    # speed improvement.
#    if arg & <uint128> 2 == <uint128> 2:
#        return False
#    if arg & <uint128> 5 == <uint128> 5:
#        return False
#    i_sqr = <uint128> (sqrt(<double> arg) + <double> 0.5)
#    if i_sqr*i_sqr == arg:
#        return True
#    return False

@cython.cdivision(True)
def points(long long bound, long long a=2, long long c=3, long long d=-5,
           long long e=7):
    """"
    Returns all solutions to (x^2-ay^2)z = l(x,y) of height <= bound.
    
    Exclude those on the line defined by z = 0.  Here, l(x,y) = cx + dy+ e.
    
    For example, we find (x, y, z) with
    (x^2 - 2y^2)z = 3x - 5y + 7 (default).

    Args:
        bound: The maximal height up to which to search.
        a, c, d, e: The parameters of the surface.

    Returns:
        A list of 3-tuples (x, y, z) of integers that are solutions
        to the equation.
    """
    sols=[]
    cdef long long z, u, u_bound
    cdef int128 disc_red, disc_1
    cdef long long quad, lin, sq_disc
    cdef long long x, x_num, y, y_num

    if c == 0:
        return points_var(bound, a, d, e)

    quad = d*d-a*c*c # The coefficient of y^2 in an equation solved later on.

    for z in prange(-bound, bound+1, nogil = True, schedule='static',
                    chunksize = 1):
        if z == 0:
            continue
        # We have cx + dy + e = 0 (mod z), i.e., cx + dy + e = zu for some u.
        # Its absoulute value is bounded by
        u_bound = ((abs(c)+abs(d))*bound + abs(e)) // abs(z)
        for u in range(-u_bound, u_bound+1):
            # Substitute x in the equation, multiply by c^2:
            # ((dy+e-zu)^2-ac^2y^2)z - zuc^2 = 0
            # We divide by z and solve the quadratic equation
            # (d^2-ac^2)y^2 + 2d(e-zu)y + (e-zu)^2 - uc^2 = 0 for y.
            # Its discriminant is 4c^2 * disc_red with
            disc_1 = e-z*u
            disc_red = a*disc_1*disc_1 + u*quad

            if not is_square(disc_red):
                continue
            sq_disc = 2*c * <long long>(sqrt(<double> disc_red)+0.5)
            lin = 2*d*(e-z*u)
            y_num = -lin + sq_disc

            if y_num % (2*quad) == 0:
                y = y_num // (2*quad)
                x_num = -d*y - e + z*u
                if x_num % c == 0:
                    x = x_num // c
                    if abs(x) <= bound and abs(y) <= bound:
                        with gil:
                            sols.append((x, y, z))

            y_num = -lin - sq_disc
            if y_num % (2*quad) == 0:
                y = y_num // (2*quad)
                x_num = -d*y - e + z*u
                if x_num % c == 0:
                    x = x_num // c
                    if abs(x) <= bound and abs(y) <= bound:
                        with gil:
                            sols.append((x, y, z))
    
    return set(sols)


@cython.cdivision(True)
def points_var(long long bound, long long a=2, long long d=1, long long e=-1):
    """"
    Find solutions to (x^2-ay^2)z = dy + e of height <= bound.
    
    Returns all solutions except those on the line defined by z = 0.
    For example, find (x, y, z) with (x^2 - 2y^2)z = y - 1 (default).

    Args:
        bound: The maximal height up to which to search.
        a, c, d, e: The parameters of the surface.

    Returns:
        A set of 3-tuples (x, y, z) of integers that are solutions
        to the equation.
    """
    sols=[]
    cdef long long z, u, u_bound
    cdef int128 x_square
    cdef long long x, y, y_num

    for z in prange(-bound, bound+1, nogil = True, schedule='static',
                    chunksize = 1):
        if z == 0:
            continue
        # We have dy + e = 0 (mod z), i.e., dy + e = zu for some u.
        # Its absoulute value is bounded by
        u_bound = (abs(d)*bound + abs(e)) // abs(z)
 
        for u in range(-u_bound, u_bound+1): # Can be optimized
            # Now y = (zu-e)/d
            y_num = z*u - e
            if y_num % d != 0:
                continue
            y = y_num//d
            
            # Now solve the equation (x^2 - ay^2) z = uz for x
            x_square = <int128> a* y*y + u
            if not is_square(x_square):
                continue
            x = <long long> (sqrt(<double> x_square) + 0.5)

            if x <= bound and abs(y) <= bound:
                with gil:
                    sols.append(( x, y, z))
                    sols.append((-x, y, z))
    
    return set(sols)

            
@cython.cdivision(True)
def points_brute(long long bound, long long a=2, long long c=3, long long d=-5,
           long long e=7):
    """
    Compute the solutions to (x^2-ay^2)z = l(x,y) of height <= bound via brute force.
    """
    sols = []
    cdef long long x, y, z_num, z_den, z

    for x in range(-bound, bound+1):
        for y in range(-bound, bound+1):
            z_num = c*x + d*y + e
            z_den = x*x - a*y*y
            if z_den == 0 or z_num == 0:
                continue 
            if z_num % z_den == 0: 
                z = z_num//z_den
                if abs(z) <= bound:
                    sols.append((x, y, z))
    return sols