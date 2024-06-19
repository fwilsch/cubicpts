"""
Computes points on diagonal cubic surfaces of the form ax^3 + ay^3 + z^3 = 1.
"""
cimport cython
from libc.math cimport pow, sqrt, log2, log
from libcpp.list cimport list as cpplist
from libcpp.vector cimport vector
from libcpp.map cimport map as cppmap
from libcpp.set cimport set as cppset
from cython.operator cimport dereference as deref, preincrement

cimport cubicpts.diagonal.wrapper as dw

cdef extern from *:
    ctypedef int int128 "__int128"
    ctypedef int uint128 "__uint128_t"

cdef class PyPrecomputedRoots:

    cdef dw.PrecomputedRoots* cPrecomputedRoots

    def __init__(self, power_bound, primelist):
        self.cPrecomputedRoots = new dw.PrecomputedRoots(power_bound, 0, primelist)
    
    def __dealloc__(self):
        del self.cPrecomputedRoots
    
    def roots_mod_pl(self, int128 p, unsigned int l):
        return self.cPrecomputedRoots.roots_mod_pl(p, l)


cdef class PyNumberGenerator:

    cdef dw.NumberGenerator* cnumbergenerator

    def __init__(self, int128 bound):
        self.cnumbergenerator = new dw.NumberGenerator(bound)

    def __dealloc__(self):
        del self.cnumbergenerator

    def next(self):
        return self.cnumbergenerator.next()

    def value(self):
        return self.cnumbergenerator.get_value()

    def roots(self):
        return deref(self.cnumbergenerator.roots())

    def divisors(self, min, max):
        return self.cnumbergenerator.divisors(min, max)


cpdef str _to_string(int128 n):
    return dw.to_string(n).decode('UTF-8')

cpdef bint _is_fourth_power(int128 n):
    return dw.is_fourth_power(n)

cpdef bint is_square(int128 n):
    return dw. is_square(n)

cpdef bint on_accumulating_curve(int128 x, int128 y, int128 z, int128 a):
    return dw.on_accumulating_curve(x, y, z, a)

cpdef int128 _pow_mod(int128 base, unsigned int exp, int128 mod):
    return dw.pow_mod(base, exp, mod)

cpdef int128 _long_pow(int128 base, unsigned int exp):
    return dw.long_pow(base, exp)

cpdef int128 _prim_root_of_unity_mod_p(int128 p):
    return dw.prim_root_of_unity_mod_p(p)

cpdef int128 _prim_root_of_unity_mod_pl(int128 p, unsigned int l):
    return dw.prim_root_of_unity_mod_pl(p, l)

cpdef vector[int128] _roots_of_unity_mod_3l(unsigned int l):
    return dw.roots_of_unity_mod_3l(l)

cpdef vector[int128] _roots_of_unity_mod_pl(int128 p, unsigned int l):
    return dw.roots_of_unity_mod_pl(p, l)

cpdef int128 _mod_inverse(int128 a, int128 m):
    return dw.mod_inverse(a, m)

cpdef int128 get_da_bound(bound, a_max):
    """Get the bound for d*a=a*(x+y) used in the computations."""
    return dw.get_da_bound(bound, a_max)

cpdef points(int128 bound, int128 a = 2):
    solutions = dw.points(bound, a)
    return [(P[0], P[1], P[2]) for P in solutions]


cpdef points_batch(int128 bound, int128 a_max = 10):
    solutions = dw.points_batch(bound, a_max)
    return [((P[0], P[1], P[2]), P[3]) for P in solutions]


cpdef points_interval(int128 bound, int128 a_max, int128 p_min, int128 p_max, bint _print = False):
    """
    Determine solutions whose largest prime factor is in [p_min, p_max).
    """
    return dw.points_batch_partial(
        bound,
        a_max,
        p_min,
        p_max,
        _print
    )


cpdef int128 divide_interval(int128 interval, int step, int total_steps):
    """
    Subdivide [0, interval) into total_steps smaller intervals.

    total_steps needs to be smaller than interval.
    """
    cdef double quotient = <double> step / total_steps
    if step == 0:
        return 0
    # constant length: interval*quotient
    # length proportional to b: interval**quotient
    # length proportional to log b 
    #return <int128> (interval*quotient * log(step)/log(total_steps))
    #return <int128> (interval*quotient*quotient * (1.0 + log(quotient)/log(interval))**2)
    #return <int128> pow(interval, sqrt(quotient))
    return <int128> pow(interval, quotient)
    #return <int128> (pow(interval, quotient)*(1.0 - log(quotient)))


cpdef points_batch_partial(int128 bound, int128 a_max, int step, int total_steps, bint _print = False):
    """
    Determine solutions to $ax^3 + ay^3 + z^3 = 1$ for several values of $a$.
    
    If this function is run for a fixed value bound, a_max, and total_steps and
    all step in [0, total_steps), then the results contain all solutions with
    |x|, |y|, |z| <= bound and all cubefree 1 < a <= a_max (and some more).
    """
    cdef int128 da_bound
    da_bound = dw.get_da_bound(bound, a_max)

    p_min = divide_interval(da_bound, step, total_steps)
    p_max = divide_interval(da_bound, step+1, total_steps)
    solutions = points_interval(bound, a_max, p_min, p_max, _print)

    return solutions