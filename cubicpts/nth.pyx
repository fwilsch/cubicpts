"""
Provides number theoretic helper functions.
"""


from sympy import jacobi_symbol, trailing

from libc.math cimport sqrt, llround
from libcpp.list cimport list as cpplist
from libcpp.vector cimport vector

cimport cython


cdef extern from *:
    ctypedef int int128 "__int128_t"
    ctypedef int uint128 "__uint128_t"

@cython.cdivision(True)
cpdef long long gcd(long long a, long long b) nogil:
    """Returns the gcd of a and b."""
    while a != 0:
        a, b = b%a, a
    if b < 0:
        return -b
    return b


@cython.cdivision(True)
cpdef long long mod_inverse(long long a, long long m) nogil:
    """
    Inverts a modulo m.
    
    Args:
        a: The residue to invert. Has to be coprime to m.
        m: The modulus, a non-zero integer.

    Returns:
        An integer between -m and m that is an inverse to a modulo m.
    """
    cdef long long m_0, t
    cdef cpplist[long long] quotients
    cdef long long factor_prev = 0
    cdef long long factor_cur = 1
    cdef long long factor_new

    if m == 1 or m == -1:
        return 0
    m_0 = m
    a = a % m
    while a != 0:
        t = a
        a = m % t
        quotients.push_back(m//t)
        m = t
    quotients.pop_back()
    for i in range(0, quotients.size()):
        factor_new = factor_prev - quotients.back() * factor_cur
        #factor_cur, factor_prev = factor_new, factor_cur
        factor_prev = factor_cur
        factor_cur = factor_new
        quotients.pop_back()
    if m < 0:
        return - (factor_cur % m_0)
    return factor_cur % m_0


@cython.cdivision(True)
cdef long long crt(long long a_1, long long a_2, long long m_1,
                   long long m_2) nogil:
    """
    Solves simultaneuous congruences for coprime moduli.
    """
    cdef long long m_1_bar, m_2_bar
    m_1_bar = mod_inverse(m_1, m_2)
    m_2_bar = mod_inverse(m_2, m_1)
    return (m_2_bar*m_2*a_1 + m_1_bar*m_1*a_2) % (m_1*m_2)


@cython.cdivision(True)
cdef vector[long long] crt_list(vector[long long] moduli, vector[vector[long long]] elements):
    """
    Solves simultaneuous congruences for coprime moduli.

    Args:
        moduli: A vector of integers, the moduli m_1, ..., m_s
        elements: A vector of length s. Elements are vectors of integers to be
                  lifted.
    """
    # This (or its usage) can probably be optimized to leverage the fact that we
    # often want to lift 1.
    cdef vector[long long] inverse_moduli, lengths, cur_index, solutions
    cdef vector[vector[long long]] lifted_elements
    cdef vector[long long] current_lifts
    cdef long long product = 1
    cdef unsigned int length = moduli.size()
    cdef unsigned int i, j
    cdef long long inverse_modulus, cur_solution
    cdef bint sols_not_exhausted = True

    if length == 0:
        return solutions

    #solutions.push_back(elements.size())
    #return solutions

    # initialize some variables
    for i in range(0, length):
        product *= moduli[i]
        lengths.push_back(elements[i].size())
        cur_index.push_back(0)


    # Compute inverses of m_i modulo the other moduli.  Use this to solve
    # the congruences k = k_i (m_i), k = 0 (m_j) for j != i.
    for i in range(0, length):
        other_moduli = product // moduli[i]
        inverse_modulus = mod_inverse(other_moduli, moduli[i])
        inverse_moduli.push_back(inverse_modulus)
        current_lifts.clear()
        for j in range(0, elements[i].size()):
            current_lifts.push_back(other_moduli*inverse_modulus*elements[i][j])
        lifted_elements.push_back(current_lifts) # Check whether this really copies
    
    #return lengths
    # The solutions are now all sums of the found lifts.
    while sols_not_exhausted:
        cur_solution = 0
        for i in range(0, length):
            cur_solution = (cur_solution + lifted_elements[i][cur_index[i]]) % product
        if cur_solution < 0:
            cur_solution += product
        solutions.push_back(cur_solution)

        # Continue to the next element:
        for i in range(0, length):
            if cur_index[i] < lengths[i]-1:
                cur_index[i] += 1
                break
            else:
                cur_index[i] = 0 
                if i == length-1:
                    sols_not_exhausted = False
                    #return solutions
                    
    return solutions


@cython.cdivision(True)
cpdef vector[long long] crt2_list(long long modulus_1, long long modulus_2,
                                  vector[long long] elements_1,
                                  vector[long long] elements_2) nogil:
    """
    Compute numbers reducing to an element of elements_i mod modulus_i
    for i=1,2 where the two moduli must be coprime.
    """
    cdef long long inverse_1 = mod_inverse(modulus_1, modulus_2)
    cdef long long total_modulus = modulus_1 * modulus_2
    cdef long long lift_2
    cdef vector[long long] solutions
    for i in range(0, elements_1.size()):
        for j in range(0, elements_2.size()):
            lift_2 = ((elements_2[j] + modulus_2 - elements_1[i]) * inverse_1) % modulus_2
            solutions.push_back((elements_1[i] + lift_2*modulus_1)% total_modulus)
    return solutions


def py_crt_list(moduli, elements):
    """
    Wrapper (unnecessarily complicated).
    """
    cdef vector[long long] c_moduli, cur_elements, c_solutions
    cdef vector[vector[long long]] c_elements
    cdef unsigned int i
    for m in moduli:
        c_moduli.push_back(<long long > m)
    for element in elements:
        cur_elements.clear()
        for e in element:
            cur_elements.push_back(<long long> e)
        c_elements.push_back(cur_elements)

    c_solutions = crt_list(c_moduli, c_elements)
    solutions = []
    for i in range(0, c_solutions.size()):
        solutions.append(c_solutions[i])
    return solutions


@cython.boundscheck(False)
cpdef bint is_square(int128 arg) nogil:
    """
    Checks whether an integer is a square.

    Args:
        arg: Integer < 2^100
    """
    cdef int128 i_sqr
    if arg < 0:
        return False
    return u_is_square(<uint128> arg)


@cython.boundscheck(False)
cpdef bint u_is_square(uint128 arg) nogil:
    """
    Checks whether an unsigned integer is a square.

    Args:
        arg: Unsigned integer < 2^100
    """
    cdef uint128 i_sqr
    # Test whether arg is a square mod 8.  Other moduli don't seem to yield a
    # speed improvement.
    if arg & <uint128> 2 == <uint128> 2:
        return False
    if arg & <uint128> 5 == <uint128> 5:
        return False
    i_sqr = <uint128> (sqrt(<double> arg) + <double> 0.5)
    if i_sqr*i_sqr == arg:
        return True
    return False


def kronecker_symbol(a,b):
    """Return the Kronecker symbol (a/b)."""
    if b<0:
        return ValueError
    if a == 1:
        return 1
    if gcd(a,b) != 1:
        return 0
    k_2 = trailing(b)
    factor_2 = 1 #the part of the result depending on the power of 2 dividing b
    if k_2 % 2 == 1 and a % 8 in (3,5):
        factor_2 = -1
    b = b // (2**k_2)
    return factor_2 * jacobi_symbol(a, b)


def dirichlet_character(a):
    """Return a list containing the values of the Dirichlet character (a/b)
    defined via the Kronecker symbol."""
    if a == 1:
        return [1, 1, 1, 1]
    chi = []
    for b in range(0, abs(4*a)):
        chi.append(kronecker_symbol(a, b))
    return chi