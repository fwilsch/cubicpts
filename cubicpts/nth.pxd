from libcpp.vector cimport vector

cdef extern from *:
    ctypedef int int128 "__int128_t"
    ctypedef int uint128 "__uint128_t"

cpdef long long gcd(long long a, long long b) nogil
cpdef long long mod_inverse(long long a, long long m) nogil
cdef long long crt(long long a_1, long long a_2, long long m_1,
                   long long m_2) nogil
cdef vector[long long] crt_list(vector[long long] moduli,
                                vector[vector[long long]] elements)
cpdef vector[long long] crt2_list(long long modulus_1, long long modulus_2,
                                  vector[long long] elements_1,
                                  vector[long long] elements_2) nogil
cpdef bint is_square(int128 arg) nogil
cpdef bint u_is_square(uint128 arg) nogil