from libcpp.list cimport list as cpplist
from libcpp.vector cimport vector
from libcpp.map cimport map as cppmap
from libcpp.set cimport set as cppset
from libcpp.string cimport string

cdef extern from *:
    ctypedef int int128 "__int128"
    ctypedef int uint128 "__uint128_t"

cdef extern from "cubicpts/diagonal/io128.h":
    string to_string(int128 n);

cdef extern from "cubicpts/diagonal/ntheory.h":
    bint is_fourth_power(int128 n);
    bint is_square(int128 n);
    int128 pow_mod(int128 base, unsigned int exp, int128 mod);
    int128 long_pow(int128 base, unsigned int exp);
    int128 prim_root_of_unity_mod_p(int128 p);
    int128 prim_root_of_unity_mod_pl(int128 p, unsigned int l);
    vector[int128] roots_of_unity_mod_3l(unsigned int l);
    vector[int128] roots_of_unity_mod_pl(int128 p, unsigned int l);
    int128 mod_inverse(int128 a, int128 m);

cdef extern from "cubicpts/diagonal/points_diagonal.h":
    vector[vector[int128]]  points(
        int128 bound,
        int128 a
    );
    vector[vector[int128]]  points_batch(
        int128 bound,
        int128 a
    );
    vector[vector[int128]] points_batch_partial(
        int128 bound,
        int128 a_max,
        int128 largest_prime_lower,
        int128 largest_prime_upper
    );
    vector[vector[int128]] points_batch_partial(
        int128 bound,
        int128 a_max,
        int128 largest_prime_lower,
        int128 largest_prime_upper,
        bint print
    );
    bint on_accumulating_curve(
        int128 x,
        int128 y,
        int128 z,
        int128 a
    )

cdef extern from "cubicpts/diagonal/precomputed_roots.h":
    cdef cppclass PrecomputedRoots:
        PrecomputedRoots(
            int128 _power_bound,
            int128 _prime_bound,
            vector[int128] &primelist
        )
        vector[int128] roots_mod_pl(int128 p, unsigned int l)

cdef extern from "cubicpts/diagonal/numbergenerator.h":
    int128 get_da_bound(int128 bound, int128 a_max);
    cdef cppclass NumberGenerator:
        int128 get_value();
        NumberGenerator(int128 bound);
        vector[int128]* roots();
        cpplist[int128] divisors(int128 min, int128 max);
        void compute_roots();
        bint next();
