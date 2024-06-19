//#pragma GCC optimize "trapv"
#include "precomputed_roots.h"
#include "ntheory.h"
#include <quadmath.h>


PrecomputedRoots::PrecomputedRoots () {
    power_bound = 0;
    prime_bound = 0;
}

PrecomputedRoots::PrecomputedRoots (
    __int128_t _power_bound,
    __int128_t _prime_bound,
    const std::vector<__int128_t> &primelist
) {
    //root_table = new cppmap[__int128_t, __int128_t]()
    __int128_t p;
    unsigned int max_exponent;
    __int128_t current_root;
    //std::cout << "created";
    //primes_gen = primerange(prime_bound)
    prime_bound = _prime_bound;
    power_bound = _power_bound;
    if (prime_bound == 0)
        prime_bound = (__int128_t) sqrtq(power_bound);
    // Populate the root table with roots modulo the maximally possible 
    // power of all primes congruent to 1 mod 3.  For the remaining ones,
    // we don't need a table.
    for (unsigned int i=0; i < primelist.size(); i++) {
        p = primelist[i];
        if (p > prime_bound)
            break;
        if (p % 3 != 1)
            continue; 
        max_exponent = (unsigned int) (log2(power_bound) / log2(p));
        current_root = prim_root_of_unity_mod_pl(p, max_exponent);
        root_table[p] = current_root;
    }
}

std::vector<__int128_t> PrecomputedRoots::roots_mod_pl (
    __int128_t p,
    unsigned int l
) const {
    std::vector<__int128_t> roots;
    __int128_t q;
    __int128_t current_root;
    if (p == 3)
        return roots_of_unity_mod_3l(l);
    roots.push_back(1);
    if (p % 3 != 1)
        return roots;
    if (p < prime_bound){
        q = long_pow(p, l);
        current_root = root_table.at(p) % q;
        roots.push_back(current_root);
        roots.push_back((current_root * current_root) % q);
        return roots;
    }
    return roots_of_unity_mod_pl(p, l);
}
