#ifndef _CUBICPTS_NUMBERGENERATOR
#define _CUBICPTS_NUMBERGENERATOR

#include <vector>
#include <map>
#include <set>
#include <list>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <exception>
#include <memory>
#include <primesieve.hpp>
#include "precomputed_roots.h"


class NumberGenerator{
    private:
        __int128_t value, length;
        std::vector<__int128_t> primes;
        std::set<unsigned int> prime_factors;
        std::map<unsigned int, unsigned int> exponents;
        const __int128_t bound, largest_prime_max;
        std::map<unsigned int, __int128_t> previous_values;
        // Roots of unity modulo $\\prod_{k>i} p_k^{e_k}$
        std::map< unsigned int, std::vector<__int128_t> > all_roots;
        std::shared_ptr<PrecomputedRoots> pre_roots;
        primesieve::iterator largest_prime;
        void increment_prime_exponent(unsigned int index, bool compute_roots=true);
        void compute_roots ();
    public:
        __int128_t get_value() const;
        NumberGenerator (__int128_t bound);
        NumberGenerator (__int128_t bound,
            __int128_t largest_prime_min,
            __int128_t largest_prime_max
        );
        ~NumberGenerator ();
        std::vector<__int128_t>* roots ();
        bool next ();
        void divisors (
            std::list<__int128_t> &divs,
            __int128_t min = 0,
            __int128_t max = -1
        ) const;
        std::list<__int128_t> divisors (
            __int128_t min = 0,
            __int128_t max = -1
        ) const;
};

#endif //_CUBICPTS_NUMBERGENERATOR