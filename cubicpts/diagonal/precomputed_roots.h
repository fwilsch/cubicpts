#ifndef _CUBICPTS_PRECOMPUTED_ROOTS
#define _CUBICPTS_PRECOMPUTED_ROOTS

#include <vector>
#include <map>
#include <set>
#include <list>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <exception>

class PrecomputedRoots{
    private:
        std::map<__int128_t, __int128_t> root_table;
        __int128_t power_bound;
        __int128_t prime_bound;
    public:
        PrecomputedRoots();
        PrecomputedRoots(
            __int128_t _power_bound,
            __int128_t _prime_bound,
            const std::vector<__int128_t> &primelist
        );
        //Return a vector of roots of unity modulo p^l
        std::vector<__int128_t> roots_mod_pl(__int128_t p, unsigned int l) const;
};

#endif //_CUBICPTS_PRECOMPUTED_ROOTS