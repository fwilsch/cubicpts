#ifndef _CUBICPTS_POINTS_DIAGONAL
#define _CUBICPTS_POINTS_DIAGONAL

#include <vector>
#include <map>
#include <set>
#include <list>
#include <cmath>
#include <algorithm>
#include <exception>
#include <memory>
#include "precomputed_roots.h"
#include "numbergenerator.h"


// Computes an upper bound for d*a that is used throughout.
__int128_t get_da_bound(__int128_t bound, __int128_t a_max);
bool yields_solution (__int128_t d, __int128_t z, __int128_t a);
bool on_accumulating_curve(__int128_t x, __int128_t y, __int128_t z,
                           __int128_t a);

std::vector< std::vector<__int128_t> > compute_solutions (
    __int128_t d,
    __int128_t z,
    __int128_t a
);

std::vector< std::vector<__int128_t> >  points (
    __int128_t bound,
    __int128_t a
);

/* Return a vector of vectors [x, y, z, a] such that each is a solution of 
$ax^3 + by^3 + z^3 = 1$, where a(x+y) is captured by gen.*/
std::vector< std::vector<__int128_t> > points_batch (
    std::unique_ptr<NumberGenerator> &gen,
    __int128_t bound,
    __int128_t a_max,
    bool print = false
);

/* Return a vector of vectors [x, y, z, a] such that each is a solution of 
$ax^3 + by^3 + z^3 = 1$*/
std::vector< std::vector<__int128_t> > points_batch (
    __int128_t bound,
    __int128_t a_max,
    bool print = false
);

std::vector< std::vector<__int128_t> > points_batch_partial (
    __int128_t bound,
    __int128_t a_max,
     __int128_t largest_prime_lower,
    __int128_t largest_prime_upper,
    bool print = false
);

#endif //_CUBICPTS_POINTS_DIAGONAL