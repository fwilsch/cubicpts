#ifndef _CUBICPTS_NTHEORY
#define _CUBICPTS_NTHEORY

#include <vector>
#include <map>
#include <set>
#include <list>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <exception>

__int128_t pow_mod (__int128_t base, unsigned int exp, __int128_t mod);
__int128_t long_pow (__int128_t base, unsigned int exp);
bool is_square (__int128_t n);
bool is_square (__uint128_t n);
bool is_cube (__int128_t n);
bool is_fourth_power(__int128_t n);
__int128_t prim_root_of_unity_mod_p (__int128_t p);
__int128_t prim_root_of_unity_mod_pl (__int128_t p, unsigned int l);
std::vector<__int128_t> roots_of_unity_mod_3l (unsigned int l);
std::vector<__int128_t> roots_of_unity_mod_pl (__int128_t p, unsigned int l);
__int128_t mod_inverse (__int128_t a, __int128_t m);
std::vector<__int128_t> crt2_list (
    __int128_t modulus_1,
    __int128_t modulus_2,
    std::vector<__int128_t> const &elements_1,
    std::vector<__int128_t> const &elements_2
);

#endif //_CUBICPTS_NTHEORY