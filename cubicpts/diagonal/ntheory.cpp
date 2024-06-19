//#pragma GCC optimize "trapv"
#include "ntheory.h"
#include <quadmath.h>


const __int128_t QUADMATH_THRESHOLD = ((__int128_t) 2) << 106;
const __uint128_t U_QUADMATH_THRESHOLD = ((__uint128_t) 2) << 106;


bool is_square(__int128_t n) {
    __int128_t candidate;
    if (n < 0)
        return 0;
    if ((n & (__int128_t) 2) == (__int128_t) 2)
        return 0;
    if ((n & (__int128_t) 5) == (__int128_t) 5)
        return 0;

    if (n < QUADMATH_THRESHOLD)
        candidate = (__int128_t) (sqrt(n));
    else    
        candidate = (__int128_t) (sqrtq(n));

    if (candidate * candidate == n)
        return 1;
    return 0;  
}

bool is_square(__uint128_t n) {
    __uint128_t candidate;
    if ((n & (__uint128_t) 2) == (__uint128_t) 2)
        return 0;
    if ((n & (__uint128_t) 5) == (__uint128_t) 5)
        return 0;
    if (n < U_QUADMATH_THRESHOLD)
        candidate = (__uint128_t) (sqrt(n));
    else
        candidate = (__uint128_t) (sqrtq(n));

    if (candidate * candidate == n)
        return 1;
    return 0;  
}

bool is_cube(__int128_t n){
    __int128_t candidate;
    if (n < 0)
        n = -n;

    if (n < QUADMATH_THRESHOLD)
        candidate = (__int128_t) pow(n, 1/3);
    else
        candidate = (__int128_t) powq(n, 1/3);
        
    if (candidate * candidate * candidate == n)
        return true;
    return false;

}

bool is_fourth_power(__int128_t n){
    __int128_t candidate, candidate_squared;
    if (n < 0)
        return false;

    if (n < QUADMATH_THRESHOLD)
        candidate = (__int128_t) pow(n, 0.25);
    else
        candidate = (__int128_t) powq(n, 0.25);
        
    candidate_squared = candidate*candidate;
    if (candidate_squared * candidate_squared == n)
        return true;
    return false;

}


// Raise n to the power of k modulo mod.
__int128_t pow_mod(__int128_t base, unsigned int exp, __int128_t mod) {  
    __int128_t solution = 1;
    __int128_t cur_power;
    if (mod <= 0)
        return 0;
    if (exp == 0)
        return 1;
    
    cur_power = base % mod;
    if (cur_power <= 0)
        cur_power += mod;
    
    if ((exp & 1) == 1)
        solution  = (solution * cur_power) % mod;   
    exp >>= 1; 

    // Repeated squaring
    while (exp != 0) { 
        cur_power = (cur_power * cur_power) % mod;
        if ((exp & 1) == 1)
            solution  = (solution * cur_power) % mod;
        exp >>= 1; 
    }
    return (__int128_t) solution;
}

// Raise base to the power exp. No overflow checking.
__int128_t long_pow(__int128_t base, unsigned int exp){ 
    __int128_t solution = 1;
    __int128_t cur_power;
    
    if (exp == 0)
        return 1;
    cur_power = base;
    
    // Repeated squaring. Only adjust power after the first iteration.
    if ((exp & 1) == 1)
        solution  *= cur_power;
    exp >>= 1;  

    while (exp != 0){
        cur_power *= cur_power;
        if ((exp & 1) == 1)
            solution  *= cur_power;
        exp >>= 1;
    }
    return solution;
}

__int128_t prim_root_of_unity_mod_p(__int128_t p){
    __int128_t test_exp, root;
    if ((p - 1) % 3 != 0)
        return 0;
    test_exp = (p - 1) / 3;
    for (int n=2; n<p; ++n){
        root = pow_mod(n, test_exp, p);
        if (root != 1)
            return root;
    }
    return 0;
}

// Return a primitive root of unity modulo p^l, where p is 1 modulo 3
__int128_t prim_root_of_unity_mod_pl (__int128_t p, unsigned int l) {
    __int128_t q_max = long_pow(p, l);
    __int128_t q_cur = p;
    __int128_t root_mod = (__int128_t) prim_root_of_unity_mod_p(p);
    __int128_t derivative, derivative_inverse;
    while (q_cur < q_max){
        q_cur *= q_cur;
        if (q_cur >= q_max)
            q_cur = q_max;

        // Newton: x_new = x - f(x)/f'(x) =  - (x^3-1) / 3x^2
        derivative = (3 * root_mod * root_mod) % q_cur;
        derivative_inverse = mod_inverse(derivative, q_cur);
        root_mod = (root_mod
                    - (pow_mod(root_mod, 3, q_cur) - 1) * derivative_inverse
                   // - (root_mod * root_mod * root_mod - 1) * derivative_inverse
                   ) % q_cur; // possibly reduce this expression after squaring
    }
    return (__int128_t) root_mod;
}

std::vector<__int128_t> roots_of_unity_mod_3l (unsigned int l) {
    std::vector<__int128_t> roots;
    roots.push_back(1);
    if (l == 1)
        return roots;
    roots.push_back(1 + long_pow(3, l-1));
    roots.push_back(1 + 2 * long_pow(3, l-1));
    return roots;
}

//Return all roots of unity modulo p^l.
std::vector<__int128_t> roots_of_unity_mod_pl (__int128_t p, unsigned int l) {
    std::vector<__int128_t> roots;
    __int128_t q, root_mod;
    if (p == 3)
        return roots_of_unity_mod_3l(l);
    roots.push_back(1);
    if ((p - 1) % 3 != 0)
        return roots;
    q = long_pow(p, l);
    root_mod = prim_root_of_unity_mod_pl(p, l);

    roots.push_back(root_mod);
    roots.push_back((root_mod * root_mod)% q);
    return roots;
}

//Inverts a modulo m.
//
//Args:
//    a: The residue to invert. Has to be coprime to m.
//    m: The modulus, a non-zero integer.
//
//Returns:
//    An integer between -m and m that is an inverse to a modulo m.
__int128_t mod_inverse (__int128_t a, __int128_t m) {
    __int128_t m_0, t;
    std::list<__int128_t> quotients;
    // TODO: check precision here. Maybe increase the size of the following
    // integers or reduce in the last loop.
    __int128_t factor_prev = 0;
    __int128_t factor_cur = 1;
    __int128_t factor_new;
    if (m == 1 || m == -1)
        return 0;
    m_0 = m;
    a = a % m;
    while (a != 0) {
        t = a;
        a = m % t;
        quotients.push_back(m/t);
        //DEBUG
        //std::cout << (long long) m << " "
        //          << (long long) t << " "
        //          << (long long) a << " "
        //          << (long long) (m/t) << "\n";
        m = t;
    }
    //std::cout << quotients.size() << "\n";
    quotients.pop_back();
    while (quotients.size() != 0) {
        factor_new = factor_prev -  factor_cur * quotients.back();
        //std::cout << quotients.size() << " "
        //          << (long long) factor_prev << " "
        //          << (long long) factor_cur << " "
        //          << (long long) factor_new << "\n";
        factor_prev = factor_cur;
        factor_cur = factor_new;
        quotients.pop_back();
    }
    if (m < 0)
        return - (factor_cur % m_0);
    return factor_cur % m_0;
}

// Compute numbers reducing to an element of elements_i mod modulus_i
// for i=1,2 where the two moduli must be coprime.
std::vector<__int128_t> crt2_list (
    __int128_t modulus_1,
    __int128_t modulus_2,
    std::vector<__int128_t> const &elements_1,
    std::vector<__int128_t> const &elements_2
) {  
    __int128_t inverse_1 = mod_inverse(modulus_1, modulus_2);
    __int128_t total_modulus = modulus_1 * modulus_2;
    __int128_t lift_2;
    __int128_t new_root;
    std::vector<__int128_t> solutions;
    for (__int128_t el_1 : elements_1) {
        for (__int128_t el_2 : elements_2) {
            lift_2 = ((el_2 - el_1) * inverse_1) % modulus_2;
            new_root = (el_1 + lift_2*modulus_1) % total_modulus;
            if (new_root < 0)
                new_root += total_modulus;
            solutions.push_back(
                new_root
            );
        }
    }
    return solutions;
}
