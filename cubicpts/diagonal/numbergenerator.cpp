//#pragma GCC optimize "trapv"
#include <quadmath.h>
#include <primesieve.hpp>
#include "numbergenerator.h"
#include "ntheory.h"
#include "precomputed_roots.h"


NumberGenerator::NumberGenerator(__int128_t bound) : 
    NumberGenerator::NumberGenerator(bound, 1, bound) { }


NumberGenerator::NumberGenerator (
    __int128_t bound,
    __int128_t largest_prime_min,
    __int128_t largest_prime_max
) : value(1), bound(bound), largest_prime_max(largest_prime_max){
    __int128_t smoothness = std::min({
        (__int128_t) sqrtq(bound),
        bound / largest_prime_min,
        largest_prime_max
    });
    primesieve::generate_primes(smoothness, &primes);
    std::sort(primes.begin(), primes.end());
    largest_prime.jump_to(std::max(smoothness + 1, largest_prime_min));
    primes.push_back(largest_prime.next_prime());
    pre_roots = std::make_shared<PrecomputedRoots>(bound, smoothness, primes);
    length = primes.size();
    if (largest_prime_min > 1) {
        auto large_prime = std::upper_bound(
            primes.begin(), primes.end(), largest_prime_min-1
        );
        unsigned int prime_index = std::distance(primes.begin(), large_prime);
        increment_prime_exponent(prime_index);
    }
}


void NumberGenerator::increment_prime_exponent(unsigned int index, bool recompute_roots) {
    if (prime_factors.find(index) == prime_factors.end()) {
        prime_factors.insert(index);
        exponents[index] = 1;
    } else {
         exponents[index] +=1;
    }
    value *= primes[index];
    if (recompute_roots)
        compute_roots();
}

NumberGenerator::~NumberGenerator () {}

__int128_t NumberGenerator::get_value() const {
    return value;
}

std::vector<__int128_t>* NumberGenerator::roots () {
    unsigned int i = *(prime_factors.begin());
    return &(all_roots[i]);
}

void NumberGenerator::compute_roots () {
    // Compute the roots of unity modulo self, using the precomputed ones of a
    // divisor of self. Use the fact that only the smallest prime can have
    // changed since the last computation.
    unsigned int smallest_index, next_index;
    __int128_t smallest_prime;
    unsigned int smallest_exponent;
    __int128_t moduli[2];

    if (prime_factors.size() == 0)
        return;
    smallest_index = *(prime_factors.begin());
    smallest_prime = primes[smallest_index];
    smallest_exponent = exponents[smallest_index];
    if (prime_factors.size() == 1) {
        // Nothing is precomputed.
        previous_values[smallest_index] = value;
        all_roots[smallest_index] = pre_roots->roots_mod_pl(
            smallest_prime,
            smallest_exponent
        );
        return;
    }
    // Use the precomputed roots of unity and combine them with 
    // the data modulo p^e for the current prime
    // Start by retrieving the precomputed data
    next_index = *(++prime_factors.begin());
    moduli[0] = previous_values[next_index];
    moduli[1] = long_pow(smallest_prime, smallest_exponent);
    all_roots[smallest_index] = crt2_list(
        moduli[0],
        moduli[1],
        all_roots[next_index], //previous roots
        pre_roots->roots_mod_pl(
            smallest_prime,
            smallest_exponent
        )
    );
    previous_values[smallest_index] = moduli[0] * moduli[1];
    return;
}

bool NumberGenerator::next () {
    // Change self to next number. Returns False if there is no new number
    // in the range; otherwise, returns True.
    __int128_t quotient = bound / value;
    __int128_t cur_power = 1;
    for (int i=0; i < length; ++i){
        if (primes[i] <= quotient) {
            // In this case, we multiply with this number once; then we're at
            // the next number in the list and compute everything.
            increment_prime_exponent(i);
            return 1;
        }
        else if (prime_factors.find(i) != prime_factors.end()) {
            // In this case, primes[i] is a factor of value, but we can't
            // multiply by it anymore. Hence, we divide by the largest possible
            // power of the prime.
            // We erase associated data, and then continue in the loop to
            // increase the next larger prime.
            cur_power = long_pow(primes[i], exponents[i]);
            value = value / cur_power;
            prime_factors.erase(i);
            previous_values.erase(i);
            exponents.erase(i);
            all_roots.erase(i);
            quotient = bound / value;
            if (i == length - 1) {
                // Step the largest prime by rewriting the last entry in primes.
                primes[length-1] = largest_prime.next_prime();
                if (primes[length-1] > largest_prime_max)
                    return 0; // We are done
                // Continue in the loop with i = length-1 in the next step to 
                // increase the exponent of the new prime and compute everything.
                i = length - 2;
                continue; 
            }
        }
        else {
        // We are at a prime that does not divide the current value and is too
        // large to multiply with.  Hence, we can move on to the next prime appearing
        // in the current value to divide by it.
            if (value != 1) {
                i = *(prime_factors.begin()) - 1;
                continue;
            }
        return 0;  
        }
    }
    // Reached the end of the generator
    return 0;
}

// Write a list of all divisors of self between min and max into divs.
// (Clear divs before.)
void NumberGenerator::divisors (
    std::list<__int128_t> &divs,
    __int128_t min,
    __int128_t max
) const {
    divs.clear();
    std::map<unsigned int, unsigned int> div_exponents;
    bool finished = false;
    if (max < 0)
        max = bound;
    __int128_t div = 1;
    if (1 >= min && 1 <= max)
        divs.push_back(1);

    while (!finished) {
        // TODO: Optimize this to look up stuff in the set less often.
        // Could write everything into a vector at the beginning
        finished = true;
        for (unsigned int i : prime_factors) {
            if (div_exponents.count(i) == 0) {
                div_exponents[i] = 1;
                div *= primes[i];
                if (div >= min && div <= max)
                    divs.push_back(div);
                finished = false;
                break;
            }
            else if (div_exponents[i] < exponents.at(i)) { // and div_exponents[i] < 3
                div_exponents[i] += 1;
                div *= primes[i];
                if (div >= min && div <= max)
                    divs.push_back(div);
                finished = false;
                break;
            }
            else if (div_exponents[i] == exponents.at(i)) {
                div /= long_pow(primes[i], div_exponents[i]);
                div_exponents.erase(i);
            }
        }
    }
}

// Return a list of all divisors between min and max of self.
std::list<__int128_t> NumberGenerator::divisors (__int128_t min, __int128_t max) const {  
    std::list<__int128_t> divs;
    divisors(divs, min, max);
    return divs;
}