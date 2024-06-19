#include "io128.h"
#include <iostream>
#include <algorithm>

#define POWER_OF_TEN 1000000000000000
#define EXPONENT_OF_TEN 15

std::string to_string(__int128_t n) {
    int sign = 1;
    std::string result, digit_string;
    long long digits;
    //return std::string(EXPONENT_OF_TEN, '0');
    if (n < 0) {
        sign = -1;
        n = -n;
    }
    do {
        digits = (long long) (n % (__int128_t) POWER_OF_TEN); 
        n /= (__int128_t) POWER_OF_TEN;
        digit_string = std::to_string((long long) digits);
        if (n != 0)
            digit_string = std::string(
                EXPONENT_OF_TEN - std::min(EXPONENT_OF_TEN, (int) digit_string.length()), '0'
            ) + digit_string;
        result = digit_string + result;
    } while (n != 0);
    if (sign == -1)
        result = '-' + result;
    return result;
}

void print_vector(
    std::vector<__int128_t> vector
) {
    for (__int128_t x : vector) {
        std::cout << to_string(x) << " ";
    }
    std::cout << "\n";
}

//std::ostream& operator<< (std::ostream& os, __int128_t n) {
//    std::string result;
//    int length;
//    if (n < 0) {
//        os << "-";
//        n = -n;
//    }
//    return os;
//}