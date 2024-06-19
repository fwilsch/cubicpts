//#pragma GCC optimize "trapv"
#include <omp.h>
#include <quadmath.h>
//#include <iostream>
#include "points_diagonal.h"
#include "numbergenerator.h"
#include "ntheory.h"
#include "io128.h"

// Get an upper bound on d*a based on upper bounds for a and the height
__int128_t get_da_bound(__int128_t bound, __int128_t a_max) {
    return std::min(
        (__int128_t) pow(4.0* (bound*bound*bound+1) * a_max*a_max, 1.0/3),
        (__int128_t) 2*a_max*bound
    );
}

// Checks whether the combination of d = x+y and z yields a point on U_a.
bool yields_solution (__int128_t d, __int128_t z, __int128_t a) {
    __int128_t lhs, lhs_num, disc, constant_num, constant_term;
    __int128_t da;

    // Write 
    // (1 - z^3) / a = x^3 + y^3, i.e., 
    // (1 - z^3) / da = d^2 - 3dx + 3x^2,
    // and solve for x (where now d=x+y can carry a negative sign).
    lhs_num = (__int128_t)1- (__int128_t)z*z*z;
    da = d*a;
    if (lhs_num % da != 0) //TODO: include a in d
        return false;
    lhs = lhs_num / da;
    
    constant_num = (__int128_t) d*d - lhs;
    if (constant_num % 3 != 0)
        return false;
    constant_term = constant_num / 3;
    disc =  (__int128_t) d*d - 4*constant_term;
    if (! is_square(disc))
        return false;
    return true;
}

// Returns a vector of solutions that d = x+y and z correspond to.
std::vector< std::vector<__int128_t> > compute_solutions (
    __int128_t d,
    __int128_t z,
    __int128_t a
) {
    std::vector< std::vector<__int128_t> > solutions;
    std::vector<__int128_t> new_solution(3);
    // Repeat the computations in yield_solution.
    __int128_t lhs, lhs_num, disc, constant_num, constant_term;
    __int128_t x, y;
    __int128_t sqrt_disc, da;
    lhs_num = (__int128_t) 1- (__int128_t) z*z*z;
    da = d*a;
    if (lhs_num % da != 0) //TODO: include a in d
        throw std::invalid_argument("lhs is not divisible by da, this should have been checked");
        //return solutions;
    lhs = lhs_num / da;
    
    constant_num = (__int128_t) d*d - lhs;
    if (constant_num % 3 != 0)
        throw std::invalid_argument("constant_num is not divisible by 3, this should have been checked");
        //return solutions;
    constant_term = constant_num / 3;
    disc =  (__int128_t) d*d - 4*constant_term;
    if (! is_square(disc))
        throw std::invalid_argument("discriminant is not a square, this should have been checked");
        //return solutions;
    sqrt_disc = (__int128_t) sqrtq(disc);
    x = (d + sqrt_disc)/2;
    y = d - x;
    if (on_accumulating_curve(x, y, z, a))
        return solutions;
    new_solution[0] = (__int128_t) x;
    new_solution[1] = (__int128_t) y;
    new_solution[2] = (__int128_t) z;
    solutions.push_back(new_solution);
    if (sqrt_disc != 0) {
        new_solution[0] = (__int128_t) y;
        new_solution[1] = (__int128_t) x;
        solutions.push_back(new_solution);
    }
    return solutions;
}


/* Checks whether the point (x,y,z) is on one of the accumulating curves 
 * of degree 4, that is, z = 9at^3 +1, x = -9at^4-3t, y = 9at^4 (or x and y replaced)*/
bool on_accumulating_curve(
    __int128_t x,
    __int128_t y,
    __int128_t z,
    __int128_t a
) {
    __int128_t t, t_cubed, x_or_y;
    // t has to be -(x+y)/3. Check whether 2 of the other equalities hold; the 
    // third one is automatic.
    if ((x+y) % 3 != 0)
        return false;
    t = - (x+y) /3;
    t_cubed = t*t*t;
    if (z != 9*a*t_cubed + 1) 
        return false;
    x_or_y = -9*a * t*t_cubed - 3*t;
    if (x == x_or_y || y == x_or_y) 
        return true;
    return false;
}

/* Compute all points on ax^3 + ay^3 +z^3 = 1 with |x|, |y|, |z| <= bound 
(and some solutions not satisfying these constraints).
Output: vector of vectors of the form [x, y, z]. 
Exclude points with z=1 or x+y = 1. */
std::vector< std::vector<__int128_t> > points (
    __int128_t bound,
    __int128_t a
) {
    std::vector< std::vector<__int128_t> > solutions;
    std::vector< std::vector<__int128_t> > new_solutions;
    /* Strategy:
    Coordinate change d = x + y with inverse y = d - x.
    As d is a factor of  x^3 + y^3, the equation becomes
    (1-z^3)/da = 3x^2 - 3dx + d^2.
    Treat da as an independent variable.  Then d is a divisor of da, and z is 
    a third root of unity modulo da.  The right hand side is bounded below by
    d^2/4, yielding some inequalities. */
    // absolute bound for d
    __int128_t d_bound; 
    __int128_t d, max_shift;
    // Write z = kd + r. Then the bound above yields a lower bound for k.
    __int128_t min_shift = std::max(
        (__int128_t) pow(a/4.0, 1.0/3)-1,
        (__int128_t) 0
    );

    // The upper bounds for d resulting from the bound for the rhs and x, y,
    // respectively.
    d_bound = std::min(
        (__int128_t) pow((4.0*bound*bound*bound + 4.0) / a, 1.0/3),
        2 * bound
    );

    // a vector of all third roots of unity mod d
    std::vector<__int128_t>* roots_of_unity;
    NumberGenerator *gen = new NumberGenerator(d_bound);
    do {
        d = gen->get_value();
        roots_of_unity = gen->roots(); 
        max_shift = bound / d;

        #pragma omp parallel for schedule(static,1)
        for (__int128_t shift = min_shift; shift <= max_shift; ++shift) {
            for (__int128_t root : *roots_of_unity) { 
                // Start with positive z; then d < 0.  Check whether this
                // combination of z and d yields solutions and if so, append
                // them.
                __int128_t z = - (shift+1)*d + root;
                if (yields_solution(d, z, a)){
                    #pragma omp critical
                    {
                        new_solutions = compute_solutions(d, z, a);
                        solutions.insert(
                            solutions.end(),
                            new_solutions.begin(),
                            new_solutions.end()
                        );
                    }
                }

                z = shift*d + root;
                if (z == 1)
                   continue;
                if (yields_solution(-d, z, a)) {
                    #pragma omp critical
                    {
                    new_solutions = compute_solutions(-d, z, a);
                    solutions.insert(
                        solutions.end(),
                        new_solutions.begin(),
                        new_solutions.end()
                    );
                    }
                }
            }
        }
    } while (gen->next());

    delete gen;
    return solutions;
}


/* Compute all points on ax^3 + ay^3 +z^3 = 1 with d = x+y in gen, with
 |x|, |y|, |z| <= bound, and a <= a_max cube free
 (and some solutions not satisfying these constraints).
 Output: vector of vectors of the form [x, y, z, a]. 
 Exclude points with z=1 or x+y = 1. */
std::vector< std::vector<__int128_t> > points_batch (
    std::unique_ptr<NumberGenerator> &gen,
    __int128_t bound,
    __int128_t a_max,
    bool print
) {
    //TODO: this shouldn't need the primelist.
    std::vector< std::vector<__int128_t> > solutions;
    std::vector< std::vector<__int128_t> > new_solutions;

    /* Strategy:
    Coordinate change d = x + y with inverse y = d - x.
    As d is a factor of  x^3 + y^3, the equation becomes
    (1-z^3)/da = 3x^2 - 3dx + d^2.
    Treat da as an independent variable.  Then d is a divisor of da, and z is 
    a third root of unity modulo da. */

    // Bounds for d depending on da, coming from the lower bound of the rhs.
    __int128_t d_lower, d_upper;
    __int128_t a, da, min_shift, max_shift;
    std::list<__int128_t> ds;

    // a vector of all third roots of unity mod d
    std::vector<__int128_t>* roots_of_unity;
    
    do {
        da = gen->get_value();
        d_lower = std::max(
            da / a_max,
            (__int128_t) 2 // d=1 is a trivial case that we exclude;
        );
        d_upper = std::min(
            (__int128_t) (2 * sqrt((1.0+ bound*bound*bound)/da)),
            (__int128_t) da - 1 // d=da corresponds to the excluded case a=1
        );
        gen->divisors(ds, d_lower, d_upper);
        roots_of_unity = gen->roots(); 
        for (__int128_t d : ds) {
            a = da / d;

            // Perfect cubes would result in more accumulating curves and clog
            // up the output.
            if (is_cube(a))
                continue;
                
            max_shift = bound / da;
            min_shift = std::max(
                (__int128_t) (pow(1.0/4 *d*d*d*a, 1.0/3)/da),
                (__int128_t) 0
            );
            
            #pragma omp parallel for schedule(static,1)
            for (__int128_t shift = min_shift; shift <= max_shift; ++shift) {
                for (__int128_t root : *roots_of_unity) { 
                    // Start with negative z; then d > 0.  Check whether this
                    // combination of z and d yields solutions and if so, append
                    // them.
                    __int128_t z = - (shift+1)*da + root;
                    if (yields_solution(d, z, a)){
                        #pragma omp critical 
                        {
                            new_solutions = compute_solutions(d, z, a);
                            for (std::vector<__int128_t> sol : new_solutions){
                                sol.push_back(a);
                                if (print)
                                    print_vector(sol);
                                else
                                    solutions.push_back(sol);
                            }
                        }
                    }

                    z = shift*da + root;
                    if (z == 1)
                       continue;
                    if (yields_solution(-d, z, a)) {
                        #pragma omp critical 
                        {
                            new_solutions = compute_solutions(-d, z, a);
                            for (std::vector<__int128_t> sol : new_solutions){
                                sol.push_back(a);
                                if (print)
                                    print_vector(sol);
                                else
                                    solutions.push_back(sol);
                            }
                        }
                    }
                }
            }
        }
    } while (gen->next());

    return solutions;
}

// Compute all points of height <= bound on surfaces with a <= a_max.
std::vector< std::vector<__int128_t> > points_batch (
    __int128_t bound,
    __int128_t a_max,
    bool print
) {
    // Upper bound for da resulting from the lower bound of the rhs and the
    // upper bounds for x, y, a, respectively.
    __int128_t da_bound = get_da_bound(bound, a_max);
    std::unique_ptr<NumberGenerator> gen = std::make_unique<NumberGenerator>(da_bound);
    
    return points_batch(gen, bound, a_max, print);
}


// Compute the set of solutions where the largest prime factor of da is in 
// [largest_prime_lower, largest_prime_upper).  If largest_prime_lower is <=1,
// da = 1 is included.
std::vector< std::vector<__int128_t> > points_batch_partial (
    __int128_t bound,
    __int128_t a_max,
    __int128_t largest_prime_lower,
    __int128_t largest_prime_upper,
    bool print
) {
    std::unique_ptr<NumberGenerator> gen;

    // Upper bound for da resulting from the lower bound of the rhs and the
    // upper bounds for x, y, a, respectively.
    __int128_t da_bound = get_da_bound(bound, a_max);
    if (largest_prime_lower >= largest_prime_upper) 
        return std::vector< std::vector<__int128_t> >();

    gen = std::make_unique<NumberGenerator>(
        da_bound, largest_prime_lower, largest_prime_upper
    );

    return points_batch(gen, bound, a_max, print);

}

