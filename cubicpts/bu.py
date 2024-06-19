"""
Compute solutions of ax^2 + by^2 +cz^2 = dxyz + 1.
"""

from math import gcd, pi, prod, sqrt, log

from sympy import legendre_symbol, jacobi_symbol, primerange, primefactors
from sympy.ntheory.factor_ import core
#import bigfloat

from cubicpts import tools


BU_PARAMS = [
    (1, 5, 5, 5), (1, 3, 6, 6), (2, 7, 14, 14),
    (2, 2, 3, 6), (6, 10, 15, 30), (1, 2, 2, 2)
]
"""The parameters treated in Baragar-Umeda."""

PRIMITIVE_SOLUTIONS = {
    (1, 5, 5, 5): [(4, 2, 1), (4, 1, 2)],
    (1, 3, 6, 6): [(2, 1, 1)],
    #(3, 4, 6, 12): [(1, 1, 1)],
    (2, 7, 14, 14): [(2, 1, 1)],
    (2, 2, 3, 6): [(1, 1, 1)],
    (6, 10, 15, 30): [(1, 1, 1)],
    (1, 2, 2, 2): [(3, 2, 2)],
    #(1, 8, 8, 8): [(3, 1, 1)]
}
"""The list of non-trivial primitive solutions for each choice of parameter."""

#with bigfloat.quadruple_precision:
#    PRIMITIVE_SOLUTIONS_F = {
#        (1, 5, 5, 5): [
#            (bigfloat.BigFloat(4.0, context=bigfloat.quadruple_precision),
#             bigfloat.BigFloat(2.0, context=bigfloat.quadruple_precision),
#             bigfloat.BigFloat(1.0, context=bigfloat.quadruple_precision)
#            ),
#            (
#             bigfloat.BigFloat(4.0, context=bigfloat.quadruple_precision),
#             bigfloat.BigFloat(1.0, context=bigfloat.quadruple_precision),
#             bigfloat.BigFloat(2.0, context=bigfloat.quadruple_precision)
#            )
#        ],
#        (1, 3, 6, 6): [(2.0, 1.0, 1.0)],
#        #(3, 4, 6, 12): [(1.0, 1.0, 1.0)],
#        (2, 7, 14, 14): [(2.0, 1.0, 1.0)],
#        (2, 2, 3, 6): [(1.0, 1.0, 1.0)],
#        (6, 10, 15, 30): [(1.0, 1.0, 1.0)],
#        (1, 2, 2, 2): [(3.0, 2.0, 2.0)],
#        #(1, 8, 8, 8): [(3.0, 1.0, 1.0)]
#    }

L_PRINCIPAL = {
    (1, 5, 5, 5): 4 * pi**4 / (9*3*15),
    (1, 3, 6, 6): pi**4 / (4*4* 3*sqrt(3) * sqrt(6)),
    #(3, 4, 6, 12): 0,
    (2, 7, 14, 14): 2 * pi**4 / (3*3*21),
    (2, 2, 3, 6): None,
    (6, 10, 15, 30): None,
    (1, 2, 2, 2): None
    #(1, 8, 8, 8): 0,
}
"""The products of the principal values of the L-functions used for
convergence."""

BU_CONSTANT = {
    (1, 5, 5, 5): 3.92062681166,
    (1, 3, 6, 6):  2.22381295435,
    #(3, 4, 6, 12): 0,
    (2, 7, 14, 14): 1.85092947320,
    (2, 2, 3, 6): 3.04230700308,
    (6, 10, 15, 30): 1.86988733010,
    (1, 2, 2, 2): 3.69061353513
    #(1, 8, 8, 8): 0
}
"""The leading constant determined in Baragar-Umeda."""


def points_mod_p(mod, a=1, b=5, c=5, d=5):
    """Compute points of the equation modulo mod."""
    results = []
    for x in range(0, p):
        for y in range(0, p):
            f_xy = (a*x*x + b*y*y - 1) % p # part of f that depends on x and y
            dxy = (d*x*y) % p
            for z in range(0, p):
                if (f_xy + c*z*z - dxy*z) % p == 0:
                    results.append((x, y, z))

    return results


def points_brute(bound, a=1, b=5, c=5, d=5):
    """Compute points of the equation with 1<= x,y,z <= bound by brute force."""
    results = []
    for x in range(1, bound+1):
        for y in range(1, bound+1):
            f_xy = a*x*x + b*y*y - 1 # part of f that depends on x and y
            dxy = d*x*y
            for z in range(1, bound+1):
                if f_xy + c*z*z - dxy*z == 0:
                    results.append((x, y, z))

    return results


def euler_factor_old(p):
    """Return sigma'_p(0) as in (6.2)"""
    chi_3 = 0
    chi_5 = 0
    if p not in (2, 3):
        chi_3 = legendre_symbol(-3, p)
    if p not in (2, 5):
        chi_5 = legendre_symbol(5, p)
    chi_15 = chi_3 * chi_5

    conv_factor = (p-chi_3)*(p-chi_15)/(p*p)
    conv_factor = conv_factor * conv_factor

    # For the bad primes, compute solutions mod p (resp. 2^3) and count them,
    # removing the trivial solutions and the odd solutions obstructed by BM.
    # For the remaining primes, use the closed formula accounting for failures
    # of strong approximation explained by the group action.  Multiply
    # everything with L-functions that make the Euler product absulutely
    # convergent and converge much faster.  (In the end: have to multiply with
    # the inverses of the principal values of these functions.)
    if p == 2:
        list_8 = points_mod_p(8)
        list_8_n = [(x, y, z) for (x, y, z) in list_8 if x%2 == 0]
        v_8 = len(list_8_n)
        sigma_0 = v_8 / (8*8)
    elif p == 3:
        list_3 = points_mod_p(3)
        list_3_n = [(x, y, z) for (x, y, z) in list_3
            if (x, y, z) not in ((1,0,0), (2,0,0))
            and x!= 0]
        v_3 = len(list_3_n)
        sigma_0 = v_3/9 # This is halved to accomodate for BMO involving 3 and 5
    elif p == 5:
        list_5 = points_mod_p(5)
        list_5_n = [P for P in list_5 if P != (1,0,0) and P !=(4,0,0)]
        v_5 = len(list_5_n)
        sigma_0 = v_5/25 # Use all solutions (except (+/-1,0,0)) here
    else:
        sigma_0 = 1 + 2*(chi_3+chi_15)/p - (3 + 2*chi_5)/(p*p)

    return conv_factor * sigma_0


def dirichlet_character(a):
    """Return a list containing the values of the Dirichlet character (a/p)."""
    chi = []
    if a == 1:
        return [1,1,1,1]
    for p in range(0,abs(4*a)):
        if gcd(p,2*a) != 1:
            chi.append(0)
        else:
            chi.append(jacobi_symbol(a, p))
    return chi


def euler_product(p_max, a=1, b=5, c=5, d=5):
    """
    Return the Euler product part of the leading constant.

    That is, the product over all primes of
        - convergence factors L_p(1,chi_1) \\cdots L_p(1,chi_4), where
          chi_1 = (k/p), chi_2 = (ka/p), ..., chi_4 = (kc/p) up to clearing
          squares,
        - the number of points in the orbit of the primitive solutions in
          \\ZZ/m\\ZZ, where m=2^3(d^2-4abc), divided by m^2, and
        - the product over the number of solutions mod p for all primes not
          dividing m, after removing solutions (x, 0, 0), (0, y, 0), (0, 0, z)
          that can never lift to non-trivial solutions, divided by p^2.
    """
    k = d*d - 4*a*b*c
    # We need squarefree representations of k, ka, kb, kc.  Note that k is
    # always negative, while a, b, c are always positive
    k_red  = - core(abs(k), 2)
    ka_red = - core(abs(k*a), 2)
    kb_red = - core(abs(k*b), 2)
    kc_red = - core(abs(k*c), 2)

    # Initialize legendre symbols by computing all values once, set up functions
    # that return the values by reducing the argument modulo 4n and looking up
    # the result.
    chi_k_list  = dirichlet_character(k_red)
    chi_a_list  = dirichlet_character(a)
    chi_b_list  = dirichlet_character(b)
    chi_c_list  = dirichlet_character(c)
    chi_ka_list = dirichlet_character(ka_red)
    chi_kb_list = dirichlet_character(kb_red)
    chi_kc_list = dirichlet_character(kc_red)

    chi_k  = lambda p : chi_k_list[p%(-4*k_red)]
    chi_a  = lambda p : chi_a_list[p%(4*a)]
    chi_b  = lambda p : chi_b_list[p%(4*b)]
    chi_c  = lambda p : chi_c_list[p%(4*c)]
    chi_ka = lambda p : chi_ka_list[p%(4*ka_red)]
    chi_kb = lambda p : chi_kb_list[p%(4*kb_red)]
    chi_kc = lambda p : chi_kc_list[p%(4*kc_red)]

    sigma = lambda p: (
        1
        + (chi_k(p) + chi_ka(p) + chi_kb(p) + chi_kc(p)) / p
        - (2 + chi_a(p) + chi_b(p) + chi_c(p))/(p*p)
    )
    convergence_factor = lambda p: (
        (1 - chi_k(p)/p) * (1 - chi_ka(p)/p)
        * (1 - chi_kb(p)/p) * (1 - chi_kc(p)/p)
    )

    mod = 8 * (4*a*b*c - d*d)
    bad_primes = set(primefactors(mod))
    solutions_mod = len(orbit_mod(mod, a, b, c, d))
    bad_conv_factors = prod((convergence_factor(p) for p in bad_primes ))
    bad_sigma_prod = bad_conv_factors * solutions_mod / (mod*mod)

    sigmas_prime = (
        sigma(p) * convergence_factor(p)
        for p in primerange(1, p_max)
        if p not in bad_primes
    )

    sigma_prod = prod(sigmas_prime)

    return bad_sigma_prod * sigma_prod

def get_constant(p_max, a=1, b=5, c=5, d=5):
    """
    Compute the conjectured leading constant after accounting for BMO.

    Do so by computing the Euler product for primes < p_max.
    """
    euler = euler_product(p_max, a, b, c, d)
    # We have to multiply with the real density and the principal value of the
    # L-function used as a convergence factor.
    return 6/d * L_PRINCIPAL[(a,b,c,d)] * euler


def get_constant_old(p_bound):
    """
    Compute the conjectured leading constant after accounting for BMO.

    Only for the specific surface x^2 + 5y^2 + 5z^2 = 5xyz + 1.
    """
    sigma = 1
    bad_sigma = prod((euler_factor_old(p) for p in [2,3,5]))
    sigma_values = (euler_factor_old(p) for p in primerange(6,p_bound))
    sigma = prod(sigma_values)
    #L_principal = 4 * pi**4 / (9*3*15)
    return 6/5 * L_PRINCIPAL[(1,5,5,5)] * sigma * bad_sigma


def points(bound, a=1, b=5, c=5, d=5):
    """
    Compute positive solutions using the Vieta involution.

    args:
        bound: The bound up to which to search.
        a,b,c,d: The parameters. They must be in ascending order and in the
                 Baragar-Umeda list
    """
    # solutions: Tuples (P, i), where P=(x,y,z) is a point and i is the
    # involution that was used to obtain it.
    #
    # Initialize with primitive solutions.
    if (a, b, c, d) in PRIMITIVE_SOLUTIONS:
        solutions = [(P, -1) for P in PRIMITIVE_SOLUTIONS[(a, b, c, d)]]
    else:
        raise ValueError("Parameters must be part of the BU list.")

    for ((x, y, z), i) in solutions:
        x_new = d//a *y*z - x
        y_new = d//b *x*z - y
        z_new = d//c *x*y - z

        # As soon as one variable is too large, the solution cannot become
        # smaller again.  Take care not to add copies of the primitive
        # solutions, which might be invariant under one involution.
        if i != 0 and x_new <= bound and x_new != x:
            solutions.append(((x_new, y, z), 0))
        if i != 1 and y_new <= bound and y_new != y:
            solutions.append(((x, y_new, z), 1))
        if i != 2 and z_new <= bound and z_new != z:
            solutions.append(((x, y, z_new), 2))

    return [P for (P,i) in solutions]


def experimental_constant(bound, a, b, c, d):
    """
    Compute an approximate leading constant by counting points of height
    <= bound.
    """
    solutions = points(bound, a, b, c, d)
    # There are C (log B)^2 points.  Multiply by 4 to account for the fact that
    # points() only yields positive solutions.
    return 4 * len(solutions) / (log(bound)**2)


def orbit_mod(mod, a=1, b=5, c=5, d=5):
    """
    Compute the orbits of the primitive solutions in the points modulo mod.
    """
    if (a, b, c, d) in PRIMITIVE_SOLUTIONS:
        orbit = [tools.reducepoint(P, mod)
                 for P in PRIMITIVE_SOLUTIONS[(a, b, c, d)]]
    else:
        raise ValueError("Parameters must be part of the Baragar-Umeda list.")

    # Use an additional set structure for faster lookup, and keep a list to
    # iterate on while expanding it.
    orbit_set = set(orbit)

    for (x, y, z) in orbit:
        #Compute the points arising from (x,y,z) by the three involutions
        x_new = (d//a * y*z - x) % mod
        y_new = (d//b * x*z - y) % mod
        z_new = (d//c * x*y - z) % mod
        vieta_points = [(x_new, y, z), (x, y_new, z), (x, y, z_new)]
        for point in vieta_points:
            # Add any solution that we haven't found yet to both the set and
            # the list.
            if point not in orbit_set:
                orbit.append(point)
                orbit_set.add(point)

    return orbit


def test_mod(p, a=1, b=5, c=5, d=5):
    """
    Compute the complement of the reductions of integral points in the Z/pZ-points.
    """
    points_mod = set(points_mod_p(p, a, b, c, d))
    orbits = [tools.reducepoint(point, p)
              for point in PRIMITIVE_SOLUTIONS[(a,b,c,d)]]

    for (x, y, z) in orbits:
        if (x, y, z) in points_mod:
            points_mod.remove((x, y, z))

    for (x, y, z) in orbits:
        x_new = (d//a * y*z - x) % p
        y_new = (d//b * x*z - y) % p
        z_new = (d//c * x*y - z) % p
        vieta_points = [(x_new, y, z), (x, y_new, z), (x, y, z_new)]
        for point in vieta_points:
            if point in points_mod:
                points_mod.remove(point)
                orbits.append(point)

    return points_mod

