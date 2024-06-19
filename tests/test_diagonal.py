"""Unit tests for the treatment of the diagonal surfaces in C++. """
import unittest
from math import sqrt
from sympy import divisors, primerange
from cypari2 import Pari

from cubicpts import diagonal
import random


class TestCountDiagonal(unittest.TestCase):
    """Tests for the counting functions in diagonal."""

    def test_small_points_py(self):
        set_brute = filter_points(points_brute(100), 100)
        set_alg = filter_points(diagonal.points(100), 100)
        self.assertEqual(set_brute, set_alg)

    def test_with_pari(self):
        set_pari = filter_points(points_pari(10**4), 10**4)
        set_alg =  filter_points(diagonal.points(10**4), 10**4)
        self.assertEqual(set_pari, set_alg)

    def test_with_sympy(self):
        set_sympy = filter_points(points_py(10**4), 10**4)
        set_alg =  filter_points(diagonal.points(10**4), 10**4)
        self.assertEqual(set_sympy, set_alg)

    def test_more_as(self):
        set_7 = filter_points(diagonal.points(10**4, 7), 10**4, 7)
        set_sympy = filter_points(points_py(10**4, 7), 10**4, 7)
        self.assertEqual(set_7, set_sympy)


class TestBatchCounts(unittest.TestCase):
    """Tests for the batch counting functions in diagonal."""

    def test_small_amax(self):
        set_all = diagonal.points_batch(10**3, 4)
        set_2 = filter_points([P for (P, a) in set_all if a == 2], 10**3, 2)
        set_4 = filter_points([P for (P, a) in set_all if a == 4], 10**3, 4)
        set_2_direct = filter_points(diagonal.points(10**4, 2), 10**3, 2)
        set_4_direct = filter_points(diagonal.points(10**4, 4), 10**3, 4)
        #self.assertEqual(len(set_2), len(set_2_direct))
        self.assertEqual(set_2, set_2_direct)
        self.assertEqual(set_4, set_4_direct)

    def test_larger_amax(self):
        bound = 10**3
        set_all = diagonal.points_batch(bound, 100)
        set_7  = filter_points([P for (P, a) in set_all if a == 7],  bound, 7)
        set_97 = filter_points([P for (P, a) in set_all if a == 97], bound, 97)
        set_7_direct  = filter_points(diagonal.points(bound, 7),  bound, 7)
        set_97_direct = filter_points(diagonal.points(bound, 97), bound, 97)
        #self.assertEqual(len(set_2), len(set_2_direct))
        self.assertEqual(set_7, set_7_direct)
        self.assertEqual(set_97, set_97_direct)


class TestParallelism(unittest.TestCase):

    def test_small_a_max(self):
        bound = 10**3
        a_max = 100
        set_direct = set(diagonal.points_batch(bound, 100))
        set_parallel = set()
        for i in range(0, 100):
            new_solutions = diagonal.points_batch_partial(bound, a_max, i, 100)
            for P in new_solutions:
                set_parallel.add(((P[0], P[1], P[2]), P[3]))
        self.assertEqual(set_direct, set_parallel)
        

class TestSanity(unittest.TestCase):
    """Some sanity tests for helper functions."""

    def test_sanity(self):
        self.assertTrue(is_square(16))
        self.assertFalse(is_square(17))
        self.assertEqual(sgn(4),1)
        self.assertEqual(sgn(-5),-1)
        set_brute = filter_points(points_brute(100), 100)
        set_pari = filter_points(points_pari(100), 100)
        set_py =  filter_points(points_py(100), 100)
        self.assertEqual(set_brute, set_pari)
        self.assertEqual(set_pari, set_py)
        self.assertEqual(set_brute, set_py)

    def test_sanity2(self):
        set_brute = filter_points(points_brute(100, 7), 100, 7)
        set_pari = filter_points(points_pari(100, 7), 100, 7)
        set_py =  filter_points(points_py(100, 7), 100, 7)
        self.assertEqual(set_brute, set_pari)
        self.assertEqual(set_pari, set_py)
        self.assertEqual(set_brute, set_py)


class TestNumberGenerator(unittest.TestCase):

    def test_constructor(self):
        gen = diagonal.PyNumberGenerator(1000)

    def test_set_small(self):
        """Check whether the generator correctly yields the first 20 numbers"""
        gen = diagonal.PyNumberGenerator(20)
        
        numbers = []
        while True:
            numbers.append(int(gen.value()))
            if not gen.next():
                break
        self.assertEqual(len(numbers), 20)
        self.assertEqual(set(numbers), set(range(1,21)))
        del gen

    def test_set(self):
        """Check whether the generator correctly yields the first 1000 numbers"""
        gen = diagonal.PyNumberGenerator(1000)
        
        numbers = []
        while True:
            numbers.append(int(gen.value()))
            if not gen.next():
                break
        self.assertEqual(len(numbers), 1000)
        self.assertEqual(set(numbers), set(range(1,1001)))
        del gen

    def test_divisors(self):
        gen = diagonal.PyNumberGenerator(10)
        for _ in range(0,5):
            gen.next()
        self.assertEqual(gen.value(), 6)
        divs = {int(d) for d in gen.divisors(0,-1)}
        self.assertEqual(divs, {1,2,3,6})
        divs = {int(d) for d in gen.divisors(2,3)}
        self.assertEqual(divs, {2,3})
    
    def test_roots_of_unity(self):
        gen = diagonal.PyNumberGenerator(10)
        for _ in range(0,5):
            gen.next()
        self.assertEqual(gen.value(), 6)
        self.assertCountEqual(gen.roots(), [1])
        
        for _ in range (0,4):
            gen.next()
        self.assertEqual(gen.value(), 7)
        self.assertCountEqual(gen.roots(), [1, 2, 4])

class TestModularArithmetic(unittest.TestCase):

    def test_pow(self):
        self.assertEqual(diagonal._pow_mod( 2, 3, 10), 8)
        self.assertEqual(diagonal._pow_mod( 2, 3,  5), 3)
        self.assertEqual(diagonal._pow_mod(-3, 3,  5), 3)
        self.assertEqual(diagonal._pow_mod( 3, 3,  9), 0)


class TestNTh(unittest.TestCase):

    def test_fourth_power(self):
        self.assertTrue(diagonal._is_fourth_power(16))
        self.assertFalse(diagonal._is_fourth_power(4))

    def test_is_square(self):
        self.assertTrue(diagonal.is_square(16))
        for _ in range(0,15):
            n = random.randint(1, 2**63)
            self.assertTrue(diagonal.is_square(n**2))

class TestIO(unittest.TestCase):

    def test_string_conversion(self):
        for i in [0, 1, -1, 10**15, 10**19, 2**100, -2**80]:
            self.assertEqual(int(diagonal._to_string(i)), i)

def filter_points(points, bound, a=2):
    """Return the set of nontrivial solutions of bounded height."""
    result = set()
    for (x, y, z) in points:
        if (x!=0 and y!=0 and z not in (-1,0,1) and x+y not in (-1, 0, 1)
                and not diagonal.on_accumulating_curve(x, y, z, a)
                and max(abs(x), abs(y), abs(z)) <= bound):
            result.add((int(x), int(y), int(z)))
    return result


def is_square(n):
    """Checks whether a number is a square."""
    if n < 0:
        return False
    candidate = int(sqrt(n))
    if candidate * candidate == n:
        return True
    return False


def sgn(n):
    """Returns the sign of a number."""
    if n > 0:
        return 1
    if n == 0:
        return 0
    return -1


def points_brute(bound, a = 2):
    """
    Compute solutions to ax^3 + ay^3 + z^3 = 1  of height at most bound
    by brute force.
    """
    sols = []
    for x in range(-bound, bound+1):
        for y in range(-bound, bound+1):
            for z in range(-bound, bound+1):
                if z == 1:
                    continue
                if a*x*x*x + a*y*y*y + z*z*z == 1:
                    sols.append((x, y, z))
    return sols


def points_pari(bound, a = 2):
    """
    Compute solutions to ax^3 + ay^3 + z^3 = 1  of height at most bound, using 
    the factoring algorithm in PARI.
    """
    # Run over z.  Write (1-z^3)/a = (x+y)(x^2 - xy + y^2) to see that d=|x+y|
    # must be a divisor of the left hand side.
    pari = Pari()
    solutions = []
    ds = []
    for z in range(-bound, bound + 1):
        if z == 1:
            continue
        lhs_0 = 1-z*z*z
        if lhs_0 % a != 0:
            continue
        lhs = lhs_0 // a
        assert lhs != 0
        ds = pari.divisors(lhs)

        for d in ds:
            # Solve the equation on the right (with y = d - x)
            disc_num = 4 * abs(lhs//d) -d*d
            if disc_num % 3 != 0:
                continue
            disc = disc_num // 3
            if not is_square(disc):
                continue
            sqrt_disc = int(sqrt(disc))
            x_num = sgn(lhs) * (d+sqrt_disc)
            if x_num % 2 != 0:
                continue
            x = x_num//2
            y = sgn(lhs)*(d-sqrt_disc)//2
            solutions.append((x, y, z))
            if sqrt_disc != 0:
                solutions.append((y, x, z))

    return solutions


def points_py(bound, a = 2):
    """
    Compute solutions to ax^3 + ay^3 + z^3 = 1  of height at most bound
    by brute force.
    """
    # Run over z.  Write (1-z^3)/a = (x+y)(x^2 - xy + y^2) to see that d=|x+y|
    # must be a divisor of the left hand side.
    solutions = []
    ds = []
    for z in range(-bound, bound + 1):
        if z == 1:
            continue
        lhs_0 = 1-z*z*z
        if lhs_0 % a != 0:
            continue
        lhs = lhs_0 // a
        ds = divisors(lhs, generator = True)
        for d in ds:
            # Solve the equation on the right (with y = d - x)
            disc_num = 4 * abs(lhs//d) -d*d
            if disc_num % 3 != 0:
                continue
            disc = disc_num // 3
            if not is_square(disc):
                continue
            sqrt_disc = int(sqrt(disc))
            x_num = sgn(lhs) * (d+sqrt_disc)
            if x_num % 2 != 0:
                continue
            x = x_num//2
            y = sgn(lhs)*(d-sqrt_disc)//2
            solutions.append((x, y, z))
            if sqrt_disc != 0:
                solutions.append((y, x, z))

    return solutions
