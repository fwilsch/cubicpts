"""Unit tests for the treatment of the Baragar-Umeda surfaces. """
import unittest

from cubicpts import tr_linear as lin

class TestCountEx4(unittest.TestCase):
    """Tests for the counting functions in ex4."""

    def test_small_points(self):
        set_brute = set(lin.points_brute(15))
        set_alg = set(lin.points(15))
        self.assertEqual(set_brute, set_alg)

    def test_small_points_2(self):
        set_brute = set(lin.points_var(1000))
        set_alg = set(lin.points(1000, 2, 0, 1, -1))
        self.assertEqual(set_brute, set_alg)


class TestCountHarpaz(unittest.TestCase):
    """Tests for the family of surfaces studied by harpaz"""
    
    def test_small_points(self):
        bound = 30
        a = 2
        set_brute = set(points_brute(bound, a))
        set_alg = set(lin.points_var(bound, a))
        self.assertEqual(set_brute, set_alg)
        a = 5 
        set_brute = set(points_brute(bound, a))
        set_alg = set(lin.points_var(bound, a))
        self.assertEqual(set_brute, set_alg)

def points_brute(bound, a):
    solutions = []
    for x in range(-bound, bound):
        for y in range(-bound, bound):
            lhs = x*x - a*y*y
            rhs = y-1
            for z in range(-bound, bound):
                if z == 0:
                    continue
                if lhs*z - rhs == 0:
                    solutions.append((x, y, z))
    return solutions


if __name__ == "__main__":
    unittest.main()
