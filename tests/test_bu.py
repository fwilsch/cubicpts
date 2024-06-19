"""Unit tests for the treatment of the Baragar-Umeda surfaces. """
import unittest

from cubicpts import bu

class TestBUCount(unittest.TestCase):
    """Tests for the counting functions."""
    
#    def test_float_small(self):
#        """
#        Check whether the float version and int version yield the same result
#        for small n.
#        """
#        n_1 = len(bu.points(1000))
#        n_2 = len(bu.points_f(1000))
#        self.assertEqual(n_1, n_2)

    def test_small_points(self):
        set_brute = set(bu.points_brute(15))
        set_vieta = set(bu.points(15))
        self.assertEqual(set_brute, set_vieta)

    def test_small_points_all(self):
        for param in bu.BU_PARAMS:
            set_brute = set(bu.points_brute(15, *param))
            set_vieta = set(bu.points(15, *param))
            self.assertEqual(set_brute, set_vieta)


if __name__ == "__main__":
    unittest.main()
