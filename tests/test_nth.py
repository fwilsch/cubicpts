"""Unit tests for the treatment of the Baragar-Umeda surfaces. """
import unittest

import pyximport

pyximport.install(setup_args = {"script_args" : ["--force"]},
                  language_level=3)

from cubicpts import nth

class TestCRTList(unittest.TestCase):
    """Test the CRT function."""
    def test_crt(self):
        assert True

class TestSquareTest(unittest.TestCase):
    """Tests for the square testing."""

    def test_small_squares(self):
        assert nth.is_square(0)
        assert nth.is_square(1)
        assert nth.is_square(4)
        assert nth.is_square(7**2)
        assert nth.is_square(12**2)
        assert nth.is_square(1363**2)

    def test_small_nonsquares(self):
        assert not nth.is_square(2)
        assert not nth.is_square(3)
        assert not nth.is_square(5)
        assert not nth.is_square(6)
        assert not nth.is_square(7)
        assert not nth.is_square(8)
        assert not nth.is_square(1363**2+1)

    def test_small_numbers(self):
        square = [False]*101
        for k in range(0, 11):
            square[k*k] = True
        for k in range(0, 101):
            assert (nth.is_square(k) == square[k])

    def test_large_squares(self):
        assert nth.is_square(10**20)
        assert nth.is_square(3**40)
        assert nth.is_square(5**40)
        assert nth.is_square((7*2)**20)
        assert nth.is_square(123456789012345**2)

    def test_large_numbers(self):
        for k in range(10**10-100, 10**10+100):
            assert nth.is_square(k*k)
            assert not nth.is_square(k**2+1)
            assert not nth.is_square(k**2-1)
        for k in range(10**12-100, 10**12+100):
            assert nth.is_square(k*k)
            assert not nth.is_square(k**2+1)
            assert not nth.is_square(k**2-1)
        for k in range(10**15-100, 10**15+100):
            assert nth.is_square(k*k)
            assert not nth.is_square(k**2+1)
            assert not nth.is_square(k**2-1)


if __name__ == "__main__":
    unittest.main()
