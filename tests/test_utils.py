#!/usr/bin/env -S uv run --script

import unittest

import numpy as np
import numpy.testing as npt

import equilx.solver.utils


class TestUtilsDistri(unittest.TestCase):
    def test_pick_exponent(self):
        res = equilx.solver.utils.make_distribution(self.initc, self.finac, self.pX, 8)
        npt.assert_allclose(self.true_result, res)

    def setUp(self):
        self.initc = np.array([0.01, 2, 0.01, -2.0], dtype=float)
        self.finac = np.array([0.29, 5.5, 0.08, 5.0], dtype=float)
        self.pX = [False, True, False, True]
        self.true_result = np.array([[1.0E-2, 1.00000000E+2, 1.0E-2, 1.0E-2],
                                     [5.0E-2, 3.16227766E+2, 2.0E-2, 1.0E-1],
                                     [9.0E-2, 1.00000000E+3, 3.0E-2, 1.0E+0],
                                     [1.3E-1, 3.16227766E+3, 4.0E-2, 1.0E+1],
                                     [1.7E-1, 1.00000000E+4, 5.0E-2, 1.0E+2],
                                     [2.1E-1, 3.16227766E+4, 6.0E-2, 1.0E+3],
                                     [2.5E-1, 1.00000000E+5, 7.0E-2, 1.0E+4],
                                     [2.9E-1, 3.16227766E+5, 8.0E-2, 1.0E+5]])

if __name__ == '__main__':
    unittest.main()
