#!/usr/bin/env -S uv run --script

import unittest

import numpy as np
import numpy.testing as npt

from equilx.data_structure import EquilibriumData
from equilx.solver.utils import make_distribution
from equilx.solver.precond import _pick_exponent


class TestPrecond(unittest.TestCase):
    def test_pick_exponent(self):
        retv = _pick_exponent(self.distri1)
        npt.assert_equal(1/np.array([1.0,2.0,4.0,6.0]), retv)

    def setUp(self):
        self.distri1 = self.__fixture1()

    def __fixture1(self):
        labels = "L Cu Im H".split()
        n_comps = len(labels)
        logbeta = np.fromiter(map(float, "10.66 20.07 28.31 35.67 42.57 48.44 7.09 4.2 7.73 10.6 13   44.36 40.54 34.41 26.97 17.21 25.18 17.97 50.43 45.82 38.97 36.3 29.79 22.28 -13.77".split()), dtype=float)
        n_equil = len(logbeta)
        stoich = np.fromiter(map(float, "1 0 0 1 1 0 0 2 1 0 0 3 1 0 0 4 1 0 0 5 1 0 0 6 0 0 1 1 0 1 1 0 0 1 2 0 0 1 3 0 0 1 4 0 1 1 0 4 1 1 0 3 1 1 0 2 1 1 0 1 1 1 0 0 1 2 0 0 1 2 0 -1 1 1 1 4 1 1 1 3 1 1 1 2 1 1 1 1 1 2 1 0 1 2 1 -1 0 0 0 -1".split()), dtype=int).reshape(n_equil, n_comps)
        analc = make_distribution(np.array([0.001, 0.001, 0.001, 2.0]),
                                  np.array([0.001, 0.001, 0.001, 11.0]),
                                  [False, False, False, True], 10)
        return EquilibriumData(analc, logbeta, stoich)


if __name__ == '__main__':
    unittest.main()
