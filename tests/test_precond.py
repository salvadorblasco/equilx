import unittest

import numpy as np

from equilx.data_structure import EquilibriumData
from equilx.solver.precond import _pick_exponent


class TestPrecond(unittest.TestCase):
    def test_pick_exponent(self):
        ...

    def setUp(self):
        self.equil1 = self.__fixture1()

    def __fixture1(self):
        labels = "L Cu Im H".split()
        logbeta = np.fromiter(map(float, "10.66 20.07 28.31 35.67 42.57 48.44 7.09 4.2 7.73 10.6 13   44.36 40.54 34.41 26.97 17.21 25.18 17.97 50.43 45.82 38.97 36.3 29.79 22.28 -13.77".split()), dtype=float)
        stoich = np.fromiter(map(float, "1   0   0   1 1   0   0   2 1   0   0   3 1   0   0   4 1   0   0   5 1   0   0   6 0   0   1   1 0   1   1   0 0   1   2   0 0   1   3   0 0   1   4   0 1   1   0   4 1   1   0   3 1   1   0   2 1   1   0   1 1   1   0   0 1   2   0   0 1   2   0  -1 1   1   1   4 1   1   1   3 1   1   1   2 1   1   1   1 1   2   1   0 1   2   1  -1 0   0   0  -1".split()), dtype=int)
        # <initial unit="mol/L">0.001 0.001 0.001 2.0</initial>
        # <final unit="mol/L">0.001 0.001 0.001 11.0</final>
        # <pX>False False False True</pX>
        # <xref p="yes">3</xref
