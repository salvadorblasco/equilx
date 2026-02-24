#!/usr/bin/env -S uv run --script

import unittest
import numpy as np
import numpy.testing as npt

from equilx.solver.ionic import calc_ionic_strength, beta_increment


class DummyEq:
    ...


class TestIonic(unittest.TestCase):
    def test_calc_is(self):
        equil = DummyEq()
        equil.chargesx = np.array([-1, 1, 0, -1], dtype=int)
        conc = np.ones((10,4), dtype=float)
        ionic = calc_ionic_strength(conc, equil)
        npt.assert_equal(np.full(10, 1.5), ionic)

if __name__ == '__main__':
    unittest.main()
