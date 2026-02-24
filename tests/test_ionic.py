#!/usr/bin/env -S uv run --script

import unittest
import numpy as np
import numpy.testing as npt

from equilx.data_structure import DebyeHuckelParameters
from equilx.solver.ionic import calc_ionic_strength, beta_increment, debye_huckel


class DummyEq:
    ...


class TestIonic(unittest.TestCase):
    def test_calc_is(self):
        equil = DummyEq()
        equil.chargesx = np.array([-1, 1, 0, -1], dtype=int)
        conc = np.ones((10,4), dtype=float)
        ionic = calc_ionic_strength(conc, equil)
        npt.assert_equal(np.full(10, 1.5), ionic)

    def test_debye_huckel(self):
        equil = DummyEq()
        equil.zstar = np.array([2,0], dtype=int)
        equil.pstar = np.array([1,-2], dtype=int)
        dbp = DebyeHuckelParameters(a=1.0,
                                    b=1.0,
                                    c0=0.0, c1=0.0,
                                    d0=0.0, d1=0.0,
                                    e0=0.0, e1=0.0)
        with self.subTest("float 0.0"):
            ionic = np.array([0.0])
            calc = debye_huckel(ionic, equil, dbp)
            npt.assert_equal(calc, np.zeros_like(calc))

        with self.subTest("float 1.0"):
            ionic = np.array([1.0])
            calc = debye_huckel(ionic, equil, dbp)
            npt.assert_equal(calc, np.array([-1.0,0.0]))

        with self.subTest("array"):
            ionic = np.linspace(0.0, 1.5, 10)
            calc = debye_huckel(ionic, equil, dbp)
            true = np.array([[ 0.        , 0.],
                             [-0.5797959 , 0.],
                             [-0.73205081, 0.],
                             [-0.82842712, 0.],
                             [-0.89897949, 0.],
                             [-0.95445115, 0.],
                             [-1.        , 0.],
                             [-1.0385186 , 0.],
                             [-1.07179677, 0.],
                             [-1.10102051, 0.]])
            npt.assert_allclose(true, calc)

if __name__ == '__main__':
    unittest.main()
