#!/usr/bin/env -S uv run --script

import unittest

import numpy as np
import numpy.testing as npt

from equilx.data_structure import EquilibriumData
from equilx.solver.utils import mass_action, species_concentration
from equilx.solver.precond import _pick_exponent, _newx, _update_theta, precondition, _calc_ratio
from equilx.solver.solver import guess_starting_point

from fixtures import synthetic_simple1_wsld, synthetic_simple1_wosld
from fixtures2 import all_fixtures, large_fixture1


class TestPrecond(unittest.TestCase):
    def test_pick_exponent(self):
        retv = _pick_exponent(large_fixture1())
        npt.assert_equal(1/np.array([1.0,2.0,4.0,6.0]), retv)

    def test_newx(self):
        n_points:int = 10
        n_compon:int = 3
        x0 = np.random.rand(n_points, n_compon)
        ratio = np.random.rand(n_points, n_compon)
        exp = np.array(n_compon*[1.0], dtype=float)
        with self.subTest(theta=0.0):
            theta = np.full((n_points,n_compon), 0.0)
            x = _newx(x0, ratio, theta, exp)
            npt.assert_equal(x0, x)
        with self.subTest(theta=1.0):
            theta = np.full((n_points,n_compon), 1.0)
            x = _newx(x0, ratio, theta, exp)
            npt.assert_equal(x0*ratio, x)
        with self.subTest(theta=0.5):
            theta = np.full((n_points,n_compon), 0.5)
            x = _newx(x0, ratio, theta, exp)
            npt.assert_equal(x0*(ratio*0.5+0.5), x)

    def test_update_theta(self):
        ratio = np.linspace(0.2, 2.0, 10)
        real = np.array([0.74, 0.58, 0.42, 0.26, 0.1, 0.23333333, 0.32857143,
                         0.4, 0.45555556, 0.5])
        calc = _update_theta(ratio)
        npt.assert_allclose(calc, real)

    def test_calc_ratio(self):
        stoichx = np.concatenate((np.eye(2, dtype=int), self.fixt4['stoich']), axis=0)
        st_sign = stoichx > 0
        st_pos = np.where(st_sign, stoichx, np.zeros_like(stoichx, dtype=int))
        st_neg = np.where(np.logical_not(st_sign), -stoichx, np.zeros_like(stoichx, dtype=int))
        abs_analc = np.abs(self.fixt4['analc'])
        calc_sump, calc_sumr = _calc_ratio(self.fixt4['conc'],
                                           st_pos, st_neg, abs_analc, self.fixt4['analc']<0.0)
        npt.assert_allclose(self.fixt4['true_sump'], calc_sump, atol=1e-10)
        npt.assert_allclose(self.fixt4['true_sumr'], calc_sumr, atol=1e-10)

    def test_precond_fixtures2(self):
        for name, fixture in self.fixtures2.items():
            with self.subTest(label=name):
                self.__success_precond(fixture)

    def __success_precond(self, fixture) -> None:
        DAMPING_THRESHOLD = 0.2
        # breakpoint()
        initial_guess = guess_starting_point(fixture)
        cex0 = species_concentration(initial_guess, fixture, full=True)
        fprev = mass_action(fixture, cex0)
        x, _ = precondition(fixture, initial_guess)
        cex1 = species_concentration(x, fixture, full=True)
        fpost = mass_action(fixture, cex1)
        # aux = np.abs(fpost/fprev)
        #print(np.argmax(aux), np.max(aux))
        self.assertTrue(np.all(DAMPING_THRESHOLD*np.abs(fprev)>np.abs(fpost)))

    def setUp(self):
        self.titr1 = self.__fixture2()
        self.titr2 = self.__fixture3()
        self.fixt4 = self.__fixture4()
        self.fixtures2 = all_fixtures

    def __fixture2(self):
        data = synthetic_simple1_wosld()
        return _read_fixture(data)

    def __fixture3(self):
        data = synthetic_simple1_wsld()
        return _read_fixture(data)

    def __fixture4(self):
        return {
            'conc': np.ones((6,4), dtype=float),
            'stoich': np.array([[ 1,1], [1,-1]], dtype=int),
            'analc': np.array([[0.0, 1.0],
                              [0.2, 0.6],
                              [0.4, 0.2],
                              [0.6, -0.2],
                              [0.8, -0.6],
                              [1.0, -1.0]]),
            'true_sumr': np.array([[ 3,  2],
                                   [3,  2],
                                   [3,  2],
                                   [3,  2.2],
                                   [3,  2.6],
                                   [3,  3]]),
            'true_sump': np.array([[0, 2],
                                   [0.2, 1.6],
                                   [0.4, 1.2],
                                   [0.6, 1],
                                   [0.8, 1],
                                   [1, 1]])
        }


def _read_fixture(data):
    return EquilibriumData(analytical_concentration=data.analyticalc,
                           log_beta=data.log_beta,
                           stoichiometry=data.stoichiometry,
                           log_ks=data.log_ks,
                           solid_stoichiometry=data.solid_stoichiometry)


if __name__ == '__main__':
    unittest.main()
