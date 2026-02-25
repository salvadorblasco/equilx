#!/usr/bin/env -S uv run --script

import unittest
import numpy as np
import numpy.testing as npt

from equilx.data_structure import EquilibriumData
from equilx.solver.nr import objective_function, jacobian, newton_raphson

from equilx.solver.utils import mass_action, species_concentration
import fixtures
from fixtures2 import large_fixture1

np.seterr(all='ignore')


class TestNR(unittest.TestCase):
    def test_nr1(self):
        real_freec = self.wsld.real_freec[:,(0,1,3)]
        equil = _read_fixture(self.wsld)
        for noise in (0.001, 0.01, 0.1):
            x0 = (real_freec + 1e-6) * np.random.normal(1.0, noise, (5,3))
            with self.subTest(noise=noise):
                calc_freec = newton_raphson(x0, equil)
                npt.assert_allclose(real_freec, calc_freec, atol=1e-10)

    def test_nr2(self):
        real_freec = self.wosld.real_freec[:,:2]
        equil = _read_fixture(self.wosld)
        for noise in (0.001, 0.01, 0.1):
            x0 = real_freec * np.random.normal(1.0, noise, (5,2))
            with self.subTest(noise=noise):
                calc_freec = newton_raphson(x0, equil)
            npt.assert_allclose(real_freec, calc_freec, atol=1e-10)

    def test_large(self):
        eq = large_fixture1()
        with np.load("_large_fixture.npz") as f:
            x0 = f['arr_0']
        x = newton_raphson(x0, eq)
        md = mass_action(eq, species_concentration(x, eq, full=True))
        npt.assert_allclose(md, np.zeros_like(md), atol=1e-15)

    def setUp(self) -> None:
        self.wsld = fixtures.synthetic_simple1_wsld()
        self.wosld = fixtures.synthetic_simple1_wosld()
        return super().setUp()


def _read_fixture(data):
    return EquilibriumData(analytical_concentration=data.analyticalc,
                           log_beta=data.log_beta,
                           stoichiometry=data.stoichiometry,
                           log_ks=data.log_ks,
                           solid_stoichiometry=data.solid_stoichiometry)


class TestOF(unittest.TestCase):
    def test_of1(self):
        equil = _read_fixture(self.wsld)
        f = objective_function(equil, self.wsld.real_freec[:,:3], self.wsld.real_freec[:,3,None])
        npt.assert_allclose(f, np.zeros_like(f), atol=5e-9)
        
    def test_of2(self):
        equil = _read_fixture(self.wosld)
        f = objective_function(equil, self.wosld.real_freec, np.zeros((5,0)))
        npt.assert_allclose(f, np.zeros_like(f), atol=5e-9)

    def setUp(self) -> None:
        self.wsld = fixtures.synthetic_simple1_wsld()
        self.wosld = fixtures.synthetic_simple1_wosld()
        return super().setUp()
    

class TestJac(unittest.TestCase):
    def test_jac1(self):
        equil = _read_fixture(self.wsld)
        caljac = jacobian(self.wsld.real_freec[:,:3],
                          self.wsld.real_freec[:, 3, None],
                          equil, self.wsld.gfun)
        npt.assert_allclose(caljac, self.wsld.real_jac)

    def test_jac2(self):
        equil = _read_fixture(self.wosld)
        caljac = jacobian(self.wosld.real_freec,
                          np.zeros((5,0)),
                          equil,
                          self.wosld.gfun)
        npt.assert_allclose(caljac, self.wosld.real_jac)

    def setUp(self) -> None:
        self.wsld = fixtures.synthetic_simple1_wsld()
        self.wosld = fixtures.synthetic_simple1_wosld()
        return super().setUp()


if __name__ == '__main__':
    unittest.main()
