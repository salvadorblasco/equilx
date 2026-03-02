#!/usr/bin/env -S uv run --script

import unittest

import numpy as np
import numpy.testing as npt

from equilx.data_structure import EquilibriumData
from equilx.solver.solver import EqSolver
from equilx.solver.utils import calc_saturation, make_titration, mass_action, species_concentration, unpack_eqsolve

from fixtures2 import case_aluminium, true_case_aluminium
from fixtures2 import case_phosphoric_calcium, case_molybdate


class TestSolver(unittest.TestCase):
    def test_phosphate_calcium(self):
        eq: EquilibriumData = case_phosphoric_calcium()
        crude_conc = EqSolver(eq, max_iterations=100, overshoot_capping=10.0)
        c_soluble = crude_conc[:,:3]
        # c_solids = crude_conc[:,4:]
        cext = species_concentration(c_soluble, eq, full=True)
        fobj = mass_action(eq, cext)
        # print(c_soluble)
        npt.assert_allclose(fobj, np.zeros_like(fobj), atol=1e-10)

    def test_aluminium(self):
        eq = case_aluminium()
        true_c_soluble, true_c_solid, true_si = true_case_aluminium()
        crude_conc = EqSolver(eq)
        d = unpack_eqsolve(crude_conc, eq)
        npt.assert_allclose(true_c_soluble[:,:eq.n_components], d['concentration components'], atol=1e-8)
        npt.assert_allclose(true_c_soluble[:,eq.n_components:], d['concentration species'], atol=1e-8)
        npt.assert_allclose(true_c_solid, d['concentration solids'], atol=1e-8)
        # calc_si = d['saturation index']
        # npt.assert_allclose(true_si, calc_si, atol=1e-8)

    def test_molybdate(self):
        eq = case_molybdate()
        crude_conc = EqSolver(eq)
        d = unpack_eqsolve(crude_conc, eq)
        y = mass_action(eq, d['concentration full'])
        npt.assert_allclose(y, np.zeros_like(y), atol=1e-11)


if __name__ == '__main__':
    unittest.main()
