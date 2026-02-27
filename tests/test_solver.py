#!/usr/bin/env -S uv run --script

import unittest

import numpy as np
import numpy.testing as npt

from equilx.data_structure import EquilibriumData
from equilx.solver.solver import EqSolver
from equilx.solver.utils import calc_saturation, make_titration, mass_action, species_concentration, unpack_eqsolve

from fixtures2 import case_aluminium, true_case_aluminium
from fixtures2 import case_phosphoric_calcium


class TestSolver(unittest.TestCase):
    @unittest.skip("broken")
    def test_phosphate_calcium(self):
        eq: EquilibriumData = case_phosphoric_calcium()
        crude_conc = EqSolver(eq)
        c_soluble = crude_conc[:,:3]
        # c_solids = crude_conc[:,4:]
        cext = species_concentration(c_soluble, eq, full=True)
        fobj = mass_action(eq, cext)
        npt.assert_allclose(fobj, np.zeros_like(fobj), atol=1e-10)

    @unittest.skip("broken")
    def test_aluminium(self):
        eq = case_aluminium()
        true_c_soluble, true_c_solid, true_si = true_case_aluminium()
        crude_conc = EqSolver(eq)
        d = unpack_eqsolve(crude_conc, eq)
        #print(np.hstack((true_c_solid, d['concentration solids'])))
        calc_si = d['saturation index']
        # suspicious = (d['concentration solids'] < 1e-9) & (calc_si > 1.01)
        # print(suspicious)
        # _trr_si = calc_saturation(true_c_soluble[:,:2], eq)
        # print(np.hstack((true_si, d['saturation index'], _trr_si)))  #, true_c_solid, d['concentration solids'])))
        #breakpoint()
        # rconc = true_c_soluble[:,:eq.n_components] / d['concentration components']
        # print(np.unravel_index(np.argmax(np.abs(1/rconc)),rconc.shape), np.max(np.abs(1/rconc)))
        print(np.hstack((true_c_soluble[7,:eq.n_components], d['concentration components'][7,:])))
        npt.assert_allclose(true_c_soluble[:,:eq.n_components], d['concentration components'], atol=1e-8)
        npt.assert_allclose(true_c_soluble[:,eq.n_components:], d['concentration species'], atol=1e-8)
        npt.assert_allclose(true_c_solid, d['concentration solids'], atol=1e-8)
        npt.assert_allclose(true_si, calc_si, atol=1e-8)


if __name__ == '__main__':
    unittest.main()
