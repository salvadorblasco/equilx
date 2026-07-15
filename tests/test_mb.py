#!/usr/bin/env -S uv run --script

import unittest

import numpy as np
import numpy.testing as npt

from equilx.data_structure import EquilibriumData
from equilx.solver.utils import mass_action, species_concentration
import equilx.solver.mb as mb

import fixtures2


class TestMassBalance(unittest.TestCase):
    def test_phosphoric(self):
        eq: EquilibriumData  = fixtures2.case_phosphoric()
        initial_guess = np.full((eq.n_points, eq.n_components), 1e-12)
        fconc = species_concentration(initial_guess, eq, full=True)
        mb1 = mass_action(eq, fconc)
        mb2 = np.empty_like(mb1, order='F')
        mb.mb_soluble_tit(mb2,
                          np.asfortranarray(fconc),
                          np.asfortranarray(eq.stoichiometry),
                          np.asfortranarray(eq.analytical_concentration))
        npt.assert_allclose(mb2, mb1)

    def test_aluminium(self):
        eq: EquilibriumData  = fixtures2.case_aluminium()
        true_c, true_s, _ = fixtures2.true_case_aluminium()
        with self.subTest('python mass balance'):
            mb1 = mass_action(eq, true_c, true_s)
            npt.assert_allclose(mb1, np.zeros_like(mb1), atol=1e-10)
        with self.subTest('fortran mass balance'):
            mb2 = np.empty_like(mb1, order='F')
            mb.mb_solids_tit(mb2,
                             np.asfortranarray(true_c),
                             np.asfortranarray(eq.stoichiometry),
                             np.asfortranarray(eq.analytical_concentration),
                             np.asfortranarray(eq.solid_stoichiometry),
                             np.asfortranarray(true_s))
            npt.assert_allclose(mb2, np.zeros_like(mb2))


if __name__ == '__main__':
    unittest.main()

