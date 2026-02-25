#!/usr/bin/env -S uv run --script

import unittest

import numpy as np
import numpy.testing as npt

from equilx.data_structure import EquilibriumData
from equilx.solver.solver import EqSolver
from equilx.solver.utils import make_titration, mass_action, species_concentration


class TestSolver(unittest.TestCase):
    def test_phosphate_calcium(self):
        analc = make_titration(initial_conc=np.array([1e-3, 1e-3, 4e-3]),
                               buret_conc=np.array([0.0, 0.0, -0.1]),
                               starting_vol=10.0,
                               titre=np.linspace(0.0, 1.0, 50))
        stoich = np.array([[1, 0, 1],
                           [1, 0, 2],
                           [1, 0, 3],
                           [1, 1, 1],
                           [0, 1,-1],
                           [0, 0,-1]], dtype=int)
        log_beta = np.array([11.7, 18.6, 20.6, 1.7, -12.6, -13.77])
        sstoich = np.array([[0, 1,-2],
                            [2, 3, 0],
                            [1, 1, 1]], dtype=int)
        log_ks = np.array([-5.0, -26.0, -7.0])
        charges = np.array([-3, 2, 1], dtype=int)
        eq = EquilibriumData(
            analytical_concentration=analc,
            log_beta=log_beta,
            stoichiometry=stoich,
            log_ks=log_ks,
            solid_stoichiometry=sstoich,
            charges=charges)
        crude_conc = EqSolver(eq)
        c_soluble = crude_conc[:,:3]
        # c_solids = crude_conc[:,4:]
        cext = species_concentration(c_soluble, eq, full=True)
        fobj = mass_action(eq, cext)
        npt.assert_allclose(fobj, np.zeros_like(fobj), atol=1e-10)


if __name__ == '__main__':
    unittest.main()
