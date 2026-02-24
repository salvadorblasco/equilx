#!/usr/bin/env -S uv run --script

import unittest
import numpy as np
import numpy.testing as npt

from equilx.data_structure import EquilibriumData, DebyeHuckelParameters
from equilx.solver.outer_fixed_point import outer_fixed_point
from equilx.solver.utils import mass_action, species_concentration

from fixtures2 import large_fixture1


class TestOuter(unittest.TestCase):
    def test_outer(self):
        x, _ = outer_fixed_point(self.initial_guess, self.eq, self.dhp)
        cex = species_concentration(x, self.eq, full=True)
        ma = mass_action(self.eq, cex)
        npt.assert_allclose(ma, np.zeros_like(ma), atol=1e-10)

    def setUp(self) -> None:
        eq: EquilibriumData = large_fixture1()
        eq.background_ions = np.linspace(0.1, 0.15, eq.n_points)
        eq.reference_ionic = 0.1
        self.eq = eq
        self.dhp = DebyeHuckelParameters()
        with np.load("_large_fixture.npz") as f:
            self.initial_guess = f['arr_0']

        return super().setUp()


if __name__ == '__main__':
    unittest.main()
