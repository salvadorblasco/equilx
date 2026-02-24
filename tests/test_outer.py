#!/usr/bin/env -S uv run --script

import unittest
import numpy as np
import numpy.testing as npt

from equilx.data_structure import EquilibriumData, DebyeHuckelParameters
from equilx.solver.precond import precondition
from equilx.solver.outer_fixed_point import outer_fixed_point
from equilx.solver.utils import mass_action

from fixtures2 import large_fixture1


class TestOuter(unittest.TestCase):
    def test_outer(self):
        x0 = precondition(self.eq)
        np.savez_compressed("_large_fixture.npz", x0)
        x = outer_fixed_point(x0, self.eq, self.dhp)
        ma = mass_action(self.eq, x)
        npt.assert_allclose(ma, np.zeros_like(ma))

    def setUp(self) -> None:
        eq = large_fixture1()
        eq.background_ionic = 0.1
        eq.reference_ionic = 0.1
        self.eq = eq
        self.dhp = DebyeHuckelParameters()
        return super().setUp()


if __name__ == '__main__':
    unittest.main()
