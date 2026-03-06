#!/usr/bin/env -S uv run --script

import unittest

import numpy as np
import numpy.testing as npt

from equilx.solver.solver import EqSolver
from equilx.solver.utils import calc_saturation, make_titration, mass_action, species_concentration, unpack_eqsolve

from fixtures2 import all_fixtures


class TestSolver(unittest.TestCase):
    def test_cases(self):
        for case, equil in self.fixtures.items():
            with self.subTest(case):
                crude_conc = EqSolver(equil)
                d = unpack_eqsolve(crude_conc, equil)
                y = mass_action(equil, d['concentration full'])
                npt.assert_allclose(y, np.zeros_like(y), atol=1e-10)

    def setUp(self) -> None:
        self.fixtures = all_fixtures
        return super().setUp()


if __name__ == '__main__':
    unittest.main()
