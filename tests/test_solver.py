#!/usr/bin/env -S uv run --script

import unittest

import numpy as np
import numpy.testing as npt

from equilx import Eq
from equilx.utils import mass_action

from fixtures2 import all_fixtures


class TestSolver(unittest.TestCase):
    def test_cases(self):
        for case, equil in self.fixtures.items():
            with self.subTest(case):
                equil, ref = self.fixtures['phosphoric']
                result = Eq(equil, freeze=ref)
                y = mass_action(equil,
                                result['free concentration'],
                                result.get('solid concentration', None))
                npt.assert_allclose(y, np.zeros_like(y), atol=1e-8)


    # def test_pp(self):
    #     equil, ref = self.fixtures['phosphoric']
    #     result = Eq(equil, freeze=ref)
    #     y = mass_action(equil,
    #                     result['free concentration'],
    #                     result.get('solid concentration', None))
    #     print("max = ", np.max(y), " min=", np.min(y))
    #     npt.assert_allclose(y, np.zeros_like(y), atol=1e-8)

    def setUp(self) -> None:
        self.fixtures = all_fixtures
        return super().setUp()


if __name__ == '__main__':
    unittest.main()
