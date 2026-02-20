#!/usr/bin/env -S uv run --script

import unittest

import numpy as np

from equilx.data_structure import EquilibriumData


class TestEquilibriumData(unittest.TestCase):
    def test_creation(self):
        stoich = np.random.randint(0, 9, size=(5,3))
        sstoich = np.random.randint(0, 9, size=(2,3))
        analc = np.random.rand(50,3)
        logb = np.random.rand(5)
        logk = np.random.rand(2)
        eq = EquilibriumData(analytical_concentration=analc, log_beta=logb, stoichiometry=stoich, log_ks=logk, solid_stoichiometry=sstoich)


if __name__ == '__main__':
    unittest.main()
