#!/usr/bin/env -S uv run --script

import unittest

import numpy as np
import numpy.testing as npt

from fixtures2 import case_aluminium, true_case_aluminium

import equilx.solver.si as silib


class TestSaturationIndex(unittest.TestCase):
    def test_aluminium(self):
        eqdata = case_aluminium()
        true_c_soluble, _, true_si =  true_case_aluminium()

        mydim = (eqdata.n_points, eqdata.n_components)
        free_c = np.array(true_c_soluble[:,:eqdata.n_components], order='F')
        solid_stoic = np.array(eqdata.solid_stoichiometry, order='F')
        assert free_c.shape == mydim
        sip = silib.saturation_product_tit(free_c, solid_stoic)
        calc_si = sip / (10**eqdata.log_ks[None,:])
        npt.assert_allclose(calc_si, true_si, atol=1e-5)


if __name__ == '__main__':
    unittest.main()
