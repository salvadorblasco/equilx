#!/usr/bin/env -S uv run --script

"Test collection for potentiometry data fitting."

import unittest
import random

import numpy as np
import numpy.testing as npt

from equilx.potentiometry import PotentiometryOpt

from fixtures_potentiometry import case_znedta


class Test_ZnEDTA(unittest.TestCase):
    """
    Test potentiometry optimisation using synthetic Zn/EDTA data
    """

    def setUp(self) -> None:
        self.data, self.model = case_znedta()
        #self.data.potentiometry_opts.weights = "calculated"
        self.true_beta = np.array([ 16.25, 19.25, 4.65])
        return super().setUp()

    # def test_first(self):
    #     potobj = PotentiometryOpt(self.data, self.model, _dummy_reporter)
    #     status, lminfo = levenberg_marquardt(potobj)
    #     result = potobj.final_result()
    #     result.update(lminfo)
    #     self.assertEqual(status, Exec.NORMAL_END)
    #     calc_beta = result['final log beta'][-4:-1]
    #     npt.assert_allclose(calc_beta, self.true_beta, rtol=0.01)
        
    def test_noisy_beta(self):
        # @pytest.mark.parametrize("error", [0.1, 0.2, 0.5, 1.0])
        error = 0.1
        sd = self.copy_data()
        sd.log_beta[-4:-1] += error*2*(np.random.rand(3)-0.5)
        result = PotentiometryOptimizer(sd, reporter=_dummy_reporter)
        calc_beta = result['final log beta'][-4:-1]
        npt.assert_allclose(calc_beta, self.true_beta, rtol=0.01)

    # @pytest.mark.parametrize("error", [0.0002, 0.0005, 0.001, 0.002])
    # def test_noisy_c0(self, error):
    #     sd = self.copy_data()
    #     true_c0 = []
    #     for tit in sd.potentiometry_opts.titrations:
    #         true_c0.append(tit.c0.copy())
    #         tit.c0_flags = [Flags.CONSTANT, Flags.CONSTANT, Flags.REFINE]
    #         tit.c0[-1] += random.uniform(-error/2, error/2) 
    #     result = PotentiometryOptimizer(sd, reporter=_dummy_reporter)
    #     calc_beta = result['final log beta'][-4:-1]
    #     npt.assert_allclose(calc_beta, self.true_beta, rtol=0.01)
    #     for (calc_c0, _), _true_c0 in zip(result['final titration parameters'], true_c0):
    #         npt.assert_allclose(calc_c0, _true_c0, atol=1e-5)


def _dummy_reporter(**kwargs):
    pass
