#!/usr/bin/env -S uv run --script

import unittest
import numpy as np
import numpy.testing as npt

from equilx.solver.nr import objective_function, jacobian, newton_raphson

np.seterr(all='ignore')


class TestNR(unittest.TestCase):
    def test_nr1(subtests, synthetic_simple1_wsld):
        real_freec = synthetic_simple1_wsld.real_freec[:,(0,1,3)]
        print(real_freec)
        for noise in (0.001, 0.01, 0.1):
            x0 = (real_freec + 1e-6) * np.random.normal(1.0, noise, (5,3))
            with subtests.test(msg=f"{noise=}"):
                calc_freec = newton_raphson(x0, synthetic_simple1_wsld.log_beta,
                                            synthetic_simple1_wsld.log_ks,
                                            synthetic_simple1_wsld.stoichiometry,
                                            synthetic_simple1_wsld.solid_stoichiometry,
                                            synthetic_simple1_wsld.analyticalc)
                npt.assert_allclose(real_freec, calc_freec, atol=1e-10)


    def test_nr2(subtests, synthetic_simple1_wosld):
        real_freec = synthetic_simple1_wosld.real_freec[:,:2]
        for noise in (0.001, 0.01, 0.1):
            x0 = real_freec * np.random.normal(1.0, noise, (5,2))
            with subtests.test(msg=f"{noise=}"):
                calc_freec = newton_raphson(x0, synthetic_simple1_wosld.log_beta,
                                            synthetic_simple1_wosld.log_ks,
                                            synthetic_simple1_wosld.stoichiometry,
                                            synthetic_simple1_wosld.solid_stoichiometry,
                                            synthetic_simple1_wosld.analyticalc)
            npt.assert_allclose(real_freec, calc_freec, atol=1e-10)


class TestOF(unittest.TestCase):
    def test_of1(synthetic_simple1_wsld):
        stoichx = np.vstack([np.eye(2), synthetic_simple1_wsld.stoichiometry])
        f = objective_function(synthetic_simple1_wsld.analyticalc,
                               synthetic_simple1_wsld.real_freec[:,:3],
                               stoichx,
                               synthetic_simple1_wsld.real_freec[:,3,None],
                               synthetic_simple1_wsld.solid_stoichiometry, 
                               10**synthetic_simple1_wsld.log_ks)
        npt.assert_allclose(f, np.zeros_like(f), atol=5e-9)
        

    def test_of2(synthetic_simple1_wosld):
        stoichx = np.vstack([np.eye(2), synthetic_simple1_wosld.stoichiometry])
        f = objective_function(synthetic_simple1_wosld.analyticalc,
                               synthetic_simple1_wosld.real_freec,
                               stoichx,
                               np.zeros((5,0)),
                               synthetic_simple1_wosld.solid_stoichiometry, 
                               10**synthetic_simple1_wosld.log_ks)
        npt.assert_allclose(f, np.zeros_like(f), atol=5e-9)
    

class TestJac(unittest.TestCase):
    def test_jac1(synthetic_simple1_wsld):
        stoichx = np.vstack([np.eye(2), synthetic_simple1_wsld.stoichiometry])
        caljac = jacobian(synthetic_simple1_wsld.real_freec[:,:3],
                          stoichx,
                          synthetic_simple1_wsld.real_freec[:, 3, None],
                          synthetic_simple1_wsld.solid_stoichiometry,
                          synthetic_simple1_wsld.gfun)
        npt.assert_allclose(synthetic_simple1_wsld.real_jac, caljac)


    def test_jac2(synthetic_simple1_wosld):
        stoichx = np.vstack([np.eye(2), synthetic_simple1_wosld.stoichiometry])
        caljac = jacobian(synthetic_simple1_wosld.real_freec,
                          stoichx,
                          np.zeros((5,0)),
                          synthetic_simple1_wosld.solid_stoichiometry,
                          synthetic_simple1_wosld.gfun)
        npt.assert_allclose(synthetic_simple1_wosld.real_jac, caljac)


if __name__ == '__main__':
    unittest.main()
