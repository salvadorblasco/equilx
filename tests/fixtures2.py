"Fixtures for ionic strength fixing"

import numpy as np

def case_phosphoric():
    labels = ('H', 'FOS')
    log_beta = np.array([12.38, 19.578, 21.726, -13.73])
    stoich = np.array([[1,1],[1,2],[1,3],[0,-1]], dtype=int)
    charge = np.array([1,-3], dtype=int)
    background_ions = 0.1
    titre_ions = 0.1
    return stoich, log_beta
