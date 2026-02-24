from equilx.mytypes import FArray
from equilx.data_structure import DebyeHuckelParameters, EquilibriumData

def calc_ionic_strength(concentration, equil: EquilibriumData) -> FArray:
    return 0.5*concentration @ equil.chargesx**2


def beta_increment(ionic_strength: FArray, reference: float, dhparams: DebyeHuckelParameters) -> FArray:
    ...
