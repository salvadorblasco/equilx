from equilx.mytypes import FArray
from equilx.data_structure import DebyeHuckelParameters, EquilibriumData

def calc_ionic_strength(concentration, equil: EquilibriumData) -> FArray:
    return 0.5*concentration @ equil.chargesx**2


def beta_increment(ionic_strength: FArray, reference: float, dhparams: DebyeHuckelParameters) -> FArray:
    ...


def debye_huckel(ionic: FArray, equil: EquilibriumData, dbparms: DebyeHuckelParameters) -> FArray:
    term1 = -equil.zstar[None,:]*dbparms.a*ionic[:,None]**0.5/(1+dbparms.b*ionic[:,None]**0.5)
    C = dbparms.c0*equil.pstar + dbparms.c1*equil.zstar
    term2 = C*ionic[:,None]
    D = dbparms.d0*equil.pstar + dbparms.d1*equil.zstar
    term3 = D*ionic[:,None]**1.5
    E = dbparms.e0*equil.pstar + dbparms.e1*equil.zstar
    term4 = E*ionic[:,None]**1.5
    total = term1+term2+term3+term4
    return total.squeeze()
