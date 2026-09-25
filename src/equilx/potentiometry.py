
"Module for potentiometry optimisation."

from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Optional
import os
import io

import numpy as np

from equilx import fsolver
from equilx.data_structure import EquilibriumData, DebyeHuckelParameters, PotentiometryTitration


def PotentiometryOpt(equil: EquilibriumData, emfdata: list[PotentiometryTitration],
                     flags_beta: list[int],
                     weighting_scheme: str,
                     temperature: float = 25.0,
                     max_its: int = 20,
                     output: io.StringIO|None = None):
    n_points = sum(p.n_points() for p in emfdata)
    len(p.junction_potential)

    with capture_fortran_stdout() as captured:
        fsolver.emf.init_module(n_points, len(emfdata), np.asfortranarray(equil.stoichiometry),
                                np.asfortranarray(equil.log_beta),
                                np.array(flags_beta, order='F'),
                                0, 6)

        for t in emfdata:
            fsolver.emf.load_titration(*t.load_data())

        if weighting_scheme == 'auto':
            fsolver.emf.set_weight_gradient()
        elif weighting_scheme == 'unit':
            fsolver.emf.set_weight_unit()
        else:
            raise ValueError("weighting_scheme unknown")

        fsolver.emf.fit_emf(max_its, temperature)

        # read variables and process them before finalizing

        fsolver.emf.finalize()

    if output is not None:
        output.write(captured)


@dataclass
class PotentiometricTitration:
    emf_read: list[float]
    titre: list[float]
    init: list[float]
    buret: list[float]
    emf0: float
    electro_active: int
    starting_volume : float
    slope: float = field(default=1.0)
    error_emf: float = field(default=0.3)
    error_vol: float = field(default=0.003)
    junction_potential: list[float] = field(default_factory=list)
    junction_component: list[int] = field(default_factory=list)
    flag_init: list[int] = field(init=False)
    flag_buret: list[int] = field(init=False)
    flag_emf0: int = field(init=False,default=0)
    flag_slope: int = field(init=False,default=0)
    flag_junction: list[int] = field(default_factory=list)

    def __post_init__(self):
        if len(self.emf_read) != len(self.titre):
            raise ValueError("emf and titre size must match")

        if len(self.init) != len(self.buret):
            raise ValueError("init and buret size must match")

    def add_junction(self, jpot: float, component: int, flag: int=0):
        self.junction_potential.append(jpot)
        self.junction_component.append(component)
        self.flag_junction.append(flag)

    def load_data(self) -> tuple:
        return (np.array(self.titre, order='F'),
                np.array(self.emf_read, order='F'),
                self.__titration_parms(),
                self.electro_active,
                self.starting_volume,
                self.emf0,
                self.slope,
                self.__titration_flags(),
                self.flag_emf0,
                self.flag_slope,
                self.error_emf,
                self.error_vol,
                np.array(self.junction_potential, order='F'),
                np.array(self.junction_component, order='F'),
                np.array(self.flag_junction, order='F'))

    def n_points(self) -> int:
        return len(self.titre)

    def __titration_parms(self):
        aux = np.array([self.init, self.buret]).T
        return np.asfortranarray(aux)

    def __titration_flags(self):
        aux = np.array([self.flag_init, self.flag_buret]).T
        return np.asfortranarray(aux)

@contextmanager
def capture_fortran_stdout():
    # Create a real pipe
    r, w = os.pipe()
    # Save original stdout
    old_stdout = os.dup(1)
    try:
        # Point FD 1 (stdout) at the write end of the pipe
        os.dup2(w, 1)
        os.close(w)          # we no longer need the original write end
        yield                # <-- call your f2py Fortran routine here
    finally:
        # Restore original stdout
        os.dup2(old_stdout, 1)
        os.close(old_stdout)
        # Read everything that was written
        with os.fdopen(r, 'r') as f:
            data = f.read()
    return data
