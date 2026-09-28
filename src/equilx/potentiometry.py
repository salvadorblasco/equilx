"Module for potentiometry optimisation."

from contextlib import contextmanager
from dataclasses import dataclass, field
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
                     **kwargs):
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

        if 'debye_parameters' in kwargs:
            if not isinstance(kwargs['debye_parameters'], DebyeHuckelParameters):
                raise TypeError

            fsolver.emf.include_ionic()

        fsolver.emf.fit_emf(max_its, temperature)

        # read variables and process them before finalizing

        fsolver.emf.finalize()

    if 'output' in kwargs:
        output: io.StringIO = kwargs['output']
        output.write(captured)


def load_data(ptit: PotentiometryTitration) -> tuple:
    return (np.array(ptit.titre, order='F'),
            np.array(ptit.emf_read, order='F'),
            ptit.__titration_parms(),
            ptit.electro_active,
            ptit.starting_volume,
            ptit.emf0,
            ptit.slope,
            ptit.__titration_flags(),
            ptit.flag_emf0,
            ptit.flag_slope,
            ptit.error_emf,
            ptit.error_vol,
            np.array(ptit.junction_potential, order='F'),
            np.array(ptit.junction_component, order='F'),
            np.array(ptit.flag_junction, order='F'))



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
