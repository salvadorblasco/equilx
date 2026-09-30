"Module for potentiometry optimisation."

from contextlib import contextmanager
import os
import io
import itertools

import numpy as np

from equilx import fsolver
from equilx.data_structure import ModelData, DebyeHuckelParameters, PotentiometryData, PotentiometryTitration


def PotentiometryOpt(equil: ModelData, emfdata: PotentiometryData,
                     weighting_scheme: str,
                     temperature: float = 25.0,
                     max_its: int = 20):

    n_points = [p.number_of_points for p in emfdata]
    total_points = sum(n_points)
    slices = list(itertools.accumulate(n_points))
    slices.pop()

    #with capture_fortran_stdout() as captured:
    print("inititating")
    fsolver.emf.init_module(total_points,
                            emfdata.number_of_titrations,
                            np.asfortranarray(equil.stoichiometry),
                            np.asfortranarray(equil.log_beta),
                            np.array(equil.log_beta_flags, order='F'),
                            0, 6)

    for n, t in enumerate(emfdata):
        print(f"Loading titration #{n}")
        fsolver.emf.load_titration(*load_data(t))

    if weighting_scheme == 'auto':
        fsolver.emf.set_weight_gradient()
    elif weighting_scheme == 'unit':
        fsolver.emf.set_weight_unit()
    else:
        raise ValueError("weighting_scheme unknown")

    #if 'debye_parameters' in kwargs:
    #    if not isinstance(kwargs['debye_parameters'], DebyeHuckelParameters):
    #        raise TypeError

    #    fsolver.emf.include_ionic()

    fsolver.emf.fit_emf(max_its, temperature)

    # read variables and process them before finalizing
    result = {
        'free_concs': np.array_split(fsolver.emf.C, slices, axis=0),
        'log_beta': fsolver.emf.lbeta,
        'titr_init': fsolver.emf.titr_init,
        'titr_buret': fsolver.emf.titr_buret,
        'electr_emf0': fsolver.emf.electr_emf0,
        'electr_slope': fsolver.emf.electr_slope,
        'electr_junction': fsolver.emf.junctionp,
        'covariance': fsolver.emf.covariance,
        'message': fsolver.emf.emf_diagnostic
    }

    fsolver.emf.finalize()

    output = io.StringIO()
    output.write(captured)
    result['output'] = output

    print(result)
    return result


def load_data(ptit: PotentiometryTitration) -> tuple:
    def titration_parms():
        return np.asfortranarray(np.vstack([ptit.initial_amount,
                                            ptit.buret_concentration]).T)
    def titration_flags():
        return np.asfortranarray(np.stack([ptit.flags_initial, ptit.flags_buret]).T)

    return (np.asfortranarray(ptit.titre_added),
            np.asfortranarray(ptit.emf),
            titration_parms(),
            ptit.electroactive,
            ptit.starting_volume,
            ptit.emf0,
            ptit.slope,
            titration_flags(),
            ptit.flag_emf0,
            ptit.flag_slope,
            ptit.emf_error,
            ptit.volume_error,
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
