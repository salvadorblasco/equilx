from typing import Protocol, Any, Callable

import numpy as np

from equilx.data_structure import ExperimentalData, ModelData
from equilx.mytypes import FArray
from equilx.consts import Flags


class OptData(Protocol):
    titrations: list[Any]


class OptimizeProtocol(Protocol):
    """
    Class to deal with the Levenberg-Marquardt fitting routine.
    """
    def __init__(self, data: OptData, model: ModelData, reporter: Callable):
        ...

    def accept_step(self) -> None:
        ...

    def final_result(self) -> dict[str, Any]:
        ...

    def matrices(self) -> tuple[FArray, FArray]:
        ...

    def reject_step(self) -> None:
        ...

    def size(self) -> tuple[int, int]:
        ...

    def trial_step(self, increments: FArray) -> None:
        ...

    def tmp_residual(self) -> FArray:
        ...

    def weights(self) -> FArray:
        ...


class OptimizeObject:
    def __init__(self, data: ExperimentalData, model: ModelData, reporter: Callable):
        self._data = data
        self._model = model
        self._reporter = reporter

        self._background_ion_concentration = np.concatenate([t.background_ions for t in data.titrations])
        self._variables: FArray
        self._freeconcentration: FArray | None = None

    def _titration_parameters(self):
        itx = iter(self._variables[self._dof_beta:].tolist())
        itd = iter(self._step[self._dof_beta:].tolist())

        def select(c, flags):
            x = c.copy()
            for n, i in enumerate(flags):
                if i == Flags.REFINE:
                    x[n] = next(itx) + next(itd)
            return x

        for titration in self._data.titrations:
            starting_concentration = select(titration.starting_concentration, titration.starting_concentration_flags)
            buret_concentration = select(titration.buret_concentration, titration.buret_concentration_flags)
            yield starting_concentration, buret_concentration

    def _analytical_concentration(self) -> FArray:
        aconc = []
        for titration, (c0, ct) in zip(self._data.titrations, self._titration_parameters()):
            aux = (c0[None, :] * titration.starting_volume + ct[None, :] * titration.get_titre[:, None]) / \
                (titration.starting_volume + titration.get_titre[:, None])
            aconc.append(aux)
        return np.concatenate(aconc, axis=0)

    def _background_concentration(self) -> FArray:
        bconc = [
            (titration.c0back * titration.v0 + titration.ctback * titration.get_titre[:, None]) / \
                (titration.v0 + titration.get_titre[:, None])
            for titration in self._data.titrations
        ]
        return np.concatenate(bconc, axis=0)

    def _calc_free_concs(self, initial=False, update=False) -> FArray:
        if initial or self._freeconcentration is None:
            _initial_guess = None
        else:
            _initial_guess = self._freeconcentration[:,:self._model.n_components]

        equil = self._model.get_equildata(self._analytical_concentration())
        equil.background_ions = self._background_ion_concentration
        equil.log_beta[...] = self._beta()

        try:
            c, *_ = EqSolver(equil, initial_guess=_initial_guess)
        except excepts.DivergedIonicStrengthWarning as divergwarn:
            raise divergwarn
        else:
            if update:
                self._freeconcentration = c
        return c


def ravel(x, y, flags):
    """Update values from one iterable with other iterable according to flags.

    This function does the opposite action than :func:`unravel`.

    Parameters:
        x (iterable): the original array values
        y (iterable): the updated values to be plugged into *x*.
        flags (sequence): flags indicating how to update *x* with *y*. Accepted
            values are:

            * 0: value is to be kept constant
            * 1: value is to be refined and the corresponding value from x
              will be substituted by the corresponding value from y.
            * >2: value is restrained. All places with the same number are
              refined together and the ratio between them is maintained.

    Yields:
        float: Raveled values.
    """
    # indices of the reference parameter for constraining
    ref_index = {i: flags.index(i) for i in range(2, 1 + max(flags))}
    ref_val = {}

    ity = iter(y)
    for i, f in enumerate(flags):
        if f == Flags.REFINE:  # refinable: return new value
            yield next(ity)
        elif f == Flags.CONSTANT:  # constant: return old value
            yield x[i]
        else:  # constrained: return or compute
            if i == ref_index[f]:
                val = next(ity)  # reference value: return new value
                ref_val[f] = val  # and store ref value
                yield val
            else:  # other: compute proportional value
                yield x[i] * ref_val[f] / x[ref_index[f]]
