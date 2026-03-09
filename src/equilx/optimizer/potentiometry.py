"Routine for potentiometry data fitting."

from typing import Any, Callable

import numpy as np

from equilx.optimizer.opt import OptimizeObject
from equilx.data_structure import ModelData, PotentiometryData
from equilx.solver.solver import EqSolver
from equilx.consts import Flags, LN10
from equilx.mytypes import FArray

from .. import excepts
from . import jacobian
from . import libemf
from . import libfit
from . import libmath

from equilx.optimizer.opt import OptimizeObject


class PotentiometryOpt(OptimizeObject):
    def __init__(self, data: PotentiometryData, model: ModelData, reporter: Callable) -> None:
        super().__init__(data, model, reporter)

        # calculate degrees of freedom
        self._dof_beta = model.degrees_of_freedom
        self._dof_conc = data.get_dof()
        self._slices = np.fromiter(data.get_slices(), dtype=float)
        self._slopes = np.fromiter(data.get_slopes(), dtype=float)
        self._emf0 = np.fromiter(data.get_emf0(), dtype=float)
        self._dof = self._dof_beta + self._dof_conc

        ##> self._hindices
        self._hindices = np.fromiter(data.get_per_point_electroactive(), dtype=int)

        ##> self._experimental_remf
        self._experimental_emf = data.get_collapsed_emf()

        self._bmatrixt = np.concatenate([jacobian.bmatrix_t(t.get_titre(), t.starting_volume, model.n_components)
                                         for t in data.titrations])

        self._bmatrixb = np.concatenate([jacobian.bmatrix_b(t.get_titre, t.starting_volume, model.n_components)
                                         for t in data.titrations])

        self._weights = self.__create_weights()

        # initial variable vector
        self._idx_refinable = []
        self._idx_refinable_beta: tuple = tuple(n for n, f in enumerate(model.log_beta_flags) if f == Flags.REFINE)
        self._slice_betas = slice(0, sum(1 for _ in self._idx_refinable_beta if _))
        self._any_beta_refined = any(self._idx_refinable_beta)
        self._idx_refinable.extend(self._idx_refinable_beta)
        beta_to_refine = np.extract(self._idx_refinable_beta, model.log_beta)
        concs_to_refine = []

        for titration in data.titrations:
            if titration.starting_concentration_flags:
                assert len(titration.starting_concentration_flags) == model.n_components
                idx_refinable_c0 = tuple(n for n, f in enumerate(titration.starting_concentration_flags) if f == Flags.REFINE)
            else:
                idx_refinable_c0 = model.n_components*[False]
            self._idx_refinable.extend(idx_refinable_c0)
            if any(idx_refinable_c0):
                concs_to_refine.extend(np.extract(idx_refinable_c0, titration.starting_concentration).tolist())

            if titration.buret_concentration:
                assert len(titration.buret_concentration_flags) == model.n_components
                idx_refinable_ct = tuple(n for n, f in enumerate(titration.buret_concentration_flags) if f == Flags.REFINE)
            else:
                idx_refinable_ct = model.n_components*[False]
            self._idx_refinable.extend(idx_refinable_ct)
            if any(idx_refinable_ct):
                concs_to_refine.append(np.extract(idx_refinable_ct, titration.buret_concentration).tolist())

        self._any_conc_refined = (len(concs_to_refine) > 0)

        # very important !! the betas are in natural logarithm!
        self._variables = np.concatenate([beta_to_refine*LN10, np.array(concs_to_refine)])
        self._step = np.zeros(self._dof, dtype=float)
        self._previous_values = np.empty_like(self._variables)

        self._initial_titration_paras = list(self._titration_parameters())

        # Calculate free concetrations initially
        self._freeconcentration = self._calc_free_concs(initial=True, update=False)

    def accept_step(self) -> None:
        "Update the variables values and reset increments to 0.0."
        self._previous_values[:] = self._variables
        self._variables += self._step
        self._step[...] = 0.0

    def reject_step(self) -> None:
        """
        The step is rejected and increments are reset to 0.0 without updating the variables.
        """
        self._step[...] = 0.0

    def final_result(self, stdev=None) -> dict[str, Any]:
        variables = self._variables.copy()
        variables[self._slice_betas] /= LN10

        if stdev is not None:
            assert len(stdev) == len(self._variables)
            assert len(stdev) == self._dof
            istd = iter(stdev.tolist())

            err_log_beta = [next(istd)/LN10 if f == Flags.REFINE else None
                            for f in self._model.log_beta_flags]

            err_titr_parms = [
                [
                    [next(istd) if c0f == Flags.REFINE else None for c0f in t.c0_flags],
                    [next(istd) if ctf == Flags.REFINE else None for ctf in t.ct_flags]
                ] for t in self._titrations() ]
        else:
            err_log_beta = None
            err_titr_parms = None

        # fvals = self.final_values()
        eactive = libemf.hselect(self._freeconcentration, self._hindices) 
        
        retval = {
            'final variables': variables,
            'final log beta': self._beta(),
            'error log beta': err_log_beta,
            'initial titration parameters': self._initial_titration_paras,
            'final titration parameters': list(self._titration_parameters()),
            'error titration parameters': err_titr_parms,
            'free concentration': self._freeconcentration,
            'slices': self._slices,
            'total concentration': self._analytical_concentration(),
            'read emf': self._experimental_emf,
            'eactive': eactive,
            'background ion concentration': self._background_concentration(),
            'weights': self._weights,
            'degrees of freedom': self._dof,
            'number of experimental points': self._data.total_experimental_points,
            'variable names': list(self._name_variables())
        }
        return retval

    def matrices(self) -> tuple[FArray, FArray]:
        """
        Compute the jacobian and the residual arrays for the accepted step.
        """
        # 1. calculate free concentrations
        freec = self._calc_free_concs(initial=False, update=True)
        assert freec.shape == (self._data.total_experimental_points, self._model.n_components + self._model.n_equils)

        # 2. calculate A
        amatrix = jacobian.amatrix(freec, self._model.stoichiometryx)
        assert amatrix.shape == (self._data.total_experimental_points, self._model.n_components, self._model.n_components)

        # 3. calculate jacobian part referring to beta
        dlogc_dlogbeta = jacobian.dlogcdlogbeta(amatrix, freec, self._model.solid_stoichiometry)
        assert dlogc_dlogbeta.shape == (self._data.total_experimental_points, self._model.n_components, self._model.n_equils)
        jbeta = dlogc_dlogbeta

        # 4. calculate jacobian part referring to titration parameters
        _jc0 = jacobian.solve_xmatrix(amatrix, self._bmatrixt)
        _jct = jacobian.solve_xmatrix(amatrix, self._bmatrixb)
        jtit = np.zeros((self._data.total_experimental_points, self._model.n_components, 2*self._data.number_of_titrations*self._model.n_components), dtype=float)
        _js = [np.concatenate((_jc0[s], _jct[s]), axis=2) for s in self._slices]
        for n, s1 in enumerate(self._slices):
            s2 = slice(n*2*self._model.n_components, (n+1)*2*self._model.n_components)
            jtit[s1, :, s2] = _js[n]

        # 5. compute the total jacobian
        jac =  np.concatenate([jbeta, jtit], axis=2)
        assert jac.shape == (self._data.total_experimental_points,
                             self._model.n_components,
                             self._model.n_equils + 2*self._model.n_components*self._data.number_of_titrations)

        # 6. remove non refined parts
        trimmed_jac1 = jac[..., self._idx_refinable]
        trimmed_jac2 = trimmed_jac1[np.arange(self._data.total_experimental_points), self._hindices]
        trimmed_jac3 = -LN10*self._slopes[:, None] * trimmed_jac2

        # 7. compute residual
        residual = self.__calculate_residual(freec)

        return trimmed_jac3, residual

    def relative_change(self, step):
        return step/self._variables

    def report_raw(self, text):
        print(text)

    def report_step(self, **kwargs):
        """
        Pass refinement parameters on each iteration to report.
        """
        kwargs['log_beta'] = self._beta()
        kwargs['previous log beta'] = self._previous_values[self._slice_betas]/LN10
        kwargs['increment'] /= LN10
        kwargs['stoichiometry'] = self._model.stoichiometry
        kwargs['any beta refined'] = self._any_beta_refined
        kwargs['any conc refined'] = self._any_conc_refined
        kwargs['titration params'] = list(self._titration_parameters())
        self._reporter(**kwargs)

    def size(self) -> tuple[int, int]:
        "Return number of points, number of variables."
        return self._data.total_experimental_points, self._dof

    def trial_step(self, increments: FArray) -> None:
        if increments.shape != self._step.shape:
            raise ValueError(f"Shape mismatch: {increments.shape} != {self._step.shape}")
        self._step[:] = increments[:]

    def tmp_residual(self) -> FArray:
        """
        Calculate and return residual for the trial step.
        """
        freec = self._calc_free_concs(initial=False, update=False)
        return self.__calculate_residual(freec)

    def weights(self) -> FArray:
        return np.diag(self._weights)

    @property
    def degrees_of_freedom(self) -> int:
        """
        The number of variables to refine.
        """
        return self._data.total_experimental_points - self._dof

    def _analytical_concentration(self) -> FArray:
        aconc = []
        for titration, (c0, ct) in zip(self._titrations(), self._titration_parameters()):
            aux = (c0[None, :] * titration.v0 + ct[None, :] * titration.get_titre[:, None]) / \
                (titration.v0 + titration.get_titre[:, None])
            aconc.append(aux)
        return np.concatenate(aconc, axis=0)

    def _background_concentration(self) -> FArray:
        bconc = [
            (titration.c0back * titration.v0 + titration.ctback * titration.get_titre[:, None]) / \
                (titration.v0 + titration.get_titre[:, None])
            for titration in self._titrations()
        ]
        return np.concatenate(bconc, axis=0)

    def _beta(self):
        """
        Return log10(beta) for the trial step.
        """
        beta = self._model.log_beta.copy()
        beta[self._idx_refinable_beta] = (self._variables[self._slice_betas] + self._step[self._slice_betas]) / LN10
        return beta

    # def _stoichiometry(self, extended=False):
    #     "Get stoichiometry array."
    #     number_components = self._model.n_components
    #     if extended:
    #         return np.vstack((np.eye(number_components, dtype=int),
    #                           np.array(self._model.stoichiometry.T)))
    #     return self._model.stoichiometry.T

    # def _titrations(self):
    #     """
    #     Iterate over the titrations.
    #     """
    #     yield from iter(self._data.potentiometry_opts.titrations)

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
        
    def _calc_free_concs(self, initial=False, update=False) -> FArray:
        if initial or self._freeconcentration is None:
            _initial_guess = None
        else:
            _initial_guess = self._freeconcentration[:,:self._model.n_components]
        log_beta = self._beta()
        total_concentration = self._analytical_concentration()

        #charges = self._data.charges
        background_ions_concentration = self._background_concentration()
        independent_component_activity = None

        try:
            c, *_ = EqSolver(self._model.get_equildata(self._analytical_concentration()))
        except excepts.DivergedIonicStrengthWarning as divergwarn:
            raise divergwarn
        else:
            if update:
                self._freeconcentration = c
        return c

    def _name_variables(self):
        if self._any_beta_refined:
            species_names = self._data.species_names
            beta_flags = self._data.potentiometry_opts.beta_flags
            for name, flag in zip(species_names, beta_flags):
                if flag == Flags.REFINE:
                    yield f"logB[{name}]"
        if self._any_conc_refined:
            components = self._data.components
            for ntit, titration in enumerate(self._titrations()):
                for comp, flag in zip(components, titration.c0_flags):
                    if flag == Flags.REFINE:
                        yield f"c0[{comp}, tit#{ntit}]"
                for comp, flag in zip(components, titration.ct_flags):
                    if flag == Flags.REFINE:
                        yield f"cT[{comp}, tit#{ntit}]"

    def __calculate_residual(self, free_concentrations):
        assert free_concentrations.shape == (self._total_points, self._nspecies + self._ncomponents)
        eactive = libemf.hselect(free_concentrations, self._hindices) 
        calculated_emf = self._emf0 + self._slopes * np.log(eactive)
        assert calculated_emf.shape == (self._total_points,)
        residual = self._experimental_emf - calculated_emf
        return residual

    def __create_weights(self):
        weights_mode = self._data.potentiometry_opts.weights
        if weights_mode == 'calculated':
            retv = np.concatenate([libemf.emf_weights(t.get_titre, t.v0_sigma, t.get_emf, t.e0_sigma)
                                   for t in self._titrations()])
        elif weights_mode == 'constants':
            retv = np.ones_like(self._experimental_emf, dtype=float)
        else:
            raise NotImplementedError(f"weight mode '{weights_mode}' not implemented")

        if np.any(np.isnan(retv)):
            raise ValueError("Some calculated weight values are NaN")

        return retv


def PotentiometryOptimizer(data: SolverData, reporter=None) -> dict[str, Any]:
    """
    Solve a potentiometry problem. Refine constants and possibly, concentrations.

    Parameters:
    -------
    data : SolverData
        The data for the refinement.

    Returns:
        x : 
            the refined data
        concs :
            the final free concentrations
        final_log_beta :
            the final refined constant values
        b_error :
            the fitting error in the concentration
        cor_matrix :
            the correlation matrix
        cov_matrix :
            the covariance matrix
        return_extra :
            additional information
    -------
    """
    bridge: Bridge = PotentiometryBridge(data, reporter)
    fit_status, fit_result = libfit.levenberg_marquardt(bridge, debug=False)

    stdev, cor, cov = fit_final_calcs(fit_result['jacobian'], fit_result['residuals'], bridge.weights()) 

    retval = bridge.final_result(stdev)
    retval.update(fit_result)
    retval['covariance'] = cov
    retval['correlation'] = cor

    retval['sigma split'] = [libmath.fit_sigma(fit_result['residuals'][s],
                                               retval['weights'][s],
                                               len(fit_result['residuals'][s]),
                                               retval['degrees of freedom'])
                             for s in retval['slices']]
    retval['sigma'] = libmath.fit_sigma(fit_result['residuals'],
                                        retval['weights'],
                                        sum(retval['number of experimental points']),
                                        retval['degrees of freedom'])

    # for compatibility in plotting
    retval['species_concentrations'] = retval['free concentration']

    return retval


def covariance_fun(J, W, F):
    """Compute covariance matrix.

    Returns the covariance matrix :math:`CV = inv(J'.W.J)*MSE`
    Where MSE is mean-square error :math:`MSE = (R'*R)/(N-p)`
    where *R* are the residuals, *N* is the number of observations and
    *p* is the number of coefficients estimated

    Parameters:
        J (:class:`numpy.ndarray`): the jacobian
        W (:class:`numpy.ndarray`): the weights matrix
        F (:class:`numpy.ndarray`): the residuals
    Returns:
        :class:`numpy.ndarray`: an (*p*, *p*)-sized array representing
            the covariance matrix.
    """
    num_params = J.shape[1] if J.ndim == 2 else 1
    mse = np.sum(F * np.diag(W) * F) / (len(F) - num_params)
    temp = np.linalg.inv(np.atleast_2d(np.dot(np.dot(J.T, W), J)))
    return temp * mse


def fit_final_calcs(jacobian, resids, weights):
    """Perform final calculations common to some routines.

    Parameters:
        jacobian (:class:`numpy.array`): the jacobian
        resids (:class:`numpy.array`): the residuals
        weights (:class:`numpy.array`): the weights
    Returns:
        * the error in beta
        * the correlation matrix
        * the covariance matrix
    """
    covariance = covariance_fun(jacobian, weights, resids)
    cov_diag = np.diag(covariance)
    stdev = np.sqrt(cov_diag)
    lenD = len(cov_diag)
    correlation = covariance / np.sqrt(
        np.dot(cov_diag.reshape((lenD, 1)), cov_diag.reshape((1, lenD)))
    )
    return stdev, correlation, covariance


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

