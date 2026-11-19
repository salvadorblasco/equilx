import itertools
from collections.abc import Sequence
from dataclasses import dataclass, field
from functools import cached_property
from enum import IntEnum, auto
from collections.abc import Iterable
from collections.abc import Iterator

import numpy as np
from numpy.typing import NDArray

from equilx.mytypes import BArray, FArray, IArray
from equilx.consts import Flags


@dataclass
class DebyeHuckelParameters:
    a: float = 0.5115
    b: float = 1.50
    c0: float | FArray = 0.1
    c1: float | FArray = 0.209
    d0: float | FArray = 0.0
    d1: float | FArray = 0.093
    e0: float | FArray = 0.0
    e1: float | FArray = 0.0


@dataclass(kw_only=True)
class _GenericTitration:
    starting_volume: float                                 # in mL
    titre_added: FArray                                    # in mL
    initial_amount: FArray                                 # in mmol
    buret_concentration: FArray                            # in mol/L
    volume_error: float = 0.0
    use_point: BArray | None = None
    background_ions: FArray = field(default_factory=lambda: np.array([], dtype=float))

    def __post_init__(self):
        # type checking
        _assert_positivefloat(self.starting_volume, "starting volume")
        _assert_floatarray(self.initial_amount, "initial amount", ndim=1)

        number_components = len(self.initial_amount)
        _assert_floatarray(self.buret_concentration, "buret concentration", ndim=1)
        if len(self.buret_concentration) != number_components:
            raise ValueError("starting concentration and buret concentration have different length")

        if 'flags_initial' not in dir(self):
            self.flags_initial = np.zeros_like(self.initial_amount)
        else:
            _assert_floatarray(self.s, "flags_initial", ndim=1)

        if 'flags_buret' not in dir(self):
            self.flags_buret = np.zeros_like(self.initial_amount)
        else:
            _assert_floatarray(self.flags_buret, "flags_buret", ndim=1)

        if 'initial_background_concentration' not in dir(self):
            self.initial_background_concentration = np.zeros_like(self.initial_amount)
        else:
            _assert_floatarray(self.initial_background_concentration, "initial_background_concentration", ndim=1)

        if 'buret_background_concentration' not in dir(self):
            self.buret_background_concentration = np.zeros_like(self.initial_amount)
        else:
            _assert_floatarray(self.buret_background_concentration, "buret_background_concentration", ndim=1)


    def get_dof(self) -> int:
        dof: int = int(np.count_nonzero(self.flags_initial == Flags.REFINE))
        dof+= int(np.count_nonzero(self.flags_buret == Flags.REFINE))
        return dof

    def get_titre(self) -> FArray:
        return np.extract(self.use_point, self.titre_added)

    @property
    def analytical_concentration(self) -> FArray:
        return self.__build_titration(self.initial_amount, self.buret_concentration)

    @property
    def ionic_strength(self) -> FArray:
        return self.__build_titration(self.starting_background_concentration,
                                      self.buret_background_concentration)

    @property
    def number_of_points(self) -> int:
        return len(self)

    def __build_titration(self, initialc: FArray, buret: FArray) -> FArray:
        return (initialc[None, :]*self.starting_volume + buret[None, :]*self.titre_added[:, None]) \
              / (self.starting_volume + self.titre_added[:, None])

    def __len__(self) -> int:
        return len(self.titre_added)


@dataclass
class _SimulationParameters:
    reference_component: int        # index of the independent component


@dataclass
class DistributionParameters(_SimulationParameters):
    starting_concentration: FArray  # in mol/L
    final_concentration: FArray     # in mol/L
    logarithmic: Sequence[bool]     # is the concentration in log units?
    number_of_points: int

    def __post_init__(self):
        _assert_floatarray(self.starting_concentration, "starting concentration", ndim=1)
        _assert_floatarray(self.final_concentration, "final concentration", ndim=1)
        _assert_same_length(self.final_concentration,
                             self.starting_concentration,
                             "final concentration")
        if not isinstance(self.logarithmic, Sequence):
            raise TypeError("logarithmic must be a sequence")
        _assert_same_length(self.logarithmic, self.starting_concentration, "logarithmic")

    def check(self) -> bool:
        "Return whether the data is complete and ready."
        return True

    @property
    def analytical_concentration(self) -> FArray:
        retval = np.empty((self.number_of_points, self.number_components), dtype=float)
        for icomp, (start, stop, islog) in enumerate(self.__auxiter()):
            if islog:
                func = np.logspace
            else:
                func = np.linspace
            retval[:, icomp] = func(start, stop, num=self.number_of_points)
        return retval

    @property
    def number_components(self) -> int:
        return len(self.starting_concentration)

    def __auxiter(self):
        yield from zip(self.starting_concentration, self.final_concentration, self.logarithmic)


@dataclass
class TitrationParameters(_SimulationParameters, _GenericTitration):
    pass


@dataclass
class PotentiometryTitration(_GenericTitration):
    emf: FArray
    emf0: float
    emf_error: float
    slope: float
    electroactive: int
    temperature: float = field(init=False, default=25.0)  # Celsius
    junction_potential: list[float] = field(default_factory=list)
    junction_component: list[int] = field(default_factory=list)
    flag_junction: list[int] = field(default_factory=list)

    flag_emf0: int = Flags.CONSTANT
    flag_slope: int = Flags.CONSTANT

    def __post_init__(self):
        super().__post_init__()
        _assert_floatarray(self.emf, "emf", ndim=1)
        _assert_float(self.emf0, "emf0")
        _assert_float(self.emf_error, "emf error")
        _assert_float(self.slope, "slope")
        _assert_same_length(self.titre_added, self.emf, "emf")

    def get_emf(self) -> FArray:
        return np.extract(self.use_point, self.emf)


@dataclass
class PotentiometryData:
    titrations: Sequence[PotentiometryTitration]

    def __post_init__(self):
        if not all(isinstance(t, PotentiometryTitration) for t in self.titrations):
            raise TypeError("all instances of titration must be PotentiometryTitration")

    def get_emf0(self) -> Iterator[float]:
        for tit in self.titrations:
            yield tit.emf0

    def get_slopes(self) -> Iterator[float]:
        for tit in self.titrations:
            yield tit.slope

    def get_electroactive(self) -> Iterator[int]:
        for tit in self.titrations:
            yield tit.electroactive

    def get_per_point_electroactive(self) -> Iterator[int]:
        for hindex, tit in zip(self.get_electroactive(), self.titrations):
            yield from itertools.repeat(hindex, tit.number_of_points)

    def get_collapsed_emf(self) -> FArray:
        return np.concatenate([t.get_emf() for t in self.titrations])

    @property
    def number_of_titrations(self) -> int:
        return len(self.titrations)

    @property
    def total_experimental_points(self) -> int:
        return sum(len(t) for t in self.titrations)

    def get_dof(self) -> int:
        return sum(t.get_dof() for t in self.titrations)

    def get_slices(self) -> list[slice]:
        slices: list[slice] = []
        counter: int = 0
        for tit in self.titrations:
            slices.append(slice(counter, counter+tit.number_of_points))
            counter += tit.number_of_points
        return slices

    def __iter__(self):
        yield from self.titrations


class PyESModes(IntEnum):
    DISTRIBUTION = auto()
    TITRATION = auto()
    POTENTIOMETRY = auto()


@dataclass(kw_only=True)
class _EquilCommon:
    log_beta: FArray
    stoichiometry: IArray
    error_log_beta: FArray | None = None
    solid_stoichiometry: IArray = field(default_factory=lambda: np.array([[]], dtype=int))
    log_ks: FArray = field(default_factory=lambda : np.array([], dtype=float))
    error_log_ks: FArray | None = None
    reference_ionic: float = 0.0
    charges: IArray | None = None
    component_labels: Sequence[str] | None = None

    def __post_init__(self):
        if self.charges is None:
            self.charges = np.zeros(self.n_components, dtype=int)
        if self.component_labels is None:
            self.labels = self.n_components * [""]
        if self.error_log_beta is None:
            self.error_log_beta = np.zeros_like(self.log_beta)
        if self.error_log_ks is None:
            self.error_log_ks = np.zeros_like(self.log_ks)

    def any_solid(self):
        return (self.n_solids > 0)

    @cached_property
    def n_components(self) -> int:
        return self.stoichiometry.shape[1]

    @cached_property
    def n_equils(self) -> int:
        return self.stoichiometry.shape[0]

    @cached_property
    def n_solids(self) -> int:
        if self.log_ks is None:
            return 0
        return len(self.log_ks)

    @cached_property
    def n_soluble(self) -> int:
        return self.n_components + self.n_equils

    @cached_property
    def stoichiometryx(self):
        return np.concatenate((np.eye(self.n_components, dtype=int),self.stoichiometry), axis=0)

    @cached_property
    def species_names(self) -> list[str]:
        return list(_species_names(self))

    @cached_property
    def chargesx(self):
        if self.charges is None:
            return None
        species_charge = self.stoichiometry @ self.charges[:,None]
        return np.concatenate((self.charges, species_charge.flatten()))

    @cached_property
    def zstar(self) -> IArray|None:
        if self.charges is None:
            return None
        part1 = self.stoichiometry @ self.charges**2
        part2 = self.stoichiometry @ self.charges
        return part1 + part2**2

    @cached_property
    def pstar(self) -> IArray|None:
        if self.charges is None:
            return None
        return np.sum(self.stoichiometry, axis=1) - 1


@dataclass(kw_only=True)
class EquilibriumData(_EquilCommon):
    analytical_concentration: FArray
    background_ions: FArray = field(default_factory=lambda: np.array([], dtype=float))

    @property
    def n_points(self) -> int:
        return self.analytical_concentration.shape[0]

    def __copy__(self):
        cls = self.__class__
        copy = cls.__new__(cls)
        copy.__dict__.update(self.__dict__)
        copy.__dict__['analytical_concentration'] = self.analytical_concentration.copy()
        return copy


@dataclass(kw_only=True)
class ModelData(_EquilCommon):
    log_beta_flags: list[Flags]
    log_solids_flags: list[Flags] = field(init=False, default_factory=list)

    @property
    def degrees_of_freedom(self) -> int:
        count1 = sum(1 for flag in self.log_beta_flags if flag == Flags.REFINE)
        count2 = sum(1 for flag in self.log_solids_flags if flag == Flags.REFINE)
        return count1 + count2

    def get_equildata(self, analytical_concentration: FArray) -> EquilibriumData:
        return EquilibriumData(
            analytical_concentration=analytical_concentration,
            log_beta=self.log_beta,
            stoichiometry=self.stoichiometry,
            log_ks=self.log_ks,
            solid_stoichiometry=self.solid_stoichiometry,
            charges=self.charges,
            reference_ionic=self.reference_ionic,
            background_ions=self.background_ions)

    def get_refinable(self) -> Iterator[int]:
        yield from map(lambda x: x == Flags.REFINE, self.log_beta_flags)


# Utilities functions

def _species_names(equil: _EquilCommon, replace_proton=True) -> Iterable[str]:
    if equil.component_labels is None:
        raise ValueError("component labels must be provided")
    for stoich in equil.stoichiometry:
        aux = []
        for label, s in zip(equil.component_labels, stoich):
            if label == 'H' and s < 0 and replace_proton:
                aux.append(f'(OH){-s}')
            elif s == 0:
                continue
            else:
                aux.append(f'({label}){s}')
        yield "".join(aux)


# Assertions

def _assert_dim(array: NDArray, dim:int):
    if array.ndim != dim:
        raise TypeError(f"Array must be {dim}-dimmensioned")


def _assert_floatarray(array, description: str, ndim:int = 1):
    if not isinstance(array, np.ndarray):
        raise TypeError(f"{description} must be an array")
    if not array.dtype.kind == "f":
        raise TypeError(f"{description} must be a float array")
    if not array.ndim == ndim:
        raise TypeError(f"{description} must have {ndim} dimension(s)")


def _assert_float(number, description: str):
    if not isinstance(number, float):
        raise TypeError(f"{description} must be a float")


def _assert_positivefloat(number, description: str):
    _assert_float(number, description)
    if number < 0.0:
        raise ValueError("{description} must no be negative")


def _assert_same_length(seq1, seq2, description:str):
    if len(seq1) != len(seq2):
        raise ValueError(f"{description} has wrong number of elements")
