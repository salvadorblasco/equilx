from collections.abc import Sequence
from dataclasses import dataclass, field, MISSING
from functools import cached_property
from enum import IntEnum, auto
from typing import Optional

import numpy as np

from .mytypes import FArray, IArray

from numpy.typing import NDArray

from equilx.consts import Flags


@dataclass
class DebyeHuckelParameters:
    a: float = 0.5115
    b: float = 1.50
    c0: float = 0.1
    c1: float = 0.209
    d0: float = 0.0
    d1: float = 0.093
    e0: float = 0.0
    e1: float = 0.0


@dataclass
class _GenericTitration:
    starting_volume: float                                         # in mL
    titre_added: FArray                                            # in mL
    starting_concentration: FArray                                 # in mol/L
    buret_concentration: FArray                                    # in mol/L
    starting_background_concentration: FArray = field(init=False, default_factory=MISSING)  # in mol/L
    buret_background_concentration: FArray = field(init=False, default_factory=MISSING)     # in mol/L

    def __post_init__(self):
        # type checking
        __assert_positivefloat(self.starting_volume, "starting volume")
        __assert_floatarray(self.starting_concentration, "starting concentration", ndim=1)

        number_components = len(self.starting_concentration)
        __assert_floatarray(self.buret_concentration, "buret concentration", ndim=1)
        if len(self.buret_concentration) != number_components:
            raise ValueError("starting concentration and buret concentration have different length")

        self.__background_init("starting_background_concentration")
        self.__background_init("buret_background_concentration")

    def __background_init(self, what):
        array = getattr(self, what)
        if array is MISSING:
            setattr(self, what, np.zeros_like(self.starting_concentration))
        else:
            __assert_floatarray(array, what, ndim=1)

    @property
    def analytical_concentration(self) -> FArray:
        return self.__build_titration(self.starting_concentration, self.buret_concentration)

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
        __assert_floatarray(self.starting_concentration, "starting concentration", ndim=1)
        __assert_floatarray(self.final_concentration, "final concentration", ndim=1)
        __assert_same_length(self.final_concentration,
                             self.starting_concentration,
                             "final concentration")
        if not isinstance(self.logarithmic, Sequence):
            raise TypeError("logarithmic must be a sequence")
        __assert_same_length(self.logarithmic, self.starting_concentration, "logarithmic")

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
    temperature: float = field(init=False, default=298.15)  # kelvin
    starting_concentration_flag: list[Flags]
    buret_concentration_flag: list[Flags]

    def __post_init__(self):
        super().__post_init__()
        __assert_floatarray(self.emf, "emf", ndim=1)
        __assert_float(self.emf0, "emf0")
        __assert_float(self.emf_error, "emf error")
        __assert_float(self.slope, "slope")
        __assert_same_length(self.titre_added, self.emf, "emf")


@dataclass
class PotentiometryData:
    titrations: Sequence[PotentiometryTitration]


class PyESModes(IntEnum):
    DISTRIBUTION = auto()
    TITRATION = auto()
    POTENTIOMETRY = auto()


@dataclass(kw_only=True)
class _EquilCommon:
    log_beta: FArray
    stoichiometry: IArray
    log_ks: Optional[FArray] = field(default_factory=lambda : np.array([], dtype=float))
    solid_stoichiometry: Optional[IArray] = field(default_factory=lambda: np.array([[]], dtype=int))
    error_log_beta: Optional[FArray] = None
    error_log_ks: Optional[FArray] = None
    reference_ionic: float = 0.0
    charges: Optional[IArray] = None
    component_labels: Optional[list[str]] = None

    def __post_init__(self):
        if self.charges is None:
            self.charges = np.zeros(self.n_components, dtype=int)
        if self.component_labels is None:
            self.labels = self.n_components * [""]
        if self.error_log_beta is None:
            self.error_log_beta = np.zeros_like(self.log_beta)
        if self.error_log_ks is None:
            self.error_log_ks = np.zeros_like(self.log_ks)

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
class ModelData(_EquilCommon):
    log_beta_flags: list[Flags]
    log_solids_flags: list[Flags]


@dataclass(kw_only=True)
class EquilibriumData(_EquilCommon):
    analytical_concentration: FArray
    background_ions: FArray = field(default_factory=lambda :np.array([], dtype=float))

    @cached_property
    def n_points(self) -> int:
        return self.analytical_concentration.shape[0]



# Assertions

def _assert_dim(array: NDArray, dim:int):
    if array.ndim != dim:
        raise TypeError(f"Array must be {dim}-dimmensioned")


def __assert_floatarray(array, description: str, ndim:int = 1):
    if not isinstance(array, np.ndarray):
        raise TypeError(f"{description} must be an array")
    if not array.dtype.kind != "f":
        raise TypeError(f"{description} must be a float array")
    if not array.ndim == ndim:
        raise TypeError(f"{description} must have {ndim} dimension(s)")


def __assert_float(number, description: str):
    if not isinstance(number, float):
        raise TypeError(f"{description} must be a float")


def __assert_positivefloat(number, description: str):
    __assert_float(number, description)
    if number < 0.0:
        raise ValueError("{description} must no be negative")


def __assert_same_length(seq1, seq2, description:str):
    if len(seq1) != len(seq2):
        raise ValueError(f"{description} has wrong number of elements")
