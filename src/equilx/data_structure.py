from collections.abc import Sequence
from dataclasses import dataclass, field, MISSING
from enum import IntEnum, auto

import numpy as np

from types import FArray, IArray


@dataclass
class DebyeHuckelParameters:
    a: float
    b: float
    c0: float = field(init=False, default=0.0)
    c1: float = field(init=False, default=0.0)
    d0: float = field(init=False, default=0.0)
    d1: float = field(init=False, default=0.0)
    e0: float = field(init=False, default=0.0)
    e1: float = field(init=False, default=0.0)


@dataclass
class _GenericTitration:
    starting_volume: float                                         # in mL
    titre_added: FArray                                            # in mL
    starting_concentration: FArray                                 # in mol/L
    buret_concentration: FArray                                    # in mol/L
    starting_background_concentration: FArray = field(init=False, default=MISSING)  # in mol/L
    buret_background_concentration: FArray = field(init=False, default=MISSING)     # in mol/L

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
    number_of_points: int


@dataclass
class DistributionParameters(_SimulationParameters):
    starting_concentration: FArray  # in mol/L
    final_concentration: FArray     # in mol/L
    logarithmic: Sequence[bool]     # is the concentration in log units?

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


@dataclass
class ModelData:
    component_labels: list[str]
    component_charges: IArray
    stoichiometry: IArray
    log_beta: FArray
    error_log_beta: FArray
    solid_stoichiometry: IArray
    log_solids: FArray
    error_log_solids: FArray


@dataclass
class PyESData:
    mode: PyESModes
    model: ModelData

    opts_distribution: DistributionParameters | None = field(init=False, default=None)
    opts_titration: TitrationParameters | None = field(init=False, default=None)
    opts_potentiometry: PotentiometryData | None = field(init=False, default=None)

    def __post_init__(self):
        ...


@dataclass
class EquilibriumData:
    analytical_concentration = FArray
    log_beta: FArray
    stoichiometry: IArray
    log_ks: FArray = np.array([], dtype=float)
    solid_stoichiometry: IArray = np.array([[]], dtype=int)

    def __post_init__(self):
        # validate arrays
        ...

# Assertions

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
