import numpy as np
import numpy.typing as npt
from typing import TypeVar

NP = TypeVar("NP", bound=int)
NC = TypeVar("NC", bound=int)
NS = TypeVar("NS", bound=int)
NE = TypeVar("NE", bound=int)

FArray = npt.NDArray[np.float64]
IArray = npt.NDArray[np.int64]
BArray = npt.NDArray[np.bool_]
