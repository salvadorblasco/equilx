import numpy as np
from typing import TypeVar

NP = TypeVar("NP", bound=int)
NC = TypeVar("NC", bound=int)
NS = TypeVar("NS", bound=int)

# # Variable but known dimensionality
# MatrixLike[N, M] = np.ndarray[tuple[N, M], np.dtype[np.float64]]
# VectorLike[N]    = np.ndarray[tuple[N],         np.dtype[np.float64]]
# 
# # Fully generic shape (any number of dimensions)
# NDArrayFloat = np.ndarray[tuple[int, ...], np.dtype[np.float64]]
