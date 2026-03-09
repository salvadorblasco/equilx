from typing import Protocol, Any, Callable

from equilx.data_structure import ExperimentalData, ModelData
from equilx.mytypes import FArray


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

        self._freeconcentration: FArray | None = None
