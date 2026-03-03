from enum import IntEnum
import math
from typing import Final

class Flags(IntEnum):
    CONSTANT = 0
    REFINE = 1

LN10: Final[float] = math.log(10)
