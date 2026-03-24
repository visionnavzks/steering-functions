from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto


class PathType(Enum):
    NONE = auto()
    CC_DUBINS = auto()
    CC00_DUBINS = auto()
    CC0PM_DUBINS = auto()
    CCPM0_DUBINS = auto()
    CCPMPM_DUBINS = auto()
    DUBINS = auto()
    CC00_RS = auto()
    HC_RS = auto()
    HC00_RS = auto()
    HC0PM_RS = auto()
    HCPM0_RS = auto()
    HCPMPM_RS = auto()
    RS = auto()


@dataclass
class State:
    x: float = 0.0
    y: float = 0.0
    theta: float = 0.0
    kappa: float = 0.0
    sigma: float = 0.0
    d: float = 0.0
    s: float = 0.0
    vel: float = 0.0
    acc: float = 0.0
    time: float = 0.0
    fork_y: float = 0.0


@dataclass
class Control:
    delta_s: float
    kappa: float
    sigma: float
