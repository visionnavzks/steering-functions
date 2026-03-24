from __future__ import annotations

from enum import Enum, auto

from .core import BaseStateSpace, Control, State
from .dubins import DubinsStateSpace
from .reeds_shepp import ReedsSheppStateSpace


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


class SteeringPath:
    def __init__(self, path_type: PathType, kappa_max: float, sigma_max: float, discretization: float) -> None:
        if kappa_max <= 0.0 or sigma_max < 0.0 or discretization <= 0.0:
            raise ValueError("Invalid parameters in SteeringPath constructor")
        self._path_type = path_type
        self._discretization = discretization
        self._kappa_max = kappa_max
        self._sigma_max = sigma_max
        self._base_state_space = self._create_planner()

    def _create_planner(self) -> BaseStateSpace:
        if self._path_type is PathType.DUBINS:
            return DubinsStateSpace(self._kappa_max, self._discretization, True)
        if self._path_type is PathType.RS:
            return ReedsSheppStateSpace(self._kappa_max, self._discretization)
        raise NotImplementedError(
            f"{self._path_type.name} is not implemented in the Python port yet; "
            "supported types are DUBINS and RS."
        )

    def compute_shortest_control_sequence(self, start: State, goal: State) -> list[Control]:
        return self._base_state_space.get_controls(start, goal)

    def compute_shortest_path(self, start: State, goal: State) -> list[State]:
        return self._base_state_space.get_path(start, goal)

    def compute_all_control_sequences(self, start: State, goal: State) -> list[list[Control]]:
        return self._base_state_space.get_all_controls(start, goal)

    def compute_all_paths(self, start: State, goal: State) -> list[list[State]]:
        return self._base_state_space.get_all_paths(start, goal)

    def get_path_type(self) -> PathType:
        return self._path_type

    def get_discretization(self) -> float:
        return self._discretization

    def get_kappa_max(self) -> float:
        return self._kappa_max

    def get_sigma_max(self) -> float:
        return self._sigma_max
