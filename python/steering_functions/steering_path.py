from __future__ import annotations

from .planners import (
    CC0PMDubinsStateSpace,
    CC00DubinsStateSpace,
    CC00ReedsSheppStateSpace,
    CCDubinsStateSpace,
    CCpm0DubinsStateSpace,
    CCpmpmDubinsStateSpace,
    DubinsStateSpace,
    HC0PMReedsSheppStateSpace,
    HC00ReedsSheppStateSpace,
    HCReedsSheppStateSpace,
    HCpm0ReedsSheppStateSpace,
    HCpmpmReedsSheppStateSpace,
    ReedsSheppStateSpace,
)
from .types import PathType, State


class SteeringPath:
    def __init__(
        self,
        path_type: PathType,
        kappa_max: float,
        sigma_max: float,
        discretization: float,
    ) -> None:
        if kappa_max <= 0 or sigma_max < 0 or discretization <= 0:
            raise ValueError("Invalid parameters in SteeringPath constructor")

        self._path_type = path_type
        self._discretization = discretization
        self._kappa_max = kappa_max
        self._sigma_max = sigma_max
        self._base_state_space = self._create_planner()

    def _create_planner(self):
        if self._path_type is PathType.NONE:
            raise ValueError("PathType.NONE is not a valid planner type")

        mapping = {
            PathType.CC_DUBINS: CCDubinsStateSpace,
            PathType.CC00_DUBINS: CC00DubinsStateSpace,
            PathType.CC0PM_DUBINS: CC0PMDubinsStateSpace,
            PathType.CCPM0_DUBINS: CCpm0DubinsStateSpace,
            PathType.CCPMPM_DUBINS: CCpmpmDubinsStateSpace,
            PathType.DUBINS: DubinsStateSpace,
            PathType.CC00_RS: CC00ReedsSheppStateSpace,
            PathType.HC_RS: HCReedsSheppStateSpace,
            PathType.HC00_RS: HC00ReedsSheppStateSpace,
            PathType.HC0PM_RS: HC0PMReedsSheppStateSpace,
            PathType.HCPM0_RS: HCpm0ReedsSheppStateSpace,
            PathType.HCPMPM_RS: HCpmpmReedsSheppStateSpace,
            PathType.RS: ReedsSheppStateSpace,
        }
        planner_type = mapping.get(self._path_type)
        if planner_type is None:
            raise ValueError(f"Unsupported planner type: {self._path_type}")
        return planner_type(self._discretization)

    @property
    def path_type(self) -> PathType:
        return self._path_type

    @property
    def discretization(self) -> float:
        return self._discretization

    @property
    def kappa_max(self) -> float:
        return self._kappa_max

    @property
    def sigma_max(self) -> float:
        return self._sigma_max

    def compute_shortest_control_sequence(self, start: State, goal: State):
        return self._base_state_space.get_controls(start, goal)

    def compute_shortest_path(self, start: State, goal: State):
        return self._base_state_space.get_path(start, goal)

    def compute_all_control_sequences(self, start: State, goal: State):
        return self._base_state_space.get_all_controls(start, goal)

    def compute_all_paths(self, start: State, goal: State):
        return self._base_state_space.get_all_paths(start, goal)
