"""High-level SteeringPath wrapper and PathType enum (Python port)."""

from enum import Enum, auto
from steering_functions.state import State, Control
from steering_functions.base_state_space import BaseStateSpace


class PathType(Enum):
    """Supported path planning algorithm types."""
    NONE = 0
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


def _create_planner(path_type, kappa_max, sigma_max, discretization):
    """Factory: create the appropriate BaseStateSpace subclass."""

    if path_type == PathType.DUBINS:
        from steering_functions.dubins_state_space import DubinsStateSpace
        return DubinsStateSpace(kappa_max, discretization, True)

    if path_type == PathType.RS:
        from steering_functions.reeds_shepp_state_space import ReedsSheppStateSpace
        return ReedsSheppStateSpace(kappa_max, discretization)

    if path_type == PathType.CC_DUBINS:
        from steering_functions.cc_dubins_state_space import CC_Dubins_State_Space
        return CC_Dubins_State_Space(kappa_max, sigma_max, discretization, True)

    if path_type == PathType.CC00_DUBINS:
        from steering_functions.cc_dubins_state_space import CC00_Dubins_State_Space
        return CC00_Dubins_State_Space(kappa_max, sigma_max, discretization, True)

    if path_type == PathType.CC0PM_DUBINS:
        from steering_functions.cc_dubins_state_space import CC0pm_Dubins_State_Space
        return CC0pm_Dubins_State_Space(kappa_max, sigma_max, discretization, True)

    if path_type == PathType.CCPM0_DUBINS:
        from steering_functions.cc_dubins_state_space import CCpm0_Dubins_State_Space
        return CCpm0_Dubins_State_Space(kappa_max, sigma_max, discretization, True)

    if path_type == PathType.CCPMPM_DUBINS:
        from steering_functions.cc_dubins_state_space import CCpmpm_Dubins_State_Space
        return CCpmpm_Dubins_State_Space(kappa_max, sigma_max, discretization, True)

    if path_type == PathType.CC00_RS:
        from steering_functions.cc_reeds_shepp_state_space import CC00_Reeds_Shepp_State_Space
        return CC00_Reeds_Shepp_State_Space(kappa_max, sigma_max, discretization)

    if path_type == PathType.HC_RS:
        from steering_functions.hc_reeds_shepp_state_space import HC_Reeds_Shepp_State_Space
        return HC_Reeds_Shepp_State_Space(kappa_max, sigma_max, discretization)

    if path_type == PathType.HC00_RS:
        from steering_functions.hc_reeds_shepp_state_space import HC00_Reeds_Shepp_State_Space
        return HC00_Reeds_Shepp_State_Space(kappa_max, sigma_max, discretization)

    if path_type == PathType.HC0PM_RS:
        from steering_functions.hc_reeds_shepp_state_space import HC0pm_Reeds_Shepp_State_Space
        return HC0pm_Reeds_Shepp_State_Space(kappa_max, sigma_max, discretization)

    if path_type == PathType.HCPM0_RS:
        from steering_functions.hc_reeds_shepp_state_space import HCpm0_Reeds_Shepp_State_Space
        return HCpm0_Reeds_Shepp_State_Space(kappa_max, sigma_max, discretization)

    if path_type == PathType.HCPMPM_RS:
        from steering_functions.hc_reeds_shepp_state_space import HCpmpm_Reeds_Shepp_State_Space
        return HCpmpm_Reeds_Shepp_State_Space(kappa_max, sigma_max, discretization)

    # Default fallback
    from steering_functions.reeds_shepp_state_space import ReedsSheppStateSpace
    return ReedsSheppStateSpace(kappa_max, discretization)


class SteeringPath:
    """High-level wrapper for steering path computation.

    Provides a unified interface to compute optimal (shortest) or all possible
    paths between two states using the specified steering function algorithm.
    """

    def __init__(self, path_type, kappa_max, sigma_max, discretization):
        if kappa_max <= 0 or sigma_max < 0 or discretization <= 0:
            raise ValueError("Invalid parameters in SteeringPath constructor")

        self._path_type = path_type
        self._discretization = discretization
        self._kappa_max = kappa_max
        self._sigma_max = sigma_max
        self._planner = _create_planner(path_type, kappa_max, sigma_max, discretization)
        if self._planner is None:
            raise RuntimeError("Failed to create path planner")

    # --- Optimal (shortest) path interface ---

    def compute_shortest_control_sequence(self, start, goal):
        """Get the control sequence for the shortest path."""
        return self._planner.get_controls(start, goal)

    def compute_shortest_path(self, start, goal):
        """Get the discretized shortest path."""
        return self._planner.get_path(start, goal)

    # --- All paths interface ---

    def compute_all_control_sequences(self, start, goal):
        """Get all possible control sequences between two states."""
        return self._planner.get_all_controls(start, goal)

    def compute_all_paths(self, start, goal):
        """Get all possible discretized paths between two states."""
        return self._planner.get_all_paths(start, goal)

    # --- Getters ---

    @property
    def path_type(self):
        return self._path_type

    @property
    def discretization(self):
        return self._discretization

    @property
    def kappa_max(self):
        return self._kappa_max

    @property
    def sigma_max(self):
        return self._sigma_max
