"""
Dubins state space – shortest paths with bounded curvature.

Ported from C++ steering_functions/src/dubins_state_space/dubins_state_space.cpp.
"""

import math
from dataclasses import dataclass, field
from typing import List, Optional

from steering_functions.state import State, Control
from steering_functions.utilities import TWO_PI, twopify
from steering_functions.base_state_space import BaseStateSpace

# ---------------------------------------------------------------------------
# Segment types
# ---------------------------------------------------------------------------
DUBINS_LEFT = 0
DUBINS_STRAIGHT = 1
DUBINS_RIGHT = 2

# Six canonical Dubins path types (LSL, RSR, RSL, LSR, RLR, LRL)
DUBINS_PATH_TYPE = [
    [DUBINS_LEFT, DUBINS_STRAIGHT, DUBINS_LEFT],
    [DUBINS_RIGHT, DUBINS_STRAIGHT, DUBINS_RIGHT],
    [DUBINS_RIGHT, DUBINS_STRAIGHT, DUBINS_LEFT],
    [DUBINS_LEFT, DUBINS_STRAIGHT, DUBINS_RIGHT],
    [DUBINS_RIGHT, DUBINS_LEFT, DUBINS_RIGHT],
    [DUBINS_LEFT, DUBINS_RIGHT, DUBINS_LEFT],
]

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------
_DUBINS_EPS = 1e-6
_DUBINS_ZERO = -1e-9


# ---------------------------------------------------------------------------
# Dubins path representation
# ---------------------------------------------------------------------------
@dataclass
class DubinsPath:
    """Compact representation of a Dubins path (three segments)."""

    type_: List[int] = field(default_factory=lambda: list(DUBINS_PATH_TYPE[0]))
    length_: List[float] = field(default_factory=lambda: [0.0, float("inf"), 0.0])

    def length(self) -> float:
        return self.length_[0] + self.length_[1] + self.length_[2]


# ---------------------------------------------------------------------------
# Six primitive Dubins word functions
# ---------------------------------------------------------------------------

def _dubins_LSL(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sa - sb))
    if tmp >= _DUBINS_ZERO:
        theta = math.atan2(cb - ca, d + sa - sb)
        t = twopify(-alpha + theta)
        p = math.sqrt(max(tmp, 0.0))
        q = twopify(beta - theta)
        return DubinsPath(list(DUBINS_PATH_TYPE[0]), [t, p, q])
    return DubinsPath()


def _dubins_RSR(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sb - sa))
    if tmp >= _DUBINS_ZERO:
        theta = math.atan2(ca - cb, d - sa + sb)
        t = twopify(alpha - theta)
        p = math.sqrt(max(tmp, 0.0))
        q = twopify(-beta + theta)
        return DubinsPath(list(DUBINS_PATH_TYPE[1]), [t, p, q])
    return DubinsPath()


def _dubins_RSL(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = d * d - 2.0 + 2.0 * (ca * cb + sa * sb - d * (sa + sb))
    if tmp >= _DUBINS_ZERO:
        p = math.sqrt(max(tmp, 0.0))
        theta = math.atan2(ca + cb, d - sa - sb) - math.atan2(2.0, p)
        t = twopify(alpha - theta)
        q = twopify(beta - theta)
        return DubinsPath(list(DUBINS_PATH_TYPE[2]), [t, p, q])
    return DubinsPath()


def _dubins_LSR(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = -2.0 + d * d + 2.0 * (ca * cb + sa * sb + d * (sa + sb))
    if tmp >= _DUBINS_ZERO:
        p = math.sqrt(max(tmp, 0.0))
        theta = math.atan2(-ca - cb, d + sa + sb) - math.atan2(-2.0, p)
        t = twopify(-alpha + theta)
        q = twopify(-beta + theta)
        return DubinsPath(list(DUBINS_PATH_TYPE[3]), [t, p, q])
    return DubinsPath()


def _dubins_RLR(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = 0.125 * (6.0 - d * d + 2.0 * (ca * cb + sa * sb + d * (sa - sb)))
    if abs(tmp) < 1.0:
        p = TWO_PI - math.acos(tmp)
        theta = math.atan2(ca - cb, d - sa + sb)
        t = twopify(alpha - theta + 0.5 * p)
        q = twopify(alpha - beta - t + p)
        return DubinsPath(list(DUBINS_PATH_TYPE[4]), [t, p, q])
    return DubinsPath()


def _dubins_LRL(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = 0.125 * (6.0 - d * d + 2.0 * (ca * cb + sa * sb - d * (sa - sb)))
    if abs(tmp) < 1.0:
        p = TWO_PI - math.acos(tmp)
        theta = math.atan2(-ca + cb, d + sa - sb)
        t = twopify(-alpha + theta + 0.5 * p)
        q = twopify(beta - alpha - t + p)
        return DubinsPath(list(DUBINS_PATH_TYPE[5]), [t, p, q])
    return DubinsPath()


def _dubins_shortest(d: float, alpha: float, beta: float) -> DubinsPath:
    """Select the shortest Dubins path among all six types."""
    if d < _DUBINS_EPS and abs(alpha - beta) < _DUBINS_EPS:
        return DubinsPath(list(DUBINS_PATH_TYPE[0]), [0.0, d, 0.0])

    path = _dubins_LSL(d, alpha, beta)
    min_length = path.length()

    tmp = _dubins_RSR(d, alpha, beta)
    length = tmp.length()
    if length < min_length:
        min_length = length
        path = tmp

    tmp = _dubins_RSL(d, alpha, beta)
    length = tmp.length()
    if length < min_length:
        min_length = length
        path = tmp

    tmp = _dubins_LSR(d, alpha, beta)
    length = tmp.length()
    if length < min_length:
        min_length = length
        path = tmp

    tmp = _dubins_RLR(d, alpha, beta)
    length = tmp.length()
    if length < min_length:
        min_length = length
        path = tmp

    tmp = _dubins_LRL(d, alpha, beta)
    if tmp.length() < min_length:
        path = tmp

    return path


# ---------------------------------------------------------------------------
# Dubins state space
# ---------------------------------------------------------------------------
class DubinsStateSpace(BaseStateSpace):
    """
    Dubins state space – shortest paths composed of arcs and straight lines
    with bounded curvature *kappa*.
    """

    def __init__(self, kappa: float, discretization: float = 0.1, forwards: bool = True):
        super().__init__(discretization)
        self._kappa = kappa
        self._kappa_inv = 1.0 / kappa
        self._forwards = forwards

    # ------------------------------------------------------------------
    # Core Dubins computation (normalised)
    # ------------------------------------------------------------------
    def dubins(self, state1: State, state2: State) -> DubinsPath:
        """Compute shortest Dubins path between two states."""
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        th = math.atan2(dy, dx)
        d = math.sqrt(dx * dx + dy * dy) * self._kappa
        alpha = twopify(state1.theta - th)
        beta = twopify(state2.theta - th)
        return _dubins_shortest(d, alpha, beta)

    # ------------------------------------------------------------------
    # Distance
    # ------------------------------------------------------------------
    def get_distance(self, state1: State, state2: State) -> float:
        """Return shortest path length from *state1* to *state2*."""
        if self._forwards:
            return self._kappa_inv * self.dubins(state1, state2).length()
        else:
            return self._kappa_inv * self.dubins(state2, state1).length()

    # ------------------------------------------------------------------
    # Controls (core helper)
    # ------------------------------------------------------------------
    def _get_controls_core(
        self, state1: State, state2: State, reverse_direction: bool = False
    ) -> List[Control]:
        if not reverse_direction:
            path = self.dubins(state1, state2)
        else:
            path = self.dubins(state2, state1)

        dubins_controls: List[Control] = []
        for i in range(3):
            control = Control()
            seg = path.type_[i]
            if seg == DUBINS_LEFT:
                control.delta_s = self._kappa_inv * path.length_[i]
                control.kappa = self._kappa
                control.sigma = 0.0
            elif seg == DUBINS_STRAIGHT:
                control.delta_s = self._kappa_inv * path.length_[i]
                control.kappa = 0.0
                control.sigma = 0.0
            elif seg == DUBINS_RIGHT:
                control.delta_s = self._kappa_inv * path.length_[i]
                control.kappa = -self._kappa
                control.sigma = 0.0
            if abs(control.delta_s) > _DUBINS_EPS:
                dubins_controls.append(control)

        if reverse_direction:
            dubins_controls.reverse()
            for c in dubins_controls:
                c.delta_s = -c.delta_s

        return dubins_controls

    # ------------------------------------------------------------------
    # Public controls interface
    # ------------------------------------------------------------------
    def get_controls(self, state1: State, state2: State) -> List[Control]:
        """Return controls for the shortest forward path."""
        return self._get_controls_core(state1, state2, False)

    def get_controls_reverse(self, state1: State, state2: State) -> List[Control]:
        """Return controls for a pure reverse driving path."""
        return self._get_controls_core(state1, state2, True)

    def get_all_controls(self, state1: State, state2: State) -> List[List[Control]]:
        """Return both forward and reverse control sequences."""
        return [
            self.get_controls(state1, state2),
            self.get_controls_reverse(state1, state2),
        ]