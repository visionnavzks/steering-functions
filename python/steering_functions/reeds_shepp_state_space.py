# Copyright (c) 2017 - for information on the respective copyright owner
# see the NOTICE file and/or the repository
#     https://github.com/hbanzhaf/steering_functions.git
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Derived from the Open Motion Planning Library (OMPL) V 1.3.1
# Copyright (c) 2010, Rice University, licensed under the BSD license.
# Original algorithm from:
#   J.A. Reeds and L.A. Shepp, "Optimal paths for a car that goes both
#   forwards and backwards," Pacific Journal of Mathematics,
#   145(2):367-393, 1990.

"""Reeds-Shepp state space: shortest paths for a car that can go forwards and backwards."""

import math
import sys
from enum import IntEnum
from typing import List, Optional, Tuple

from .base_state_space import BaseStateSpace
from .state import Control, State
from .utilities import PI, pify, polar

_RS_EPS = 1e-6
_RS_ZERO = 10.0 * sys.float_info.epsilon


class _SegType(IntEnum):
    NOP = 0
    LEFT = 1
    STRAIGHT = 2
    RIGHT = 3


# 18 Reeds-Shepp path types (5 segments each; RS_NOP pads unused slots)
_PATH_TYPES = [
    (_SegType.LEFT, _SegType.RIGHT, _SegType.LEFT, _SegType.NOP, _SegType.NOP),        # 0
    (_SegType.RIGHT, _SegType.LEFT, _SegType.RIGHT, _SegType.NOP, _SegType.NOP),       # 1
    (_SegType.LEFT, _SegType.RIGHT, _SegType.LEFT, _SegType.RIGHT, _SegType.NOP),      # 2
    (_SegType.RIGHT, _SegType.LEFT, _SegType.RIGHT, _SegType.LEFT, _SegType.NOP),      # 3
    (_SegType.LEFT, _SegType.RIGHT, _SegType.STRAIGHT, _SegType.LEFT, _SegType.NOP),   # 4
    (_SegType.RIGHT, _SegType.LEFT, _SegType.STRAIGHT, _SegType.RIGHT, _SegType.NOP),  # 5
    (_SegType.LEFT, _SegType.STRAIGHT, _SegType.RIGHT, _SegType.LEFT, _SegType.NOP),   # 6
    (_SegType.RIGHT, _SegType.STRAIGHT, _SegType.LEFT, _SegType.RIGHT, _SegType.NOP),  # 7
    (_SegType.LEFT, _SegType.RIGHT, _SegType.STRAIGHT, _SegType.RIGHT, _SegType.NOP),  # 8
    (_SegType.RIGHT, _SegType.LEFT, _SegType.STRAIGHT, _SegType.LEFT, _SegType.NOP),   # 9
    (_SegType.RIGHT, _SegType.STRAIGHT, _SegType.RIGHT, _SegType.LEFT, _SegType.NOP),  # 10
    (_SegType.LEFT, _SegType.STRAIGHT, _SegType.LEFT, _SegType.RIGHT, _SegType.NOP),   # 11
    (_SegType.LEFT, _SegType.STRAIGHT, _SegType.RIGHT, _SegType.NOP, _SegType.NOP),    # 12
    (_SegType.RIGHT, _SegType.STRAIGHT, _SegType.LEFT, _SegType.NOP, _SegType.NOP),    # 13
    (_SegType.LEFT, _SegType.STRAIGHT, _SegType.LEFT, _SegType.NOP, _SegType.NOP),     # 14
    (_SegType.RIGHT, _SegType.STRAIGHT, _SegType.RIGHT, _SegType.NOP, _SegType.NOP),   # 15
    (_SegType.LEFT, _SegType.RIGHT, _SegType.STRAIGHT, _SegType.LEFT, _SegType.RIGHT), # 16
    (_SegType.RIGHT, _SegType.LEFT, _SegType.STRAIGHT, _SegType.RIGHT, _SegType.LEFT), # 17
]


class ReedsSheppPath:
    """Complete description of a Reeds-Shepp path (up to five segments)."""

    def __init__(
        self,
        path_type: tuple = _PATH_TYPES[0],
        t: float = sys.float_info.max,
        u: float = 0.0,
        v: float = 0.0,
        w: float = 0.0,
        x: float = 0.0,
    ) -> None:
        self.type = path_type
        self.length_ = [t, u, v, w, x]
        self.total_length_ = abs(t) + abs(u) + abs(v) + abs(w) + abs(x)

    def length(self) -> float:
        return self.total_length_

    def get_length(self, index: int) -> float:
        return self.length_[index]


# ---------------------------------------------------------------------------
# Internal helper functions (mirror of the C++ anonymous namespace)
# ---------------------------------------------------------------------------


def _tau_omega(
    u: float, v: float, xi: float, eta: float, phi: float
) -> Tuple[float, float]:
    delta = pify(u - v)
    A = math.sin(u) - math.sin(delta)
    B = math.cos(u) - math.cos(delta) - 1.0
    t1 = math.atan2(eta * A - xi * B, xi * A + eta * B)
    t2 = 2.0 * (math.cos(delta) - math.cos(v) - math.cos(u)) + 3.0
    tau = pify(t1 + PI) if t2 < 0 else pify(t1)
    omega = pify(tau - u + v - phi)
    return tau, omega


def _LpSpLp(
    x: float, y: float, phi: float
) -> Optional[Tuple[float, float, float]]:
    """Formula 8.1."""
    u, t = polar(x - math.sin(phi), y - 1.0 + math.cos(phi))
    if t >= -_RS_ZERO:
        v = pify(phi - t)
        if v >= -_RS_ZERO:
            return t, u, v
    return None


def _LpSpRp(
    x: float, y: float, phi: float
) -> Optional[Tuple[float, float, float]]:
    """Formula 8.2."""
    u1, t1 = polar(x + math.sin(phi), y - 1.0 - math.cos(phi))
    u1_sq = u1 * u1
    if u1_sq >= 4.0:
        u = math.sqrt(u1_sq - 4.0)
        theta = math.atan2(2.0, u)
        t = pify(t1 + theta)
        v = pify(t - phi)
        if t >= -_RS_ZERO and v >= -_RS_ZERO:
            return t, u, v
    return None


def _CSC(x: float, y: float, phi: float) -> List[ReedsSheppPath]:
    paths: List[ReedsSheppPath] = []

    res = _LpSpLp(x, y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[14], t, u, v))

    res = _LpSpLp(-x, y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[14], -t, -u, -v))

    res = _LpSpLp(x, -y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[15], t, u, v))

    res = _LpSpLp(-x, -y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[15], -t, -u, -v))

    res = _LpSpRp(x, y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[12], t, u, v))

    res = _LpSpRp(-x, y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[12], -t, -u, -v))

    res = _LpSpRp(x, -y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[13], t, u, v))

    res = _LpSpRp(-x, -y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[13], -t, -u, -v))

    return paths


def _LpRmL(
    x: float, y: float, phi: float
) -> Optional[Tuple[float, float, float]]:
    """Formula 8.3 / 8.4 (typo in paper)."""
    xi = x - math.sin(phi)
    eta = y - 1.0 + math.cos(phi)
    u1, theta = polar(xi, eta)
    if u1 <= 4.0:
        u = -2.0 * math.asin(0.25 * u1)
        t = pify(theta + 0.5 * u + PI)
        v = pify(phi - t + u)
        if t >= -_RS_ZERO and u <= _RS_ZERO:
            return t, u, v
    return None


def _CCC(x: float, y: float, phi: float) -> List[ReedsSheppPath]:
    paths: List[ReedsSheppPath] = []

    res = _LpRmL(x, y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[0], t, u, v))

    res = _LpRmL(-x, y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[0], -t, -u, -v))

    res = _LpRmL(x, -y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[1], t, u, v))

    res = _LpRmL(-x, -y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[1], -t, -u, -v))

    # backwards
    xb = x * math.cos(phi) + y * math.sin(phi)
    yb = x * math.sin(phi) - y * math.cos(phi)

    res = _LpRmL(xb, yb, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[0], v, u, t))

    res = _LpRmL(-xb, yb, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[0], -v, -u, -t))

    res = _LpRmL(xb, -yb, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[1], v, u, t))

    res = _LpRmL(-xb, -yb, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[1], -v, -u, -t))

    return paths


def _LpRupLumRm(
    x: float, y: float, phi: float
) -> Optional[Tuple[float, float, float]]:
    """Formula 8.7."""
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho = 0.25 * (2.0 + math.sqrt(xi * xi + eta * eta))
    if rho <= 1.0:
        u = math.acos(rho)
        t, v = _tau_omega(u, -u, xi, eta, phi)
        if t >= -_RS_ZERO and v <= _RS_ZERO:
            return t, u, v
    return None


def _LpRumLumRp(
    x: float, y: float, phi: float
) -> Optional[Tuple[float, float, float]]:
    """Formula 8.8."""
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho = (20.0 - xi * xi - eta * eta) / 16.0
    if 0 <= rho <= 1:
        u = -math.acos(rho)
        if u >= -0.5 * PI:
            t, v = _tau_omega(u, u, xi, eta, phi)
            if t >= -_RS_ZERO and v >= -_RS_ZERO:
                return t, u, v
    return None


def _CCCC(x: float, y: float, phi: float) -> List[ReedsSheppPath]:
    paths: List[ReedsSheppPath] = []

    res = _LpRupLumRm(x, y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[2], t, u, -u, v))

    res = _LpRupLumRm(-x, y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[2], -t, -u, u, -v))

    res = _LpRupLumRm(x, -y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[3], t, u, -u, v))

    res = _LpRupLumRm(-x, -y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[3], -t, -u, u, -v))

    res = _LpRumLumRp(x, y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[2], t, u, u, v))

    res = _LpRumLumRp(-x, y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[2], -t, -u, -u, -v))

    res = _LpRumLumRp(x, -y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[3], t, u, u, v))

    res = _LpRumLumRp(-x, -y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[3], -t, -u, -u, -v))

    return paths


def _LpRmSmLm(
    x: float, y: float, phi: float
) -> Optional[Tuple[float, float, float]]:
    """Formula 8.9."""
    xi = x - math.sin(phi)
    eta = y - 1.0 + math.cos(phi)
    rho, theta = polar(xi, eta)
    if rho >= 2.0:
        r = math.sqrt(rho * rho - 4.0)
        u = 2.0 - r
        t = pify(theta + math.atan2(r, -2.0))
        v = pify(phi - 0.5 * PI - t)
        if t >= -_RS_ZERO and u <= _RS_ZERO and v <= _RS_ZERO:
            return t, u, v
    return None


def _LpRmSmRm(
    x: float, y: float, phi: float
) -> Optional[Tuple[float, float, float]]:
    """Formula 8.10."""
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho, theta = polar(-eta, xi)
    if rho >= 2.0:
        t = theta
        u = 2.0 - rho
        v = pify(t + 0.5 * PI - phi)
        if t >= -_RS_ZERO and u <= _RS_ZERO and v <= _RS_ZERO:
            return t, u, v
    return None


def _CCSC(x: float, y: float, phi: float) -> List[ReedsSheppPath]:
    paths: List[ReedsSheppPath] = []
    HALF_PI = PI / 2.0

    res = _LpRmSmLm(x, y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[4], t, -HALF_PI, u, v))

    res = _LpRmSmLm(-x, y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[4], -t, HALF_PI, -u, -v))

    res = _LpRmSmLm(x, -y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[5], t, -HALF_PI, u, v))

    res = _LpRmSmLm(-x, -y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[5], -t, HALF_PI, -u, -v))

    res = _LpRmSmRm(x, y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[8], t, -HALF_PI, u, v))

    res = _LpRmSmRm(-x, y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[8], -t, HALF_PI, -u, -v))

    res = _LpRmSmRm(x, -y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[9], t, -HALF_PI, u, v))

    res = _LpRmSmRm(-x, -y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[9], -t, HALF_PI, -u, -v))

    # backwards
    xb = x * math.cos(phi) + y * math.sin(phi)
    yb = x * math.sin(phi) - y * math.cos(phi)

    res = _LpRmSmLm(xb, yb, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[6], v, u, -HALF_PI, t))

    res = _LpRmSmLm(-xb, yb, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[6], -v, -u, HALF_PI, -t))

    res = _LpRmSmLm(xb, -yb, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[7], v, u, -HALF_PI, t))

    res = _LpRmSmLm(-xb, -yb, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[7], -v, -u, HALF_PI, -t))

    res = _LpRmSmRm(xb, yb, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[10], v, u, -HALF_PI, t))

    res = _LpRmSmRm(-xb, yb, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[10], -v, -u, HALF_PI, -t))

    res = _LpRmSmRm(xb, -yb, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[11], v, u, -HALF_PI, t))

    res = _LpRmSmRm(-xb, -yb, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[11], -v, -u, HALF_PI, -t))

    return paths


def _LpRmSLmRp(
    x: float, y: float, phi: float
) -> Optional[Tuple[float, float, float]]:
    """Formula 8.11 (typo in paper)."""
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho, theta = polar(xi, eta)
    if rho >= 2.0:
        u = 4.0 - math.sqrt(rho * rho - 4.0)
        if u <= _RS_ZERO:
            t = pify(math.atan2((4 - u) * xi - 2 * eta, -2 * xi + (u - 4) * eta))
            v = pify(t - phi)
            if t >= -_RS_ZERO and v >= -_RS_ZERO:
                return t, u, v
    return None


def _CCSCC(x: float, y: float, phi: float) -> List[ReedsSheppPath]:
    paths: List[ReedsSheppPath] = []
    HALF_PI = PI / 2.0

    res = _LpRmSLmRp(x, y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[16], t, -HALF_PI, u, -HALF_PI, v))

    res = _LpRmSLmRp(-x, y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[16], -t, HALF_PI, -u, HALF_PI, -v))

    res = _LpRmSLmRp(x, -y, -phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[17], t, -HALF_PI, u, -HALF_PI, v))

    res = _LpRmSLmRp(-x, -y, phi)
    if res:
        t, u, v = res
        paths.append(ReedsSheppPath(_PATH_TYPES[17], -t, HALF_PI, -u, HALF_PI, -v))

    return paths


def _get_all_rs_paths_normalised(
    x: float, y: float, phi: float
) -> List[ReedsSheppPath]:
    """Return all Reeds-Shepp paths in kappa=1 normalised coordinates."""
    paths: List[ReedsSheppPath] = []
    paths.extend(_CSC(x, y, phi))
    paths.extend(_CCC(x, y, phi))
    paths.extend(_CCCC(x, y, phi))
    paths.extend(_CCSC(x, y, phi))
    paths.extend(_CCSCC(x, y, phi))
    return paths


def _reeds_shepp_normalised(x: float, y: float, phi: float) -> ReedsSheppPath:
    """Return the shortest Reeds-Shepp path in kappa=1 coordinates."""
    best = ReedsSheppPath()
    min_len = best.length()
    for p in _get_all_rs_paths_normalised(x, y, phi):
        l = p.length()
        if l < min_len:
            min_len = l
            best = p
    return best


# ---------------------------------------------------------------------------
# Public class
# ---------------------------------------------------------------------------


class ReedsSheppStateSpace(BaseStateSpace):
    """SE(2) state space with distance measured by the length of Reeds-Shepp curves.

    Args:
        kappa: Maximum curvature (1/m).
        discretization: Step size for path integration (m).
    """

    def __init__(self, kappa: float, discretization: float = 0.1) -> None:
        super().__init__(discretization)
        assert kappa > 0.0
        self._kappa = kappa
        self._kappa_inv = 1.0 / kappa

    # ------------------------------------------------------------------
    # Core RS computation
    # ------------------------------------------------------------------

    def reeds_shepp(self, state1: State, state2: State) -> ReedsSheppPath:
        """Return the shortest RS path from *state1* to *state2* with kappa=1."""
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        dth = state2.theta - state1.theta
        c = math.cos(state1.theta)
        s = math.sin(state1.theta)
        x = c * dx + s * dy
        y = -s * dx + c * dy
        return _reeds_shepp_normalised(x * self._kappa, y * self._kappa, dth)

    def get_all_rs_paths(self, state1: State, state2: State) -> List[ReedsSheppPath]:
        """Return all RS paths between *state1* and *state2* with kappa=1."""
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        dth = state2.theta - state1.theta
        c = math.cos(state1.theta)
        s = math.sin(state1.theta)
        x = c * dx + s * dy
        y = -s * dx + c * dy
        return _get_all_rs_paths_normalised(x * self._kappa, y * self._kappa, dth)

    def get_distance(self, state1: State, state2: State) -> float:
        """Shortest RS path length from *state1* to *state2*."""
        return self._kappa_inv * self.reeds_shepp(state1, state2).length()

    # ------------------------------------------------------------------
    # Control extraction
    # ------------------------------------------------------------------

    def extract_controls_from_path(self, path: ReedsSheppPath) -> List[Control]:
        """Convert a :class:`ReedsSheppPath` to a list of :class:`Control` objects."""
        controls: List[Control] = []
        for i in range(5):
            seg = path.type[i]
            if seg == _SegType.NOP:
                break
            length = path.length_[i]
            if seg == _SegType.LEFT:
                ctrl = Control(self._kappa_inv * length, self._kappa, 0.0)
            elif seg == _SegType.RIGHT:
                ctrl = Control(self._kappa_inv * length, -self._kappa, 0.0)
            else:  # STRAIGHT
                ctrl = Control(self._kappa_inv * length, 0.0, 0.0)
            if abs(ctrl.delta_s) > _RS_EPS:
                controls.append(ctrl)
        return controls

    def get_controls(self, state1: State, state2: State) -> List[Control]:
        """Return controls for the shortest RS path."""
        return self.extract_controls_from_path(self.reeds_shepp(state1, state2))

    def get_all_controls(
        self, state1: State, state2: State
    ) -> List[List[Control]]:
        """Return all valid RS control sequences between *state1* and *state2*."""
        # Trivial case: same state
        if (
            abs(state1.x - state2.x) < 1e-6
            and abs(state1.y - state2.y) < 1e-6
            and abs(state1.theta - state2.theta) < 1e-6
        ):
            return [[Control(0.0, self._kappa, 0.0)]]

        all_controls: List[List[Control]] = []
        for rs_path in self.get_all_rs_paths(state1, state2):
            ctrl = self.extract_controls_from_path(rs_path)
            if ctrl:
                all_controls.append(ctrl)
        return all_controls
