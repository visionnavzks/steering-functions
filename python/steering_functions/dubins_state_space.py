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
#   A.M. Shkel and V. Lumelsky, "Classification of the Dubins set,"
#   Robotics and Autonomous Systems, 34(4):179-202, 2001.

"""Dubins state space: shortest forward-only curved paths for a car-like robot."""

import math
import sys
from enum import IntEnum
from typing import List, Optional

from .base_state_space import BaseStateSpace
from .state import Control, State
from .utilities import TWO_PI, twopify

_DUBINS_EPS = 1e-6
_DUBINS_ZERO = -1e-9


class _SegType(IntEnum):
    LEFT = 0
    STRAIGHT = 1
    RIGHT = 2


# 6 Dubins path types: LSL, RSR, RSL, LSR, RLR, LRL
_PATH_TYPES = [
    (_SegType.LEFT, _SegType.STRAIGHT, _SegType.LEFT),
    (_SegType.RIGHT, _SegType.STRAIGHT, _SegType.RIGHT),
    (_SegType.RIGHT, _SegType.STRAIGHT, _SegType.LEFT),
    (_SegType.LEFT, _SegType.STRAIGHT, _SegType.RIGHT),
    (_SegType.RIGHT, _SegType.LEFT, _SegType.RIGHT),
    (_SegType.LEFT, _SegType.RIGHT, _SegType.LEFT),
]


class DubinsPath:
    """Complete description of a Dubins path (three segments)."""

    def __init__(
        self,
        path_type: tuple = _PATH_TYPES[0],
        t: float = 0.0,
        p: float = sys.float_info.max,
        q: float = 0.0,
    ) -> None:
        self.type = path_type
        self.length_ = [t, p, q]
        assert t >= 0.0
        assert p >= 0.0
        assert q >= 0.0

    def length(self) -> float:
        return self.length_[0] + self.length_[1] + self.length_[2]


# ---------------------------------------------------------------------------
# Helper functions (closed-form Dubins primitives)
# ---------------------------------------------------------------------------


def _dubinsLSL(d: float, alpha: float, beta: float) -> Optional[DubinsPath]:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sa - sb))
    if tmp >= _DUBINS_ZERO:
        theta = math.atan2(cb - ca, d + sa - sb)
        t = twopify(-alpha + theta)
        p = math.sqrt(max(tmp, 0.0))
        q = twopify(beta - theta)
        return DubinsPath(_PATH_TYPES[0], t, p, q)
    return None


def _dubinsRSR(d: float, alpha: float, beta: float) -> Optional[DubinsPath]:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sb - sa))
    if tmp >= _DUBINS_ZERO:
        theta = math.atan2(ca - cb, d - sa + sb)
        t = twopify(alpha - theta)
        p = math.sqrt(max(tmp, 0.0))
        q = twopify(-beta + theta)
        return DubinsPath(_PATH_TYPES[1], t, p, q)
    return None


def _dubinsRSL(d: float, alpha: float, beta: float) -> Optional[DubinsPath]:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = d * d - 2.0 + 2.0 * (ca * cb + sa * sb - d * (sa + sb))
    if tmp >= _DUBINS_ZERO:
        p = math.sqrt(max(tmp, 0.0))
        theta = math.atan2(ca + cb, d - sa - sb) - math.atan2(2.0, p)
        t = twopify(alpha - theta)
        q = twopify(beta - theta)
        return DubinsPath(_PATH_TYPES[2], t, p, q)
    return None


def _dubinsLSR(d: float, alpha: float, beta: float) -> Optional[DubinsPath]:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = -2.0 + d * d + 2.0 * (ca * cb + sa * sb + d * (sa + sb))
    if tmp >= _DUBINS_ZERO:
        p = math.sqrt(max(tmp, 0.0))
        theta = math.atan2(-ca - cb, d + sa + sb) - math.atan2(-2.0, p)
        t = twopify(-alpha + theta)
        q = twopify(-beta + theta)
        return DubinsPath(_PATH_TYPES[3], t, p, q)
    return None


def _dubinsRLR(d: float, alpha: float, beta: float) -> Optional[DubinsPath]:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = 0.125 * (6.0 - d * d + 2.0 * (ca * cb + sa * sb + d * (sa - sb)))
    if abs(tmp) < 1.0:
        p = TWO_PI - math.acos(tmp)
        theta = math.atan2(ca - cb, d - sa + sb)
        t = twopify(alpha - theta + 0.5 * p)
        q = twopify(alpha - beta - t + p)
        return DubinsPath(_PATH_TYPES[4], t, p, q)
    return None


def _dubinsLRL(d: float, alpha: float, beta: float) -> Optional[DubinsPath]:
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    tmp = 0.125 * (6.0 - d * d + 2.0 * (ca * cb + sa * sb - d * (sa - sb)))
    if abs(tmp) < 1.0:
        p = TWO_PI - math.acos(tmp)
        theta = math.atan2(-ca + cb, d + sa - sb)
        t = twopify(-alpha + theta + 0.5 * p)
        q = twopify(beta - alpha - t + p)
        return DubinsPath(_PATH_TYPES[5], t, p, q)
    return None


def _dubins(d: float, alpha: float, beta: float) -> DubinsPath:
    """Return the shortest Dubins path in normalised (kappa=1) coordinates."""
    if d < _DUBINS_EPS and abs(alpha - beta) < _DUBINS_EPS:
        return DubinsPath(_PATH_TYPES[0], 0.0, d, 0.0)

    candidates = [
        _dubinsLSL(d, alpha, beta),
        _dubinsRSR(d, alpha, beta),
        _dubinsRSL(d, alpha, beta),
        _dubinsLSR(d, alpha, beta),
        _dubinsRLR(d, alpha, beta),
        _dubinsLRL(d, alpha, beta),
    ]
    best = DubinsPath()
    min_len = best.length()
    for cand in candidates:
        if cand is not None:
            l = cand.length()
            if l < min_len:
                min_len = l
                best = cand
    return best


# ---------------------------------------------------------------------------
# Public class
# ---------------------------------------------------------------------------


class DubinsStateSpace(BaseStateSpace):
    """SE(2) state space with distance measured by the length of Dubins curves.

    Args:
        kappa: Maximum curvature (1/m).
        discretization: Step size for path integration (m).
        forwards: If ``True`` (default) the car drives forward.
    """

    def __init__(
        self, kappa: float, discretization: float = 0.1, forwards: bool = True
    ) -> None:
        super().__init__(discretization)
        assert kappa > 0.0
        self._kappa = kappa
        self._kappa_inv = 1.0 / kappa
        self._forwards = forwards

    # ------------------------------------------------------------------
    # Core Dubins computation
    # ------------------------------------------------------------------

    def dubins(self, state1: State, state2: State) -> DubinsPath:
        """Return Dubins path from *state1* to *state2* with kappa = 1."""
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        th = math.atan2(dy, dx)
        d = math.sqrt(dx * dx + dy * dy) * self._kappa
        alpha = twopify(state1.theta - th)
        beta = twopify(state2.theta - th)
        return _dubins(d, alpha, beta)

    def get_distance(self, state1: State, state2: State) -> float:
        """Shortest Dubins path length from *state1* to *state2*."""
        if self._forwards:
            return self._kappa_inv * self.dubins(state1, state2).length()
        else:
            return self._kappa_inv * self.dubins(state2, state1).length()

    # ------------------------------------------------------------------
    # Control extraction
    # ------------------------------------------------------------------

    def _get_controls_core(
        self, state1: State, state2: State, reverse: bool = False
    ) -> List[Control]:
        """Core helper that optionally reverses the path direction."""
        if not reverse:
            path = self.dubins(state1, state2)
        else:
            path = self.dubins(state2, state1)

        controls: List[Control] = []
        for i in range(3):
            seg = path.type[i]
            length = path.length_[i]
            if seg == _SegType.LEFT:
                ctrl = Control(self._kappa_inv * length, self._kappa, 0.0)
            elif seg == _SegType.STRAIGHT:
                ctrl = Control(self._kappa_inv * length, 0.0, 0.0)
            else:  # RIGHT
                ctrl = Control(self._kappa_inv * length, -self._kappa, 0.0)
            if abs(ctrl.delta_s) > _DUBINS_EPS:
                controls.append(ctrl)

        if reverse:
            controls.reverse()
            for ctrl in controls:
                ctrl.delta_s = -ctrl.delta_s

        return controls

    def get_controls(self, state1: State, state2: State) -> List[Control]:
        """Return controls for the forward shortest path."""
        return self._get_controls_core(state1, state2, reverse=False)

    def get_controls_reverse(self, state1: State, state2: State) -> List[Control]:
        """Return controls for the pure reverse shortest path."""
        return self._get_controls_core(state1, state2, reverse=True)

    def get_all_controls(
        self, state1: State, state2: State
    ) -> List[List[Control]]:
        return [
            self.get_controls(state1, state2),
            self.get_controls_reverse(state1, state2),
        ]

    # ------------------------------------------------------------------
    # Path generation
    # ------------------------------------------------------------------

    def get_path(  # type: ignore[override]
        self, state1: State, state2: State, controls: List[Control] | None = None
    ) -> List[State]:
        """Return the shortest path, choosing between forward and reverse."""
        forward_len = self.get_distance(state1, state2)
        reverse_len = self._kappa_inv * self.dubins(state2, state1).length()

        if forward_len <= reverse_len:
            ctrl = self._get_controls_core(state1, state2, reverse=False)
        else:
            ctrl = self._get_controls_core(state1, state2, reverse=True)

        if controls is not None:
            controls.clear()
            controls.extend(ctrl)
        return self.integrate(state1, ctrl)

    def interpolate_multi(
        self,
        state: State,
        controls: List[Control],
        ts: List[float],
    ) -> List[State]:
        """Return states at multiple normalised distances *ts* along the path."""
        state_inters: List[State] = []
        if not controls or not ts:
            return state_inters

        s_path = sum(abs(c.delta_s) for c in controls)
        s_targets = [min(max(t, 0.0), 1.0) * s_path for t in ts]

        state_curr = State(**vars(state))
        s = 0.0
        idx = 0
        for control in controls:
            if idx >= len(s_targets):
                break
            abs_delta_s = abs(control.delta_s)
            while idx < len(s_targets):
                s_inter = s_targets[idx]
                if s_inter - s > abs_delta_s:
                    break
                state_inters.append(
                    self.integrate_ODE(state_curr, control, s_inter - s)
                )
                idx += 1
            state_curr = self.integrate_ODE(state_curr, control, abs_delta_s)
            s += abs_delta_s

        return state_inters
