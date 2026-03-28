# Copyright (c) 2017 - for information on the respective copyright
# owner see the NOTICE file and/or the repository
#
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
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing
# permissions and limitations under the License.
#
# This source code is derived from the Open Motion Planning Library
# (OMPL) V 1.3.1 (https://github.com/ompl/ompl).
# Copyright (c) 2010, Rice University, licensed under the BSD license,
# cf. 3rd-party-licenses.txt file in the root directory of this source
# tree.

import math
import sys

from steering_functions.utilities import pify, polar, PI
from steering_functions.base_state_space import BaseStateSpace
from steering_functions.state import State, Control

# Segment types
RS_NOP = 0
RS_LEFT = 1
RS_STRAIGHT = 2
RS_RIGHT = 3

# 18 Reeds-Shepp path type definitions (each is a tuple of 5 segment types)
reeds_shepp_path_type = [
    (RS_LEFT, RS_RIGHT, RS_LEFT, RS_NOP, RS_NOP),        # 0
    (RS_RIGHT, RS_LEFT, RS_RIGHT, RS_NOP, RS_NOP),       # 1
    (RS_LEFT, RS_RIGHT, RS_LEFT, RS_RIGHT, RS_NOP),      # 2
    (RS_RIGHT, RS_LEFT, RS_RIGHT, RS_LEFT, RS_NOP),      # 3
    (RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_NOP),   # 4
    (RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_NOP),  # 5
    (RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_LEFT, RS_NOP),   # 6
    (RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_RIGHT, RS_NOP),  # 7
    (RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP),  # 8
    (RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_NOP),   # 9
    (RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_LEFT, RS_NOP),  # 10
    (RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_RIGHT, RS_NOP),   # 11
    (RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_NOP, RS_NOP),    # 12
    (RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_NOP, RS_NOP),    # 13
    (RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_NOP, RS_NOP),     # 14
    (RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP, RS_NOP),   # 15
    (RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_RIGHT), # 16
    (RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_LEFT), # 17
]

# Helper constants
_RS_EPS = 1e-6
_RS_ZERO = 10.0 * sys.float_info.epsilon


class ReedsSheppPath:
    """Complete description of a Reeds-Shepp path."""

    def __init__(self, type_=None, t=float('inf'), u=0.0, v=0.0, w=0.0, x=0.0):
        if type_ is None:
            type_ = reeds_shepp_path_type[0]
        self.type_ = list(type_)
        self.length_ = [t, u, v, w, x]
        self.total_length_ = abs(t) + abs(u) + abs(v) + abs(w) + abs(x)

    def length(self):
        return self.total_length_

    def get_length(self, index):
        return self.length_[index]


# ---------------------------------------------------------------------------
# Helper functions (formulas from the Reeds & Shepp paper)
# ---------------------------------------------------------------------------

def _tauOmega(u, v, xi, eta, phi):
    """Compute tau and omega helper values."""
    delta = pify(u - v)
    A = math.sin(u) - math.sin(delta)
    B = math.cos(u) - math.cos(delta) - 1.0
    t1 = math.atan2(eta * A - xi * B, xi * A + eta * B)
    t2 = 2.0 * (math.cos(delta) - math.cos(v) - math.cos(u)) + 3.0
    tau = pify(t1 + PI) if t2 < 0 else pify(t1)
    omega = pify(tau - u + v - phi)
    return tau, omega


def _LpSpLp(x, y, phi):
    """Formula 8.1 in Reeds-Shepp paper."""
    u, t = polar(x - math.sin(phi), y - 1.0 + math.cos(phi))
    if t >= -_RS_ZERO:
        v = pify(phi - t)
        if v >= -_RS_ZERO:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0


def _LpSpRp(x, y, phi):
    """Formula 8.2."""
    u1, t1 = polar(x + math.sin(phi), y - 1.0 - math.cos(phi))
    u1 = u1 * u1
    if u1 >= 4.0:
        u = math.sqrt(u1 - 4.0)
        theta = math.atan2(2.0, u)
        t = pify(t1 + theta)
        v = pify(t - phi)
        if t >= -_RS_ZERO and v >= -_RS_ZERO:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0


def _LpRmL(x, y, phi):
    """Formula 8.3 / 8.4 (typo in paper)."""
    xi = x - math.sin(phi)
    eta = y - 1.0 + math.cos(phi)
    u1, theta = polar(xi, eta)
    if u1 <= 4.0:
        u = -2.0 * math.asin(0.25 * u1)
        t = pify(theta + 0.5 * u + PI)
        v = pify(phi - t + u)
        if t >= -_RS_ZERO and u <= _RS_ZERO:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0


def _LpRupLumRm(x, y, phi):
    """Formula 8.7."""
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho = 0.25 * (2.0 + math.sqrt(xi * xi + eta * eta))
    if rho <= 1.0:
        u = math.acos(rho)
        t, v = _tauOmega(u, -u, xi, eta, phi)
        if t >= -_RS_ZERO and v <= _RS_ZERO:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0


def _LpRumLumRp(x, y, phi):
    """Formula 8.8."""
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho = (20.0 - xi * xi - eta * eta) / 16.0
    if 0.0 <= rho <= 1.0:
        u = -math.acos(rho)
        if u >= -0.5 * PI:
            t, v = _tauOmega(u, u, xi, eta, phi)
            if t >= -_RS_ZERO and v >= -_RS_ZERO:
                return True, t, u, v
    return False, 0.0, 0.0, 0.0


def _LpRmSmLm(x, y, phi):
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
            return True, t, u, v
    return False, 0.0, 0.0, 0.0


def _LpRmSmRm(x, y, phi):
    """Formula 8.10."""
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho, theta = polar(-eta, xi)
    if rho >= 2.0:
        t = theta
        u = 2.0 - rho
        v = pify(t + 0.5 * PI - phi)
        if t >= -_RS_ZERO and u <= _RS_ZERO and v <= _RS_ZERO:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0


def _LpRmSLmRp(x, y, phi):
    """Formula 8.11 (typo in paper)."""
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho, theta = polar(xi, eta)
    if rho >= 2.0:
        u = 4.0 - math.sqrt(rho * rho - 4.0)
        if u <= _RS_ZERO:
            t = pify(math.atan2((4.0 - u) * xi - 2.0 * eta,
                                -2.0 * xi + (u - 4.0) * eta))
            v = pify(t - phi)
            if t >= -_RS_ZERO and v >= -_RS_ZERO:
                return True, t, u, v
    return False, 0.0, 0.0, 0.0


# ---------------------------------------------------------------------------
# Path collectors
# ---------------------------------------------------------------------------

def _CSC(x, y, phi, paths):
    """Collect all CSC paths."""
    ok, t, u, v = _LpSpLp(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[14], t, u, v))

    ok, t, u, v = _LpSpLp(-x, y, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[14], -t, -u, -v))

    ok, t, u, v = _LpSpLp(x, -y, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[15], t, u, v))

    ok, t, u, v = _LpSpLp(-x, -y, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[15], -t, -u, -v))

    ok, t, u, v = _LpSpRp(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[12], t, u, v))

    ok, t, u, v = _LpSpRp(-x, y, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[12], -t, -u, -v))

    ok, t, u, v = _LpSpRp(x, -y, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[13], t, u, v))

    ok, t, u, v = _LpSpRp(-x, -y, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[13], -t, -u, -v))


def _CCC(x, y, phi, paths):
    """Collect all CCC paths."""
    ok, t, u, v = _LpRmL(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[0], t, u, v))

    ok, t, u, v = _LpRmL(-x, y, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[0], -t, -u, -v))

    ok, t, u, v = _LpRmL(x, -y, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[1], t, u, v))

    ok, t, u, v = _LpRmL(-x, -y, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[1], -t, -u, -v))

    # backwards
    xb = x * math.cos(phi) + y * math.sin(phi)
    yb = x * math.sin(phi) - y * math.cos(phi)

    ok, t, u, v = _LpRmL(xb, yb, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[0], v, u, t))

    ok, t, u, v = _LpRmL(-xb, yb, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[0], -v, -u, -t))

    ok, t, u, v = _LpRmL(xb, -yb, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[1], v, u, t))

    ok, t, u, v = _LpRmL(-xb, -yb, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[1], -v, -u, -t))


def _CCCC(x, y, phi, paths):
    """Collect all CCCC paths."""
    ok, t, u, v = _LpRupLumRm(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[2], t, u, -u, v))

    ok, t, u, v = _LpRupLumRm(-x, y, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[2], -t, -u, u, -v))

    ok, t, u, v = _LpRupLumRm(x, -y, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[3], t, u, -u, v))

    ok, t, u, v = _LpRupLumRm(-x, -y, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[3], -t, -u, u, -v))

    ok, t, u, v = _LpRumLumRp(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[2], t, u, u, v))

    ok, t, u, v = _LpRumLumRp(-x, y, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[2], -t, -u, -u, -v))

    ok, t, u, v = _LpRumLumRp(x, -y, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[3], t, u, u, v))

    ok, t, u, v = _LpRumLumRp(-x, -y, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[3], -t, -u, -u, -v))


def _CCSC(x, y, phi, paths):
    """Collect all CCSC paths."""
    ok, t, u, v = _LpRmSmLm(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[4], t, -0.5 * PI, u, v))

    ok, t, u, v = _LpRmSmLm(-x, y, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[4], -t, 0.5 * PI, -u, -v))

    ok, t, u, v = _LpRmSmLm(x, -y, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[5], t, -0.5 * PI, u, v))

    ok, t, u, v = _LpRmSmLm(-x, -y, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[5], -t, 0.5 * PI, -u, -v))

    ok, t, u, v = _LpRmSmRm(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[8], t, -0.5 * PI, u, v))

    ok, t, u, v = _LpRmSmRm(-x, y, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[8], -t, 0.5 * PI, -u, -v))

    ok, t, u, v = _LpRmSmRm(x, -y, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[9], t, -0.5 * PI, u, v))

    ok, t, u, v = _LpRmSmRm(-x, -y, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[9], -t, 0.5 * PI, -u, -v))

    # backwards
    xb = x * math.cos(phi) + y * math.sin(phi)
    yb = x * math.sin(phi) - y * math.cos(phi)

    ok, t, u, v = _LpRmSmLm(xb, yb, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[6], v, u, -0.5 * PI, t))

    ok, t, u, v = _LpRmSmLm(-xb, yb, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[6], -v, -u, 0.5 * PI, -t))

    ok, t, u, v = _LpRmSmLm(xb, -yb, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[7], v, u, -0.5 * PI, t))

    ok, t, u, v = _LpRmSmLm(-xb, -yb, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[7], -v, -u, 0.5 * PI, -t))

    ok, t, u, v = _LpRmSmRm(xb, yb, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[10], v, u, -0.5 * PI, t))

    ok, t, u, v = _LpRmSmRm(-xb, yb, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[10], -v, -u, 0.5 * PI, -t))

    ok, t, u, v = _LpRmSmRm(xb, -yb, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[11], v, u, -0.5 * PI, t))

    ok, t, u, v = _LpRmSmRm(-xb, -yb, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[11], -v, -u, 0.5 * PI, -t))


def _CCSCC(x, y, phi, paths):
    """Collect all CCSCC paths."""
    ok, t, u, v = _LpRmSLmRp(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[16],
                                    t, -0.5 * PI, u, -0.5 * PI, v))

    ok, t, u, v = _LpRmSLmRp(-x, y, -phi)  # timeflip
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[16],
                                    -t, 0.5 * PI, -u, 0.5 * PI, -v))

    ok, t, u, v = _LpRmSLmRp(x, -y, -phi)  # reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[17],
                                    t, -0.5 * PI, u, -0.5 * PI, v))

    ok, t, u, v = _LpRmSLmRp(-x, -y, phi)  # timeflip + reflect
    if ok:
        paths.append(ReedsSheppPath(reeds_shepp_path_type[17],
                                    -t, 0.5 * PI, -u, 0.5 * PI, -v))


# ---------------------------------------------------------------------------
# Solvers (module-level, operating on normalized coordinates with kappa=1)
# ---------------------------------------------------------------------------

def _get_all_rs_paths(x, y, phi):
    """Return all valid Reeds-Shepp paths for normalized coordinates."""
    paths = []
    _CSC(x, y, phi, paths)
    _CCC(x, y, phi, paths)
    _CCCC(x, y, phi, paths)
    _CCSC(x, y, phi, paths)
    _CCSCC(x, y, phi, paths)
    return paths


def _reeds_shepp(x, y, phi):
    """Find the shortest Reeds-Shepp path for normalized coordinates."""
    path = ReedsSheppPath()
    all_paths = _get_all_rs_paths(x, y, phi)
    Lmin = path.length()
    for p in all_paths:
        L = p.length()
        if L < Lmin:
            path = p
            Lmin = L
    return path


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class ReedsSheppStateSpace(BaseStateSpace):
    """An SE(2) state space where distance is measured by the length of
    Reeds-Shepp curves.

    The notation and solutions are taken from:
    J.A. Reeds and L.A. Shepp, "Optimal paths for a car that goes both
    forwards and backwards," Pacific Journal of Mathematics,
    145(2):367-393, 1990.
    """

    def __init__(self, kappa, discretization=0.1):
        assert kappa > 0.0 and discretization > 0.0
        super().__init__(discretization)
        self.kappa_ = kappa
        self.kappa_inv_ = 1.0 / kappa

    def reeds_shepp(self, state1, state2):
        """Return type and length of segments of path from state1 to state2
        with curvature = 1.0."""
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        dth = state2.theta - state1.theta
        c = math.cos(state1.theta)
        s = math.sin(state1.theta)
        x = c * dx + s * dy
        y = -s * dx + c * dy
        return _reeds_shepp(x * self.kappa_, y * self.kappa_, dth)

    def get_all_rs_paths(self, state1, state2):
        """Return all valid Reeds-Shepp paths from state1 to state2."""
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        dth = state2.theta - state1.theta
        c = math.cos(state1.theta)
        s = math.sin(state1.theta)
        x = c * dx + s * dy
        y = -s * dx + c * dy
        return _get_all_rs_paths(x * self.kappa_, y * self.kappa_, dth)

    def get_distance(self, state1, state2):
        """Return shortest path length from state1 to state2 with
        curvature = kappa_."""
        return self.kappa_inv_ * self.reeds_shepp(state1, state2).length()

    def get_controls(self, state1, state2):
        """Return controls of the shortest path from state1 to state2 with
        curvature = kappa_."""
        path = self.reeds_shepp(state1, state2)
        return self.extract_controls_from_path(path)

    def extract_controls_from_path(self, path):
        """Extract a list of Control objects from a ReedsSheppPath."""
        controls = []
        for i in range(5):
            seg_type = path.type_[i]
            if seg_type == RS_NOP:
                return controls

            control = Control()
            if seg_type == RS_LEFT:
                control.delta_s = self.kappa_inv_ * path.length_[i]
                control.kappa = self.kappa_
                control.sigma = 0.0
            elif seg_type == RS_RIGHT:
                control.delta_s = self.kappa_inv_ * path.length_[i]
                control.kappa = -self.kappa_
                control.sigma = 0.0
            elif seg_type == RS_STRAIGHT:
                control.delta_s = self.kappa_inv_ * path.length_[i]
                control.kappa = 0.0
                control.sigma = 0.0

            if abs(control.delta_s) > _RS_EPS:
                controls.append(control)
        return controls

    def get_all_controls(self, state1, state2):
        """Return controls for all valid Reeds-Shepp paths from state1 to
        state2."""
        if (abs(state1.x - state2.x) < 1e-6 and
                abs(state1.y - state2.y) < 1e-6 and
                abs(state1.theta - state2.theta) < 1e-6):
            return [[Control(0.0, self.kappa_, 0.0)]]

        all_rs_paths = self.get_all_rs_paths(state1, state2)
        all_controls = []
        for rs_path in all_rs_paths:
            control = self.extract_controls_from_path(rs_path)
            if len(control) > 0:
                all_controls.append(control)
        return all_controls
