from __future__ import annotations

from dataclasses import dataclass
import math
from typing import List, Sequence

from .core import BaseStateSpace, Control, State, PI, pify, polar

RS_EPS = 1e-6
RS_ZERO = 10.0 * float.fromhex("0x1.0000000000000p-52")

NOP = "N"
LEFT = "L"
STRAIGHT = "S"
RIGHT = "R"

_REEDS_SHEPP_PATH_TYPES: tuple[tuple[str, str, str, str, str], ...] = (
    (LEFT, RIGHT, LEFT, NOP, NOP),
    (RIGHT, LEFT, RIGHT, NOP, NOP),
    (LEFT, RIGHT, LEFT, RIGHT, NOP),
    (RIGHT, LEFT, RIGHT, LEFT, NOP),
    (LEFT, RIGHT, STRAIGHT, LEFT, NOP),
    (RIGHT, LEFT, STRAIGHT, RIGHT, NOP),
    (LEFT, STRAIGHT, RIGHT, LEFT, NOP),
    (RIGHT, STRAIGHT, LEFT, RIGHT, NOP),
    (LEFT, RIGHT, STRAIGHT, RIGHT, NOP),
    (RIGHT, LEFT, STRAIGHT, LEFT, NOP),
    (RIGHT, STRAIGHT, RIGHT, LEFT, NOP),
    (LEFT, STRAIGHT, LEFT, RIGHT, NOP),
    (LEFT, STRAIGHT, RIGHT, NOP, NOP),
    (RIGHT, STRAIGHT, LEFT, NOP, NOP),
    (LEFT, STRAIGHT, LEFT, NOP, NOP),
    (RIGHT, STRAIGHT, RIGHT, NOP, NOP),
    (LEFT, RIGHT, STRAIGHT, LEFT, RIGHT),
    (RIGHT, LEFT, STRAIGHT, RIGHT, LEFT),
)


@dataclass
class ReedsSheppPath:
    type_: Sequence[str] = _REEDS_SHEPP_PATH_TYPES[0]
    lengths: tuple[float, float, float, float, float] = (math.inf, 0.0, 0.0, 0.0, 0.0)

    def length(self) -> float:
        return sum(abs(length) for length in self.lengths)


def _tau_omega(u: float, v: float, xi: float, eta: float, phi: float) -> tuple[float, float]:
    delta = pify(u - v)
    a = math.sin(u) - math.sin(delta)
    b = math.cos(u) - math.cos(delta) - 1.0
    t1 = math.atan2(eta * a - xi * b, xi * a + eta * b)
    t2 = 2.0 * (math.cos(delta) - math.cos(v) - math.cos(u)) + 3.0
    tau = pify(t1 + PI) if t2 < 0.0 else pify(t1)
    omega = pify(tau - u + v - phi)
    return tau, omega


def _lp_sp_lp(x: float, y: float, phi: float) -> tuple[bool, float, float, float]:
    u, t = polar(x - math.sin(phi), y - 1.0 + math.cos(phi))
    if t >= -RS_ZERO:
        v = pify(phi - t)
        if v >= -RS_ZERO:
            return True, t, u, v
    return False, 0.0, 0.0, 0.0


def _lp_sp_rp(x: float, y: float, phi: float) -> tuple[bool, float, float, float]:
    u1, t1 = polar(x + math.sin(phi), y - 1.0 - math.cos(phi))
    u1 = u1 * u1
    if u1 >= 4.0:
        u = math.sqrt(u1 - 4.0)
        theta = math.atan2(2.0, u)
        t = pify(t1 + theta)
        v = pify(t - phi)
        return t >= -RS_ZERO and v >= -RS_ZERO, t, u, v
    return False, 0.0, 0.0, 0.0


def _csc(x: float, y: float, phi: float, paths: List[ReedsSheppPath]) -> None:
    ok, t, u, v = _lp_sp_lp(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[14], (t, u, v, 0.0, 0.0)))
    ok, t, u, v = _lp_sp_lp(-x, y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[14], (-t, -u, -v, 0.0, 0.0)))
    ok, t, u, v = _lp_sp_lp(x, -y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[15], (t, u, v, 0.0, 0.0)))
    ok, t, u, v = _lp_sp_lp(-x, -y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[15], (-t, -u, -v, 0.0, 0.0)))

    ok, t, u, v = _lp_sp_rp(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[12], (t, u, v, 0.0, 0.0)))
    ok, t, u, v = _lp_sp_rp(-x, y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[12], (-t, -u, -v, 0.0, 0.0)))
    ok, t, u, v = _lp_sp_rp(x, -y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[13], (t, u, v, 0.0, 0.0)))
    ok, t, u, v = _lp_sp_rp(-x, -y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[13], (-t, -u, -v, 0.0, 0.0)))


def _lp_rm_l(x: float, y: float, phi: float) -> tuple[bool, float, float, float]:
    xi = x - math.sin(phi)
    eta = y - 1.0 + math.cos(phi)
    u1, theta = polar(xi, eta)
    if u1 <= 4.0:
        u = -2.0 * math.asin(0.25 * u1)
        t = pify(theta + 0.5 * u + PI)
        v = pify(phi - t + u)
        return t >= -RS_ZERO and u <= RS_ZERO, t, u, v
    return False, 0.0, 0.0, 0.0


def _ccc(x: float, y: float, phi: float, paths: List[ReedsSheppPath]) -> None:
    ok, t, u, v = _lp_rm_l(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[0], (t, u, v, 0.0, 0.0)))
    ok, t, u, v = _lp_rm_l(-x, y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[0], (-t, -u, -v, 0.0, 0.0)))
    ok, t, u, v = _lp_rm_l(x, -y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[1], (t, u, v, 0.0, 0.0)))
    ok, t, u, v = _lp_rm_l(-x, -y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[1], (-t, -u, -v, 0.0, 0.0)))

    xb = x * math.cos(phi) + y * math.sin(phi)
    yb = x * math.sin(phi) - y * math.cos(phi)
    ok, t, u, v = _lp_rm_l(xb, yb, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[0], (v, u, t, 0.0, 0.0)))
    ok, t, u, v = _lp_rm_l(-xb, yb, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[0], (-v, -u, -t, 0.0, 0.0)))
    ok, t, u, v = _lp_rm_l(xb, -yb, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[1], (v, u, t, 0.0, 0.0)))
    ok, t, u, v = _lp_rm_l(-xb, -yb, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[1], (-v, -u, -t, 0.0, 0.0)))


def _lp_rup_lum_rm(x: float, y: float, phi: float) -> tuple[bool, float, float, float]:
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho = 0.25 * (2.0 + math.sqrt(xi * xi + eta * eta))
    if rho <= 1.0:
        u = math.acos(rho)
        t, v = _tau_omega(u, -u, xi, eta, phi)
        return t >= -RS_ZERO and v <= RS_ZERO, t, u, v
    return False, 0.0, 0.0, 0.0


def _lp_rum_lum_rp(x: float, y: float, phi: float) -> tuple[bool, float, float, float]:
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho = (20.0 - xi * xi - eta * eta) / 16.0
    if 0.0 <= rho <= 1.0:
        u = -math.acos(rho)
        if u >= -0.5 * PI:
            t, v = _tau_omega(u, u, xi, eta, phi)
            return t >= -RS_ZERO and v >= -RS_ZERO, t, u, v
    return False, 0.0, 0.0, 0.0


def _cccc(x: float, y: float, phi: float, paths: List[ReedsSheppPath]) -> None:
    ok, t, u, v = _lp_rup_lum_rm(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[2], (t, u, -u, v, 0.0)))
    ok, t, u, v = _lp_rup_lum_rm(-x, y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[2], (-t, -u, u, -v, 0.0)))
    ok, t, u, v = _lp_rup_lum_rm(x, -y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[3], (t, u, -u, v, 0.0)))
    ok, t, u, v = _lp_rup_lum_rm(-x, -y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[3], (-t, -u, u, -v, 0.0)))

    ok, t, u, v = _lp_rum_lum_rp(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[2], (t, u, u, v, 0.0)))
    ok, t, u, v = _lp_rum_lum_rp(-x, y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[2], (-t, -u, -u, -v, 0.0)))
    ok, t, u, v = _lp_rum_lum_rp(x, -y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[3], (t, u, u, v, 0.0)))
    ok, t, u, v = _lp_rum_lum_rp(-x, -y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[3], (-t, -u, -u, -v, 0.0)))


def _lp_rm_sm_lm(x: float, y: float, phi: float) -> tuple[bool, float, float, float]:
    xi = x - math.sin(phi)
    eta = y - 1.0 + math.cos(phi)
    rho, theta = polar(xi, eta)
    if rho >= 2.0:
        r = math.sqrt(rho * rho - 4.0)
        u = 2.0 - r
        t = pify(theta + math.atan2(r, -2.0))
        v = pify(phi - 0.5 * PI - t)
        return t >= -RS_ZERO and u <= RS_ZERO and v <= RS_ZERO, t, u, v
    return False, 0.0, 0.0, 0.0


def _lp_rm_sm_rm(x: float, y: float, phi: float) -> tuple[bool, float, float, float]:
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho, theta = polar(-eta, xi)
    if rho >= 2.0:
        t = theta
        u = 2.0 - rho
        v = pify(t + 0.5 * PI - phi)
        return t >= -RS_ZERO and u <= RS_ZERO and v <= RS_ZERO, t, u, v
    return False, 0.0, 0.0, 0.0


def _ccsc(x: float, y: float, phi: float, paths: List[ReedsSheppPath]) -> None:
    ok, t, u, v = _lp_rm_sm_lm(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[4], (t, -0.5 * PI, u, v, 0.0)))
    ok, t, u, v = _lp_rm_sm_lm(-x, y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[4], (-t, 0.5 * PI, -u, -v, 0.0)))
    ok, t, u, v = _lp_rm_sm_lm(x, -y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[5], (t, -0.5 * PI, u, v, 0.0)))
    ok, t, u, v = _lp_rm_sm_lm(-x, -y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[5], (-t, 0.5 * PI, -u, -v, 0.0)))

    ok, t, u, v = _lp_rm_sm_rm(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[8], (t, -0.5 * PI, u, v, 0.0)))
    ok, t, u, v = _lp_rm_sm_rm(-x, y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[8], (-t, 0.5 * PI, -u, -v, 0.0)))
    ok, t, u, v = _lp_rm_sm_rm(x, -y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[9], (t, -0.5 * PI, u, v, 0.0)))
    ok, t, u, v = _lp_rm_sm_rm(-x, -y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[9], (-t, 0.5 * PI, -u, -v, 0.0)))

    xb = x * math.cos(phi) + y * math.sin(phi)
    yb = x * math.sin(phi) - y * math.cos(phi)
    ok, t, u, v = _lp_rm_sm_lm(xb, yb, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[6], (v, u, -0.5 * PI, t, 0.0)))
    ok, t, u, v = _lp_rm_sm_lm(-xb, yb, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[6], (-v, -u, 0.5 * PI, -t, 0.0)))
    ok, t, u, v = _lp_rm_sm_lm(xb, -yb, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[7], (v, u, -0.5 * PI, t, 0.0)))
    ok, t, u, v = _lp_rm_sm_lm(-xb, -yb, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[7], (-v, -u, 0.5 * PI, -t, 0.0)))
    ok, t, u, v = _lp_rm_sm_rm(xb, yb, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[10], (v, u, -0.5 * PI, t, 0.0)))
    ok, t, u, v = _lp_rm_sm_rm(-xb, yb, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[10], (-v, -u, 0.5 * PI, -t, 0.0)))
    ok, t, u, v = _lp_rm_sm_rm(xb, -yb, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[11], (v, u, -0.5 * PI, t, 0.0)))
    ok, t, u, v = _lp_rm_sm_rm(-xb, -yb, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[11], (-v, -u, 0.5 * PI, -t, 0.0)))


def _lp_rm_slm_rp(x: float, y: float, phi: float) -> tuple[bool, float, float, float]:
    xi = x + math.sin(phi)
    eta = y - 1.0 - math.cos(phi)
    rho, _ = polar(xi, eta)
    if rho >= 2.0:
        u = 4.0 - math.sqrt(rho * rho - 4.0)
        if u <= RS_ZERO:
            t = pify(math.atan2((4.0 - u) * xi - 2.0 * eta, -2.0 * xi + (u - 4.0) * eta))
            v = pify(t - phi)
            return t >= -RS_ZERO and v >= -RS_ZERO, t, u, v
    return False, 0.0, 0.0, 0.0


def _ccscc(x: float, y: float, phi: float, paths: List[ReedsSheppPath]) -> None:
    ok, t, u, v = _lp_rm_slm_rp(x, y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[16], (t, -0.5 * PI, u, -0.5 * PI, v)))
    ok, t, u, v = _lp_rm_slm_rp(-x, y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[16], (-t, 0.5 * PI, -u, 0.5 * PI, -v)))
    ok, t, u, v = _lp_rm_slm_rp(x, -y, -phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[17], (t, -0.5 * PI, u, -0.5 * PI, v)))
    ok, t, u, v = _lp_rm_slm_rp(-x, -y, phi)
    if ok:
        paths.append(ReedsSheppPath(_REEDS_SHEPP_PATH_TYPES[17], (-t, 0.5 * PI, -u, 0.5 * PI, -v)))


def _get_all_rs_paths(x: float, y: float, phi: float) -> List[ReedsSheppPath]:
    paths: List[ReedsSheppPath] = []
    _csc(x, y, phi, paths)
    _ccc(x, y, phi, paths)
    _cccc(x, y, phi, paths)
    _ccsc(x, y, phi, paths)
    _ccscc(x, y, phi, paths)
    return paths


def _select_rs_path(x: float, y: float, phi: float) -> ReedsSheppPath:
    return min(_get_all_rs_paths(x, y, phi), key=lambda path: path.length(), default=ReedsSheppPath())


class ReedsSheppStateSpace(BaseStateSpace):
    def __init__(self, kappa: float, discretization: float = 0.1) -> None:
        super().__init__(discretization)
        if kappa <= 0.0:
            raise ValueError("kappa must be positive")
        self._kappa = kappa
        self._kappa_inv = 1.0 / kappa

    def reeds_shepp(self, state1: State, state2: State) -> ReedsSheppPath:
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        dth = state2.theta - state1.theta
        c = math.cos(state1.theta)
        s = math.sin(state1.theta)
        x = c * dx + s * dy
        y = -s * dx + c * dy
        return _select_rs_path(x * self._kappa, y * self._kappa, dth)

    def get_all_rs_paths(self, state1: State, state2: State) -> List[ReedsSheppPath]:
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        dth = state2.theta - state1.theta
        c = math.cos(state1.theta)
        s = math.sin(state1.theta)
        x = c * dx + s * dy
        y = -s * dx + c * dy
        return _get_all_rs_paths(x * self._kappa, y * self._kappa, dth)

    def get_distance(self, state1: State, state2: State) -> float:
        return self._kappa_inv * self.reeds_shepp(state1, state2).length()

    def extract_controls_from_path(self, path: ReedsSheppPath) -> List[Control]:
        controls: List[Control] = []
        for segment_type, length in zip(path.type_, path.lengths):
            if segment_type == NOP:
                return controls
            if segment_type == LEFT:
                control = Control(self._kappa_inv * length, self._kappa, 0.0)
            elif segment_type == RIGHT:
                control = Control(self._kappa_inv * length, -self._kappa, 0.0)
            else:
                control = Control(self._kappa_inv * length, 0.0, 0.0)
            if abs(control.delta_s) > RS_EPS:
                controls.append(control)
        return controls

    def get_controls(self, state1: State, state2: State) -> List[Control]:
        return self.extract_controls_from_path(self.reeds_shepp(state1, state2))

    def get_all_controls(self, state1: State, state2: State) -> List[List[Control]]:
        if (
            abs(state1.x - state2.x) < RS_EPS
            and abs(state1.y - state2.y) < RS_EPS
            and abs(state1.theta - state2.theta) < RS_EPS
        ):
            return [[Control(0.0, self._kappa, 0.0)]]
        controls = [self.extract_controls_from_path(path) for path in self.get_all_rs_paths(state1, state2)]
        return [control for control in controls if control]
