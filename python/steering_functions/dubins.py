from __future__ import annotations

from dataclasses import dataclass
import math
from typing import List, Sequence, Tuple

from .core import BaseStateSpace, Control, State, twopify

DUBINS_EPS = 1e-6
DUBINS_ZERO = -1e-9

LEFT = "L"
STRAIGHT = "S"
RIGHT = "R"

_DUBINS_PATH_TYPES: tuple[tuple[str, str, str], ...] = (
    (LEFT, STRAIGHT, LEFT),
    (RIGHT, STRAIGHT, RIGHT),
    (RIGHT, STRAIGHT, LEFT),
    (LEFT, STRAIGHT, RIGHT),
    (RIGHT, LEFT, RIGHT),
    (LEFT, RIGHT, LEFT),
)


@dataclass
class DubinsPath:
    type_: Sequence[str] = _DUBINS_PATH_TYPES[0]
    lengths: Tuple[float, float, float] = (0.0, math.inf, 0.0)

    def length(self) -> float:
        return sum(self.lengths)


def _dubins_lsl(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa, cb, sb = math.cos(alpha), math.sin(alpha), math.cos(beta), math.sin(beta)
    tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sa - sb))
    if tmp >= DUBINS_ZERO:
        theta = math.atan2(cb - ca, d + sa - sb)
        t = twopify(-alpha + theta)
        p = math.sqrt(max(tmp, 0.0))
        q = twopify(beta - theta)
        return DubinsPath(_DUBINS_PATH_TYPES[0], (t, p, q))
    return DubinsPath()


def _dubins_rsr(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa, cb, sb = math.cos(alpha), math.sin(alpha), math.cos(beta), math.sin(beta)
    tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sb - sa))
    if tmp >= DUBINS_ZERO:
        theta = math.atan2(ca - cb, d - sa + sb)
        t = twopify(alpha - theta)
        p = math.sqrt(max(tmp, 0.0))
        q = twopify(-beta + theta)
        return DubinsPath(_DUBINS_PATH_TYPES[1], (t, p, q))
    return DubinsPath()


def _dubins_rsl(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa, cb, sb = math.cos(alpha), math.sin(alpha), math.cos(beta), math.sin(beta)
    tmp = d * d - 2.0 + 2.0 * (ca * cb + sa * sb - d * (sa + sb))
    if tmp >= DUBINS_ZERO:
        p = math.sqrt(max(tmp, 0.0))
        theta = math.atan2(ca + cb, d - sa - sb) - math.atan2(2.0, p)
        t = twopify(alpha - theta)
        q = twopify(beta - theta)
        return DubinsPath(_DUBINS_PATH_TYPES[2], (t, p, q))
    return DubinsPath()


def _dubins_lsr(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa, cb, sb = math.cos(alpha), math.sin(alpha), math.cos(beta), math.sin(beta)
    tmp = -2.0 + d * d + 2.0 * (ca * cb + sa * sb + d * (sa + sb))
    if tmp >= DUBINS_ZERO:
        p = math.sqrt(max(tmp, 0.0))
        theta = math.atan2(-ca - cb, d + sa + sb) - math.atan2(-2.0, p)
        t = twopify(-alpha + theta)
        q = twopify(-beta + theta)
        return DubinsPath(_DUBINS_PATH_TYPES[3], (t, p, q))
    return DubinsPath()


def _dubins_rlr(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa, cb, sb = math.cos(alpha), math.sin(alpha), math.cos(beta), math.sin(beta)
    tmp = 0.125 * (6.0 - d * d + 2.0 * (ca * cb + sa * sb + d * (sa - sb)))
    if abs(tmp) < 1.0:
        p = 2.0 * math.pi - math.acos(tmp)
        theta = math.atan2(ca - cb, d - sa + sb)
        t = twopify(alpha - theta + 0.5 * p)
        q = twopify(alpha - beta - t + p)
        return DubinsPath(_DUBINS_PATH_TYPES[4], (t, p, q))
    return DubinsPath()


def _dubins_lrl(d: float, alpha: float, beta: float) -> DubinsPath:
    ca, sa, cb, sb = math.cos(alpha), math.sin(alpha), math.cos(beta), math.sin(beta)
    tmp = 0.125 * (6.0 - d * d + 2.0 * (ca * cb + sa * sb - d * (sa - sb)))
    if abs(tmp) < 1.0:
        p = 2.0 * math.pi - math.acos(tmp)
        theta = math.atan2(-ca + cb, d + sa - sb)
        t = twopify(-alpha + theta + 0.5 * p)
        q = twopify(beta - alpha - t + p)
        return DubinsPath(_DUBINS_PATH_TYPES[5], (t, p, q))
    return DubinsPath()


def _select_dubins_path(d: float, alpha: float, beta: float) -> DubinsPath:
    if d < DUBINS_EPS and abs(alpha - beta) < DUBINS_EPS:
        return DubinsPath(_DUBINS_PATH_TYPES[0], (0.0, d, 0.0))

    candidates = [
        _dubins_lsl(d, alpha, beta),
        _dubins_rsr(d, alpha, beta),
        _dubins_rsl(d, alpha, beta),
        _dubins_lsr(d, alpha, beta),
        _dubins_rlr(d, alpha, beta),
        _dubins_lrl(d, alpha, beta),
    ]
    return min(candidates, key=lambda path: path.length())


class DubinsStateSpace(BaseStateSpace):
    def __init__(self, kappa: float, discretization: float = 0.1, forwards: bool = True) -> None:
        super().__init__(discretization)
        if kappa <= 0.0:
            raise ValueError("kappa must be positive")
        self._kappa = kappa
        self._kappa_inv = 1.0 / kappa
        self._forwards = forwards

    def dubins(self, state1: State, state2: State) -> DubinsPath:
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        theta = math.atan2(dy, dx)
        d = math.hypot(dx, dy) * self._kappa
        alpha = twopify(state1.theta - theta)
        beta = twopify(state2.theta - theta)
        return _select_dubins_path(d, alpha, beta)

    def get_distance(self, state1: State, state2: State) -> float:
        if self._forwards:
            return self._kappa_inv * self.dubins(state1, state2).length()
        return self._kappa_inv * self.dubins(state2, state1).length()

    def _get_controls_core(
        self,
        state1: State,
        state2: State,
        reverse_direction: bool = False,
    ) -> List[Control]:
        path = self.dubins(state2, state1) if reverse_direction else self.dubins(state1, state2)
        controls: List[Control] = []

        for segment_type, length in zip(path.type_, path.lengths):
            if segment_type == LEFT:
                control = Control(self._kappa_inv * length, self._kappa, 0.0)
            elif segment_type == RIGHT:
                control = Control(self._kappa_inv * length, -self._kappa, 0.0)
            else:
                control = Control(self._kappa_inv * length, 0.0, 0.0)

            if abs(control.delta_s) > DUBINS_EPS:
                controls.append(control)

        if reverse_direction:
            controls.reverse()
            controls = [Control(-control.delta_s, control.kappa, control.sigma) for control in controls]

        return controls

    def get_controls(self, state1: State, state2: State) -> List[Control]:
        return self._get_controls_core(state1, state2, False)

    def get_controls_reverse(self, state1: State, state2: State) -> List[Control]:
        return self._get_controls_core(state1, state2, True)

    def get_all_controls(self, state1: State, state2: State) -> List[List[Control]]:
        return [self.get_controls(state1, state2), self.get_controls_reverse(state1, state2)]

    def get_path(self, state1: State, state2: State) -> List[State]:
        forward_length = self.get_distance(state1, state2)
        reverse_length = self._kappa_inv * self.dubins(state2, state1).length()
        controls = (
            self._get_controls_core(state1, state2, False)
            if forward_length <= reverse_length
            else self._get_controls_core(state1, state2, True)
        )
        return self.integrate(state1, controls)
