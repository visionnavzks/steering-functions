from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
import math
from typing import Iterable, List

EPSILON = 1e-6
PI = math.pi
TWO_PI = 2.0 * math.pi


def sgn(value: float) -> float:
    return -1.0 if value < 0.0 else 1.0


def point_distance(x1: float, y1: float, x2: float, y2: float) -> float:
    return math.hypot(x2 - x1, y2 - y1)


def polar(x: float, y: float) -> tuple[float, float]:
    return math.hypot(x, y), math.atan2(y, x)


def twopify(alpha: float) -> float:
    return alpha - TWO_PI * math.floor(alpha / TWO_PI)


def pify(alpha: float) -> float:
    value = math.fmod(alpha, TWO_PI)
    if value < -PI:
        value += TWO_PI
    elif value > PI:
        value -= TWO_PI
    return value


def end_of_circular_arc(
    x_i: float,
    y_i: float,
    theta_i: float,
    kappa: float,
    direction: float,
    length: float,
) -> tuple[float, float, float]:
    x_f = x_i + (1.0 / kappa) * (-math.sin(theta_i) + math.sin(theta_i + direction * length * kappa))
    y_f = y_i + (1.0 / kappa) * (math.cos(theta_i) - math.cos(theta_i + direction * length * kappa))
    theta_f = pify(theta_i + kappa * direction * length)
    return x_f, y_f, theta_f


def end_of_straight_line(
    x_i: float,
    y_i: float,
    theta: float,
    direction: float,
    length: float,
) -> tuple[float, float]:
    x_f = x_i + direction * length * math.cos(theta)
    y_f = y_i + direction * length * math.sin(theta)
    return x_f, y_f


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

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, State):
            return False
        return all(
            abs(lhs - rhs) < EPSILON
            for lhs, rhs in (
                (self.x, other.x),
                (self.y, other.y),
                (self.theta, other.theta),
                (self.kappa, other.kappa),
                (self.sigma, other.sigma),
                (self.d, other.d),
                (self.s, other.s),
                (self.vel, other.vel),
                (self.acc, other.acc),
                (self.time, other.time),
                (self.fork_y, other.fork_y),
            )
        )

    def to_string(self) -> str:
        return (
            "State("
            f"x: {self.x}, y: {self.y}, theta: {self.theta}, kappa: {self.kappa}, "
            f"sigma: {self.sigma}, d: {self.d}, s: {self.s}, vel: {self.vel}, "
            f"acc: {self.acc}, time: {self.time}, fork_y: {self.fork_y})"
        )


@dataclass
class Control:
    delta_s: float
    kappa: float
    sigma: float

    def to_str(self) -> str:
        return (
            "Control Segment "
            f"[delta_s: {self.delta_s}, curvature: {self.kappa}, sharpness: {self.sigma}]"
        )


class BaseStateSpace(ABC):
    def __init__(self, discretization: float = 0.1) -> None:
        if discretization <= 0.0:
            raise ValueError("discretization must be positive")
        self._discretization = discretization

    @staticmethod
    def integrate_ode(state: State, control: Control, integration_step: float) -> State:
        next_state = State()
        direction = sgn(control.delta_s)
        if abs(control.sigma) > EPSILON:
            raise NotImplementedError(
                "The Python port currently supports only sigma == 0 controls."
            )

        if abs(control.kappa) > EPSILON:
            next_state.x, next_state.y, next_state.theta = end_of_circular_arc(
                state.x,
                state.y,
                state.theta,
                control.kappa,
                direction,
                integration_step,
            )
            next_state.kappa = control.kappa
        else:
            next_state.x, next_state.y = end_of_straight_line(
                state.x,
                state.y,
                state.theta,
                direction,
                integration_step,
            )
            next_state.theta = state.theta

        next_state.sigma = control.sigma
        next_state.d = direction
        next_state.s = state.s + integration_step
        return next_state

    @abstractmethod
    def get_controls(self, state1: State, state2: State) -> List[Control]:
        raise NotImplementedError

    def get_path(self, state1: State, state2: State) -> List[State]:
        controls = self.get_controls(state1, state2)
        return self.integrate(state1, controls)

    @abstractmethod
    def get_all_controls(self, state1: State, state2: State) -> List[List[Control]]:
        raise NotImplementedError

    def get_all_paths(self, state1: State, state2: State) -> List[List[State]]:
        return [self.integrate(state1, controls) for controls in self.get_all_controls(state1, state2)]

    def interpolate(self, state: State, controls: Iterable[Control], t: float) -> State:
        controls = list(controls)
        if not controls:
            return state

        total_length = sum(abs(control.delta_s) for control in controls)
        target_s = max(0.0, min(1.0, t)) * total_length
        traversed = 0.0
        current = state
        for control in controls:
            segment_length = abs(control.delta_s)
            if target_s - traversed > segment_length:
                current = self.integrate_ode(current, control, segment_length)
                traversed += segment_length
                continue
            return self.integrate_ode(current, control, target_s - traversed)
        return current

    def integrate(
        self,
        state: State,
        controls: Iterable[Control],
        discretization: float | None = None,
    ) -> List[State]:
        controls = list(controls)
        if not controls:
            return []

        step = discretization or self._discretization
        path: List[State] = []
        current = State(
            x=state.x,
            y=state.y,
            theta=state.theta,
            kappa=controls[0].kappa,
            sigma=controls[0].sigma,
            d=sgn(controls[0].delta_s),
            s=state.s,
            vel=state.vel,
            acc=state.acc,
            time=state.time,
            fork_y=state.fork_y,
        )
        path.append(current)

        for control in controls:
            abs_delta_s = abs(control.delta_s)
            n_steps = max(2, math.ceil(abs_delta_s / step))
            step_size = abs_delta_s / n_steps
            for _ in range(n_steps):
                current = self.integrate_ode(current, control, step_size)
                path.append(current)
        return path
