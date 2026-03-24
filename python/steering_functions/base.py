from __future__ import annotations

import abc
from dataclasses import replace
import math
from typing import Iterable, List

from .types import Control, State


class BaseStateSpace(abc.ABC):
    def __init__(self, discretization: float = 0.1) -> None:
        if discretization <= 0.0:
            raise ValueError("discretization must be > 0")
        self._discretization = discretization

    @staticmethod
    def integrate_ode(state: State, control: Control, integration_step: float) -> State:
        if integration_step < 0.0:
            raise ValueError("integration_step must be >= 0")
        d = 1.0 if control.delta_s >= 0 else -1.0
        x = state.x + math.cos(state.theta) * d * integration_step
        y = state.y + math.sin(state.theta) * d * integration_step
        theta = state.theta + control.kappa * d * integration_step
        return State(
            x=x,
            y=y,
            theta=theta,
            kappa=control.kappa,
            sigma=control.sigma,
            d=d,
            s=state.s + integration_step,
            vel=state.vel,
            acc=state.acc,
            time=state.time,
            fork_y=state.fork_y,
        )

    @abc.abstractmethod
    def get_controls(self, state1: State, state2: State) -> List[Control]:
        raise NotImplementedError

    @abc.abstractmethod
    def get_all_controls(self, state1: State, state2: State) -> List[List[Control]]:
        raise NotImplementedError

    def get_path(self, state1: State, state2: State) -> List[State]:
        return self.integrate(state1, self.get_controls(state1, state2))

    def integrate(
        self,
        state: State,
        controls: Iterable[Control],
        discretization: float | None = None,
    ) -> List[State]:
        controls_list = list(controls)
        if not controls_list:
            return []

        step = self._discretization if discretization is None else discretization
        if step <= 0.0:
            raise ValueError("discretization must be > 0")

        curr = replace(state)
        curr.kappa = controls_list[0].kappa
        curr.sigma = controls_list[0].sigma
        curr.d = 1.0 if controls_list[0].delta_s >= 0 else -1.0
        path = [curr]

        for control in controls_list:
            seg_length = abs(control.delta_s)
            n = max(2, math.ceil(seg_length / step))
            seg_step = seg_length / n
            for _ in range(n):
                curr = self.integrate_ode(curr, control, seg_step)
                path.append(curr)
        return path

    def get_all_paths(self, state1: State, state2: State) -> List[List[State]]:
        return [self.integrate(state1, controls) for controls in self.get_all_controls(state1, state2)]
