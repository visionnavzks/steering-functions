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

"""Base state space with ODE integration for the steering functions library."""

import math
from abc import ABC, abstractmethod
from typing import List

from .state import Control, State
from .utilities import (
    end_of_circular_arc,
    end_of_clothoid,
    end_of_straight_line,
    get_epsilon,
    sgn,
)


class BaseStateSpace(ABC):
    """Abstract base class for all steering state spaces.

    Provides common path integration and interpolation.
    Subclasses must implement :meth:`get_controls` and :meth:`get_all_controls`.
    """

    def __init__(self, discretization: float = 0.1) -> None:
        assert discretization > 0.0
        self._discretization = discretization

    # ------------------------------------------------------------------
    # ODE integration
    # ------------------------------------------------------------------

    @staticmethod
    def integrate_ODE(state: State, control: Control, integration_step: float) -> State:
        """Single-step ODE integration for state propagation.

        Handles clothoid (variable curvature), circular arc (constant
        curvature), and straight-line segments.
        """
        state_next = State()
        kappa = control.kappa
        sigma = control.sigma
        d = sgn(control.delta_s)
        eps = get_epsilon()

        if abs(sigma) > eps:
            x_f, y_f, theta_f, kappa_f = end_of_clothoid(
                state.x,
                state.y,
                state.theta,
                state.kappa,
                sigma,
                d,
                integration_step,
            )
            state_next.x = x_f
            state_next.y = y_f
            state_next.theta = theta_f
            state_next.kappa = kappa_f
            state_next.sigma = sigma
        elif abs(kappa) > eps:
            x_f, y_f, theta_f = end_of_circular_arc(
                state.x, state.y, state.theta, kappa, d, integration_step
            )
            state_next.x = x_f
            state_next.y = y_f
            state_next.theta = theta_f
            state_next.kappa = kappa
        else:
            x_f, y_f = end_of_straight_line(
                state.x, state.y, state.theta, d, integration_step
            )
            state_next.x = x_f
            state_next.y = y_f
            state_next.theta = state.theta

        state_next.d = d
        state_next.s = state.s + integration_step
        return state_next

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    @abstractmethod
    def get_controls(self, state1: State, state2: State) -> List[Control]:
        """Return controls for the shortest path from *state1* to *state2*."""

    @abstractmethod
    def get_all_controls(
        self, state1: State, state2: State
    ) -> List[List[Control]]:
        """Return all valid control sequences between *state1* and *state2*."""

    # ------------------------------------------------------------------
    # Path generation helpers
    # ------------------------------------------------------------------

    def get_path(
        self, state1: State, state2: State, controls: List[Control] | None = None
    ) -> List[State]:
        """Return the discretised path from *state1* to *state2*.

        If *controls* is a list it will be populated with the control
        sequence used to generate the path.
        """
        ctrl = self.get_controls(state1, state2)
        if controls is not None:
            controls.clear()
            controls.extend(ctrl)
        return self.integrate(state1, ctrl)

    def interpolate(
        self, state: State, controls: List[Control], t: float
    ) -> State:
        """Interpolate state along the path at normalised parameter *t* ∈ [0, 1]."""
        state_curr = State(**vars(state))
        state_inter = State(**vars(state))

        if not controls:
            return state_inter

        s_path = sum(abs(c.delta_s) for c in controls)
        t = min(max(t, 0.0), 1.0)
        s_inter = t * s_path
        s_accumulated = 0.0

        for control in controls:
            abs_delta_s = abs(control.delta_s)
            if s_inter - s_accumulated > abs_delta_s:
                state_curr = self.integrate_ODE(state_curr, control, abs_delta_s)
                s_accumulated += abs_delta_s
            else:
                state_inter = self.integrate_ODE(
                    state_curr, control, s_inter - s_accumulated
                )
                break

        return state_inter

    def integrate(
        self,
        state: State,
        controls: List[Control],
        discretization: float | None = None,
    ) -> List[State]:
        """Integrate *controls* from *state* to produce a discretised path.

        Args:
            state: Starting state.
            controls: Control sequence to integrate.
            discretization: Step size.  Defaults to the instance value.

        Returns:
            List of :class:`State` objects representing the path.
        """
        if discretization is None:
            discretization = self._discretization

        path: List[State] = []
        if not controls:
            return path

        state_curr = State()
        state_curr.x = state.x
        state_curr.y = state.y
        state_curr.theta = state.theta
        state_curr.kappa = controls[0].kappa
        state_curr.sigma = controls[0].sigma
        state_curr.d = sgn(controls[0].delta_s)
        path.append(State(**vars(state_curr)))

        for control in controls:
            abs_delta_s = abs(control.delta_s)
            n = max(2, math.ceil(abs_delta_s / discretization))
            step_size = abs_delta_s / n
            s_seg = 0.0
            for i in range(n):
                s_seg += step_size
                if s_seg > abs_delta_s:
                    integration_step = step_size - (s_seg - abs_delta_s)
                    s_seg = abs_delta_s
                else:
                    integration_step = step_size
                state_next = self.integrate_ODE(state_curr, control, integration_step)
                path.append(State(**vars(state_next)))
                state_curr = state_next

        return path

    def get_all_paths(
        self, state1: State, state2: State
    ) -> List[List[State]]:
        """Return all possible discretised paths between *state1* and *state2*."""
        all_controls = self.get_all_controls(state1, state2)
        return [self.integrate(state1, ctrl) for ctrl in all_controls]
