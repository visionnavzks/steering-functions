"""
Base class for all steering state spaces.

Ported from C++ steering_functions/src/base_state_space.cpp.
"""

import math
from abc import ABC, abstractmethod

from steering_functions.state import State, Control
from steering_functions.utilities import (
    sgn,
    get_epsilon,
    end_of_clothoid,
    end_of_circular_arc,
    end_of_straight_line,
)


class BaseStateSpace(ABC):
    """
    Base class for all steering state spaces.

    Provides common path integration, interpolation, and ODE integration.
    Subclasses must implement :meth:`get_controls` and :meth:`get_all_controls`.
    """

    def __init__(self, discretization=0.1):
        assert discretization > 0.0
        self._discretization = discretization

    # ------------------------------------------------------------------
    # ODE integration (single step)
    # ------------------------------------------------------------------
    @staticmethod
    def integrate_ODE(state, control, integration_step):
        """
        Single-step ODE integration for state propagation.

        Handles clothoid (variable curvature), circular arc (constant curvature),
        and straight line segments.
        """
        state_next = State()
        kappa = control.kappa
        sigma = control.sigma
        d = sgn(control.delta_s)

        if abs(sigma) > get_epsilon():
            x_f, y_f, theta_f, kappa_f = end_of_clothoid(
                state.x, state.y, state.theta, state.kappa,
                sigma, d, integration_step,
            )
            state_next.x = x_f
            state_next.y = y_f
            state_next.theta = theta_f
            state_next.kappa = kappa_f
            state_next.sigma = sigma
        else:
            if abs(kappa) > get_epsilon():
                x_f, y_f, theta_f = end_of_circular_arc(
                    state.x, state.y, state.theta,
                    kappa, d, integration_step,
                )
                state_next.x = x_f
                state_next.y = y_f
                state_next.theta = theta_f
                state_next.kappa = kappa
            else:
                x_f, y_f = end_of_straight_line(
                    state.x, state.y, state.theta,
                    d, integration_step,
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
    def get_controls(self, state1, state2):
        """Get controls for the shortest path from *state1* to *state2*."""
        pass

    @abstractmethod
    def get_all_controls(self, state1, state2):
        """Get all possible control sequences between two states."""
        pass

    # ------------------------------------------------------------------
    # Path computation
    # ------------------------------------------------------------------
    def get_path(self, state1, state2, controls_out=None):
        """
        Get path between two states.

        If *controls_out* is provided (a list), it is cleared and filled with
        the control sequence used.
        """
        controls = self.get_controls(state1, state2)
        if controls_out is not None:
            controls_out.clear()
            controls_out.extend(controls)
        return self.integrate(state1, controls)

    # ------------------------------------------------------------------
    # Interpolation
    # ------------------------------------------------------------------
    def interpolate(self, state, controls, t):
        """
        Interpolate state along a path at normalised parameter *t* in [0, 1].
        """
        state_curr = State(
            x=state.x, y=state.y, theta=state.theta,
            kappa=state.kappa, sigma=state.sigma,
            d=state.d, s=state.s, vel=state.vel, acc=state.acc,
            time=state.time, fork_y=state.fork_y,
        )
        state_inter = State(
            x=state.x, y=state.y, theta=state.theta,
            kappa=state.kappa, sigma=state.sigma,
            d=state.d, s=state.s, vel=state.vel, acc=state.acc,
            time=state.time, fork_y=state.fork_y,
        )

        if not controls:
            return state_inter

        # total path length
        s_path = 0.0
        for control in controls:
            s_path += abs(control.delta_s)

        t = min(max(t, 0.0), 1.0)
        s_inter = t * s_path
        s_accumulated = 0.0

        for control in controls:
            abs_delta_s = abs(control.delta_s)
            if s_inter - s_accumulated > abs_delta_s:
                state_curr = BaseStateSpace.integrate_ODE(state_curr, control, abs_delta_s)
                s_accumulated += abs_delta_s
            else:
                state_inter = BaseStateSpace.integrate_ODE(
                    state_curr, control, s_inter - s_accumulated
                )
                break

        return state_inter

    # ------------------------------------------------------------------
    # Integration (discretised path)
    # ------------------------------------------------------------------
    def integrate(self, state, controls, discretization=None):
        """
        Integrate controls to generate a discretised path.
        """
        if discretization is None:
            discretization = self._discretization

        path = []
        if not controls:
            return path

        state_curr = State()
        # initialise first state
        state_curr.x = state.x
        state_curr.y = state.y
        state_curr.theta = state.theta
        state_curr.kappa = controls[0].kappa
        state_curr.sigma = controls[0].sigma
        state_curr.d = sgn(controls[0].delta_s)
        path.append(State(
            x=state_curr.x, y=state_curr.y, theta=state_curr.theta,
            kappa=state_curr.kappa, sigma=state_curr.sigma,
            d=state_curr.d, s=state_curr.s,
        ))

        for control in controls:
            abs_delta_s = abs(control.delta_s)
            s_seg = 0.0

            n = max(2, math.ceil(abs_delta_s / discretization))
            step_size = abs_delta_s / n
            assert step_size > 0

            for _ in range(n):
                s_seg += step_size
                if s_seg > abs_delta_s:
                    integration_step = step_size - (s_seg - abs_delta_s)
                    s_seg = abs_delta_s
                else:
                    integration_step = step_size
                state_next = BaseStateSpace.integrate_ODE(state_curr, control, integration_step)
                path.append(state_next)
                state_curr = state_next

        return path

    # ------------------------------------------------------------------
    # All paths
    # ------------------------------------------------------------------
    def get_all_paths(self, state1, state2):
        """Get all possible discretised paths between two states."""
        all_controls = self.get_all_controls(state1, state2)
        return [self.integrate(state1, cs) for cs in all_controls]
