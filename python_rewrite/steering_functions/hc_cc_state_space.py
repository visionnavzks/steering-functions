"""
HC/CC state space base class.

Ported from C++ hc_cc_state_space.hpp / hc_cc_state_space.cpp.
"""

import math

from steering_functions.base_state_space import BaseStateSpace
from steering_functions.state import State, Control
from steering_functions.hc_cc_circle import HC_CC_Circle_Param
from steering_functions.paths import _build_controls
from steering_functions.utilities import (
    get_epsilon,
    end_of_clothoid,
    point_distance,
)


class HC_CC_StateSpace(BaseStateSpace):
    """
    Base state space for HC- and CC-based steering.

    Computes :class:`HC_CC_Circle_Param` from *kappa* and *sigma*, overrides
    :meth:`get_path` to filter negligible control segments, and provides a
    default :meth:`get_all_controls` that wraps :meth:`get_controls`.
    """

    def __init__(self, kappa, sigma, discretization):
        assert kappa > 0.0 and sigma > 0.0 and discretization > 0.0
        super().__init__(discretization)
        self.kappa_ = kappa
        self.sigma_ = sigma
        self.hc_cc_circle_param_ = HC_CC_Circle_Param()

        # intermediate configuration after first clothoid
        length_min = kappa / sigma
        if length_min > get_epsilon():
            x_i, y_i, theta_i, _kappa_i = end_of_clothoid(
                0.0, 0.0, 0.0, 0.0, sigma, 1, length_min,
            )
        else:
            x_i = 0.0
            y_i = 0.0
            theta_i = 0.0

        # radius
        xc = x_i - math.sin(theta_i) / kappa
        yc = y_i + math.cos(theta_i) / kappa
        radius = point_distance(xc, yc, 0.0, 0.0)

        # mu
        mu = math.atan2(abs(xc), abs(yc))
        sin_mu = math.sin(mu)
        cos_mu = math.cos(mu)

        # delta_min
        delta_min = 0.5 * kappa ** 2 / sigma

        self.hc_cc_circle_param_.set_param(
            kappa, sigma, radius, mu, sin_mu, cos_mu, delta_min,
        )

    # ------------------------------------------------------------------
    # Default get_distance / get_controls via _find_path
    # ------------------------------------------------------------------
    def get_distance(self, state1, state2):
        """Return shortest-path distance between two states."""
        path = self._find_path(state1, state2)
        return path.length

    def get_controls(self, state1, state2):
        """Return controls for the shortest path between two states."""
        path = self._find_path(state1, state2)
        return _build_controls(path, self._CONTROL_TABLE)

    # ------------------------------------------------------------------
    # Path computation (filters negligible segments)
    # ------------------------------------------------------------------
    def get_path(self, state1, state2, controls_out=None):
        """Return path from *state1* to *state2*, filtering negligible segments."""
        unfiltered = self.get_controls(state1, state2)
        controls = [c for c in unfiltered if abs(c.delta_s) > 1e-6]
        if controls_out is not None:
            controls_out.clear()
            controls_out.extend(controls)
        return self.integrate(state1, controls)

    # ------------------------------------------------------------------
    # All controls (single-path wrapper)
    # ------------------------------------------------------------------
    def get_all_controls(self, state1, state2):
        """Return a list containing the single control sequence."""
        return [self.get_controls(state1, state2)]
