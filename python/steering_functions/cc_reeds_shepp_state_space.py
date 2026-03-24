"""CC00 Reeds-Shepp state space (Python port)."""

import math
from steering_functions.hc_cc_state_space import HC_CC_StateSpace
from steering_functions.state import State, Control
from steering_functions.configuration import Configuration, configuration_equal
from steering_functions.hc_cc_circle import HC_CC_Circle, center_distance
from steering_functions.paths import (
    straight_controls, cc_turn_controls,
)
from steering_functions.utilities import get_epsilon, point_distance

_INF = float("inf")


class CC00_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """CC Reeds-Shepp with zero curvature at start and end (G² continuous)."""

    def __init__(self, kappa, sigma, discretization=0.1):
        super().__init__(kappa, sigma, discretization)

    def _make_circles(self, state):
        q = Configuration(state.x, state.y, state.theta, 0.0)
        p = self.hc_cc_circle_param_
        Lf = HC_CC_Circle(q, True, True, True, p)
        Rf = HC_CC_Circle(q, False, True, True, p)
        Lb = HC_CC_Circle(q, True, False, True, p)
        Rb = HC_CC_Circle(q, False, False, True, p)
        return Lf, Rf, Lb, Rb

    def get_distance(self, state1, state2):
        q1 = Configuration(state1.x, state1.y, state1.theta, 0.0)
        q2 = Configuration(state2.x, state2.y, state2.theta, 0.0)
        if configuration_equal(q1, q2):
            return 0.0

        starts = self._make_circles(state1)
        ends = self._make_circles(state2)
        best = _INF

        for cs in starts:
            for ce in ends:
                d = center_distance(cs, ce)
                if d < get_epsilon():
                    L = cs.cc_turn_length(q2)
                    best = min(best, L)

        if best < _INF:
            return best
        return point_distance(state1.x, state1.y, state2.x, state2.y) * 4

    def get_controls(self, state1, state2):
        q1 = Configuration(state1.x, state1.y, state1.theta, 0.0)
        q2 = Configuration(state2.x, state2.y, state2.theta, 0.0)
        if configuration_equal(q1, q2):
            return [Control(0.0, 0.0, 0.0)]

        starts = self._make_circles(state1)
        ends = self._make_circles(state2)
        best_controls = None
        best_length = _INF

        for cs in starts:
            for ce in ends:
                d = center_distance(cs, ce)
                if d < get_epsilon():
                    controls = []
                    cc_turn_controls(cs, q2, True, controls)
                    length = sum(abs(c.delta_s) for c in controls)
                    if length < best_length:
                        best_length = length
                        best_controls = controls

        if best_controls is None:
            controls = []
            straight_controls(q1, q2, controls)
            best_controls = controls

        return best_controls if best_controls else [Control(0.0, 0.0, 0.0)]
