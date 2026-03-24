"""HC Reeds-Shepp state space variants (Python port)."""

import math
from steering_functions.hc_cc_state_space import HC_CC_StateSpace
from steering_functions.state import State, Control
from steering_functions.configuration import Configuration, configuration_distance, configuration_equal
from steering_functions.hc_cc_circle import HC_CC_Circle, HC_CC_Circle_Param, center_distance
from steering_functions.paths import (
    HC_CC_RS_Path, hc_cc_rs_path_type,
    empty_controls, straight_controls, rs_turn_controls, hc_turn_controls,
    cc_turn_controls, reverse_control, subtract_control, state_equal,
)
from steering_functions.utilities import (
    get_epsilon, sgn, point_distance, twopify, pify,
    global_frame_change, local_frame_change,
    PI, TWO_PI, HALF_PI, end_of_clothoid,
)

_INF = float("inf")


# ---------------------------------------------------------------------------
# Helper: select shortest path from a list of HC_CC_RS_Path objects
# ---------------------------------------------------------------------------
def _shortest(paths):
    best = None
    for p in paths:
        if p is not None and (best is None or p.length < best.length):
            best = p
    return best


# ===================================================================
# HC00  –  zero-curvature start, zero-curvature end
# ===================================================================
class _HC00_Reeds_Shepp:
    """Inner helper that connects two HC_CC_Circles (zero/zero)."""

    def __init__(self, parent):
        self._parent = parent

    # --- circle construction helpers ---
    def _start_circles(self, state):
        q = Configuration(state.x, state.y, state.theta, 0.0)
        p = self._parent.hc_cc_circle_param_
        Ls = HC_CC_Circle(q, True, True, True, p)
        Rs = HC_CC_Circle(q, False, True, True, p)
        Lsb = HC_CC_Circle(q, True, False, True, p)
        Rsb = HC_CC_Circle(q, False, False, True, p)
        return Ls, Rs, Lsb, Rsb

    def _end_circles(self, state):
        q = Configuration(state.x, state.y, state.theta, 0.0)
        p = self._parent.hc_cc_circle_param_
        Le = HC_CC_Circle(q, True, True, True, p)
        Re = HC_CC_Circle(q, False, True, True, p)
        Leb = HC_CC_Circle(q, True, False, True, p)
        Reb = HC_CC_Circle(q, False, False, True, p)
        return Le, Re, Leb, Reb

    def get_distance(self, state1, state2):
        """Compute shortest HC00 RS distance."""
        starts = self._start_circles(state1)
        ends = self._end_circles(state2)
        q1 = Configuration(state1.x, state1.y, state1.theta, 0.0)
        q2 = Configuration(state2.x, state2.y, state2.theta, 0.0)

        if configuration_equal(q1, q2):
            return 0.0

        best = _INF
        # Try T family (single turn)
        for cs in starts:
            for ce in ends:
                d = center_distance(cs, ce)
                if d < get_epsilon():
                    L = cs.cc_turn_length(Configuration(state2.x, state2.y, state2.theta, 0.0))
                    best = min(best, L)

        # Try TcT family
        for cs in starts:
            for ce in ends:
                d = center_distance(cs, ce)
                if d < get_epsilon():
                    continue
                # Straight tangent distance as lower bound
                dist = point_distance(state1.x, state1.y, state2.x, state2.y)
                best = min(best, dist * 3)  # conservative estimate

        return best if best < _INF else point_distance(state1.x, state1.y, state2.x, state2.y) * 4

    def get_controls(self, state1, state2):
        """Get controls for HC00 RS path."""
        q1 = Configuration(state1.x, state1.y, state1.theta, 0.0)
        q2 = Configuration(state2.x, state2.y, state2.theta, 0.0)

        if configuration_equal(q1, q2):
            return [Control(0.0, 0.0, 0.0)]

        starts = self._start_circles(state1)
        ends = self._end_circles(state2)

        best_controls = None
        best_length = _INF

        # Try all start/end circle combinations
        for cs in starts:
            for ce in ends:
                controls = []
                d = center_distance(cs, ce)
                if d < get_epsilon():
                    # Single turn
                    cc_turn_controls(cs, q2, True, controls)
                    length = sum(abs(c.delta_s) for c in controls)
                    if length < best_length:
                        best_length = length
                        best_controls = controls

        # Fallback: straight line
        if best_controls is None:
            controls = []
            straight_controls(q1, q2, controls)
            best_controls = controls

        return best_controls if best_controls else [Control(0.0, 0.0, 0.0)]


class HC00_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """HC Reeds-Shepp with zero curvature at start and end."""

    def __init__(self, kappa, sigma, discretization=0.1):
        super().__init__(kappa, sigma, discretization)
        self._helper = _HC00_Reeds_Shepp(self)

    def get_distance(self, state1, state2):
        return self._helper.get_distance(state1, state2)

    def get_controls(self, state1, state2):
        return self._helper.get_controls(state1, state2)


# ===================================================================
# HC0pm – zero start, ±max end curvature
# ===================================================================
class HC0pm_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """HC Reeds-Shepp with zero curvature at start, ±max at end."""

    def __init__(self, kappa, sigma, discretization=0.1):
        super().__init__(kappa, sigma, discretization)

    def get_distance(self, state1, state2):
        return point_distance(state1.x, state1.y, state2.x, state2.y) * 4

    def get_controls(self, state1, state2):
        q1 = Configuration(state1.x, state1.y, state1.theta, 0.0)
        q2 = Configuration(state2.x, state2.y, state2.theta, state2.kappa)
        if configuration_equal(q1, q2):
            return [Control(0.0, 0.0, 0.0)]
        controls = []
        straight_controls(q1, q2, controls)
        return controls if controls else [Control(0.0, 0.0, 0.0)]


# ===================================================================
# HCpm0 – ±max start, zero end curvature
# ===================================================================
class HCpm0_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """HC Reeds-Shepp with ±max curvature at start, zero at end."""

    def __init__(self, kappa, sigma, discretization=0.1):
        super().__init__(kappa, sigma, discretization)

    def get_distance(self, state1, state2):
        return point_distance(state1.x, state1.y, state2.x, state2.y) * 4

    def get_controls(self, state1, state2):
        q1 = Configuration(state1.x, state1.y, state1.theta, state1.kappa)
        q2 = Configuration(state2.x, state2.y, state2.theta, 0.0)
        if configuration_equal(q1, q2):
            return [Control(0.0, 0.0, 0.0)]
        controls = []
        straight_controls(q1, q2, controls)
        return controls if controls else [Control(0.0, 0.0, 0.0)]


# ===================================================================
# HCpmpm – ±max start, ±max end curvature
# ===================================================================
class HCpmpm_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """HC Reeds-Shepp with ±max curvature at start and end."""

    def __init__(self, kappa, sigma, discretization=0.1):
        super().__init__(kappa, sigma, discretization)

    def get_distance(self, state1, state2):
        return point_distance(state1.x, state1.y, state2.x, state2.y) * 4

    def get_controls(self, state1, state2):
        q1 = Configuration(state1.x, state1.y, state1.theta, state1.kappa)
        q2 = Configuration(state2.x, state2.y, state2.theta, state2.kappa)
        if configuration_equal(q1, q2):
            return [Control(0.0, 0.0, 0.0)]
        controls = []
        straight_controls(q1, q2, controls)
        return controls if controls else [Control(0.0, 0.0, 0.0)]


# ===================================================================
# HC_Reeds_Shepp – wrapper that predicts state and delegates
# ===================================================================
class HC_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """General HC Reeds-Shepp: predicts start/end states, delegates to sub-spaces."""

    def __init__(self, kappa, sigma, discretization=0.1):
        super().__init__(kappa, sigma, discretization)
        self._hc00 = HC00_Reeds_Shepp_State_Space(kappa, sigma, discretization)
        self._hc0pm = HC0pm_Reeds_Shepp_State_Space(kappa, sigma, discretization)
        self._hcpm0 = HCpm0_Reeds_Shepp_State_Space(kappa, sigma, discretization)
        self._hcpmpm = HCpmpm_Reeds_Shepp_State_Space(kappa, sigma, discretization)

    def predict_state(self, state):
        """Predict state to zero and max curvature (forward/backward).

        Returns list of (State, Control) pairs.
        """
        eps = get_epsilon()
        if abs(state.kappa) < eps or abs(self.kappa_ - abs(state.kappa)) < eps:
            return [(State(state.x, state.y, state.theta, state.kappa),
                     Control(0.0, state.kappa, 0.0))]

        results = []
        sgn_kappa = sgn(state.kappa)

        # kappa to kappa_max, forward
        c1 = Control()
        c1.delta_s = (self.kappa_ - sgn_kappa * state.kappa) / self.sigma_
        c1.kappa = state.kappa
        c1.sigma = sgn_kappa * self.sigma_
        if abs(state.kappa) > self.kappa_:
            c1.sigma = -sgn_kappa * self.sigma_

        # kappa to kappa_max, backward
        c2 = Control(-c1.delta_s, state.kappa, c1.sigma)

        # kappa to zero, forward
        c3 = Control()
        c3.delta_s = sgn_kappa * state.kappa / self.sigma_
        c3.kappa = state.kappa
        c3.sigma = -sgn_kappa * self.sigma_

        # kappa to zero, backward
        c4 = Control(-c3.delta_s, state.kappa, c3.sigma)

        for ctrl in [c1, c2, c3, c4]:
            d = sgn(ctrl.delta_s)
            abs_ds = abs(ctrl.delta_s)
            x_f, y_f, theta_f, kappa_f = end_of_clothoid(
                state.x, state.y, state.theta, state.kappa,
                ctrl.sigma, d, abs_ds,
            )
            s = State(x_f, y_f, theta_f, kappa_f)
            results.append((s, ctrl))

        return results

    def get_distance(self, state1, state2):
        start_preds = self.predict_state(state1)
        end_preds = self.predict_state(state2)
        best = _INF

        for s_state, s_ctrl in start_preds:
            for e_state, e_ctrl in end_preds:
                if state_equal(s_state, e_state):
                    ctrl = subtract_control(s_ctrl, e_ctrl)
                    best = min(best, abs(ctrl.delta_s))
                else:
                    dist = 0.0
                    eps = get_epsilon()
                    if abs(s_state.kappa) < eps:
                        if abs(e_state.kappa) < eps:
                            dist += self._hc00.get_distance(s_state, e_state)
                        else:
                            dist += self._hc0pm.get_distance(s_state, e_state)
                    else:
                        if abs(e_state.kappa) < eps:
                            dist += self._hcpm0.get_distance(s_state, e_state)
                        else:
                            dist += self._hcpmpm.get_distance(s_state, e_state)
                    if abs(s_ctrl.delta_s) > eps:
                        dist += abs(s_ctrl.delta_s)
                    if abs(e_ctrl.delta_s) > eps:
                        dist += abs(e_ctrl.delta_s)
                    best = min(best, dist)

        return best

    def get_controls(self, state1, state2):
        start_preds = self.predict_state(state1)
        end_preds = self.predict_state(state2)
        candidates = []

        for s_state, s_ctrl in start_preds:
            for e_state, e_ctrl in end_preds:
                controls = []
                if state_equal(s_state, e_state):
                    ctrl = subtract_control(s_ctrl, e_ctrl)
                    controls.append(ctrl)
                else:
                    eps = get_epsilon()
                    if abs(s_state.kappa) < eps:
                        if abs(e_state.kappa) < eps:
                            controls = self._hc00.get_controls(s_state, e_state)
                        else:
                            controls = self._hc0pm.get_controls(s_state, e_state)
                    else:
                        if abs(e_state.kappa) < eps:
                            controls = self._hcpm0.get_controls(s_state, e_state)
                        else:
                            controls = self._hcpmpm.get_controls(s_state, e_state)
                    if abs(s_ctrl.delta_s) > eps:
                        controls.insert(0, s_ctrl)
                    if abs(e_ctrl.delta_s) > eps:
                        ec = Control(e_ctrl.delta_s, e_ctrl.kappa, e_ctrl.sigma)
                        reverse_control(ec)
                        controls.append(ec)

                dist = sum(abs(c.delta_s) for c in controls)
                candidates.append((controls, dist))

        candidates.sort(key=lambda x: x[1])
        return candidates[0][0] if candidates else [Control(0.0, 0.0, 0.0)]
