"""HC Reeds-Shepp state space variants (Python port)."""

import math
from steering_functions.hc_cc_state_space import HC_CC_StateSpace
from steering_functions.state import State, Control
from steering_functions.configuration import Configuration, configuration_distance, configuration_equal, configuration_aligned
from steering_functions.hc_cc_circle import HC_CC_Circle, HC_CC_Circle_Param, center_distance, configuration_on_hc_cc_circle
from steering_functions.paths import (
    HC_CC_RS_Path, hc_cc_rs_path_type, nb_hc_cc_rs_paths,
    empty_controls, straight_controls, rs_turn_controls, hc_turn_controls,
    cc_turn_controls, reverse_control, subtract_control, state_equal,
    _build_controls, _evaluate_rs_families,
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
_CC_REGULAR = False


class _HC00_Reeds_Shepp:
    """Inner helper that evaluates all HC00 path families between two circles."""

    def __init__(self, parent):
        self._parent = parent
        self.distance = 0.0
        self.angle = 0.0

    # ---- helpers --------------------------------------------------------
    def _p(self):
        return self._parent.hc_cc_circle_param_

    def _rsp(self):
        return self._parent.rs_circle_param_

    # ---- TT -------------------------------------------------------------
    def TT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return abs(self.distance - 2 * c1.radius) < get_epsilon()

    def TT_tangent_circles(self, c1, c2):
        x = (c1.xc + c2.xc) / 2
        y = (c1.yc + c2.yc) / 2
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        if c1.left:
            theta = angle + HALF_PI - c1.mu if c1.forward else angle + HALF_PI + c1.mu
        else:
            theta = angle - HALF_PI + c1.mu if c1.forward else angle - HALF_PI - c1.mu
        return Configuration(x, y, theta, 0)

    def TT_path(self, c1, c2):
        q = self.TT_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = cstart.cc_turn_length(q) + cend.cc_turn_length(q)
        return length, cstart, cend, q

    # ---- TcT ------------------------------------------------------------
    def TcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return abs(self.distance - 2 * abs(c1.kappa_inv)) < get_epsilon()

    def TcT_tangent_circles(self, c1, c2):
        dist = center_distance(c1, c2)
        delta_x = 0.5 * dist
        delta_y = 0.0
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        if c1.left:
            theta = angle + HALF_PI
            x, y = global_frame_change(c1.xc, c1.yc, angle, delta_x,
                                       delta_y if c1.forward else -delta_y)
        else:
            theta = angle - HALF_PI
            x, y = global_frame_change(c1.xc, c1.yc, angle, delta_x,
                                       -delta_y if c1.forward else delta_y)
        return Configuration(x, y, theta, c1.kappa)

    def TcT_path(self, c1, c2):
        q = self.TcT_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        length = cstart.hc_turn_length(q) + cend.hc_turn_length(q)
        return length, cstart, cend, q

    # ---- TcTcT ----------------------------------------------------------
    def TcTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance <= 4 * abs(c1.kappa_inv)

    def TcTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r = 2 * abs(c1.kappa_inv)
        delta_x = 0.5 * self.distance
        delta_y = math.sqrt(max(0.0, r ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2 = self.TcT_tangent_circles(tgt1, c2)
        q3 = self.TcT_tangent_circles(c1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TcTcT_path(self, c1, c2):
        qa, qb, qc, qd = self.TcTcT_tangent_circles(c1, c2)
        p = self._p()
        rsp = self._rsp()
        middle1 = HC_CC_Circle(qa, not c1.left, not c1.forward, True, rsp)
        middle2 = HC_CC_Circle(qc, not c1.left, not c1.forward, True, rsp)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        l1 = cstart.hc_turn_length(qa) + middle1.rs_turn_length(qb) + cend.hc_turn_length(qb)
        l2 = cstart.hc_turn_length(qc) + middle2.rs_turn_length(qd) + cend.hc_turn_length(qd)
        if l1 < l2:
            return l1, cstart, cend, qa, qb, middle1
        return l2, cstart, cend, qc, qd, middle2

    # ---- TcTT -----------------------------------------------------------
    def TcTT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return ((self.distance <= 2 * c1.radius + 2 * abs(c1.kappa_inv)) and
                (self.distance >= 2 * c1.radius - 2 * abs(c1.kappa_inv)))

    def TcTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * abs(c1.kappa_inv)
        r2 = 2 * c1.radius
        delta_x = (r1 ** 2 + self.distance ** 2 - r2 ** 2) / (2 * self.distance)
        delta_y = math.sqrt(max(0.0, r1 ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2 = self.TT_tangent_circles(tgt1, c2)
        q3 = self.TcT_tangent_circles(c1, tgt2)
        q4 = self.TT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TcTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TcTT_tangent_circles(c1, c2)
        p = self._p()
        middle1 = HC_CC_Circle(qb, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qd, not c1.left, c1.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        l1 = cstart.hc_turn_length(qa) + middle1.hc_turn_length(qa) + cend.cc_turn_length(qb)
        l2 = cstart.hc_turn_length(qc) + middle2.hc_turn_length(qc) + cend.cc_turn_length(qd)
        if l1 < l2:
            return l1, cstart, cend, qa, qb, middle1
        return l2, cstart, cend, qc, qd, middle2

    # ---- TTcT -----------------------------------------------------------
    def TTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return ((self.distance <= 2 * c1.radius + 2 * abs(c1.kappa_inv)) and
                (self.distance >= 2 * c1.radius - 2 * abs(c1.kappa_inv)))

    def TTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * c1.radius
        r2 = 2 * abs(c1.kappa_inv)
        delta_x = (r1 ** 2 + self.distance ** 2 - r2 ** 2) / (2 * self.distance)
        delta_y = math.sqrt(max(0.0, r1 ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        q1 = self.TT_tangent_circles(c1, tgt1)
        q2 = self.TcT_tangent_circles(tgt1, c2)
        q3 = self.TT_tangent_circles(c1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TTcT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTcT_tangent_circles(c1, c2)
        p = self._p()
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        l1 = cstart.cc_turn_length(qa) + middle1.hc_turn_length(qb) + cend.hc_turn_length(qb)
        l2 = cstart.cc_turn_length(qc) + middle2.hc_turn_length(qd) + cend.hc_turn_length(qd)
        if l1 < l2:
            return l1, cstart, cend, qa, qb, middle1
        return l2, cstart, cend, qc, qd, middle2

    # ---- TST (TiST / TeST) ---------------------------------------------
    def TiST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 2 * c1.radius

    def TeST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 2 * c1.radius * c1.sin_mu

    def TST_exists(self, c1, c2):
        return self.TiST_exists(c1, c2) or self.TeST_exists(c1, c2)

    def TiST_tangent_circles(self, c1, c2):
        dist = center_distance(c1, c2)
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        alpha = math.asin(2 * c1.radius * c1.cos_mu / dist)
        dx = c1.radius * c1.sin_mu
        dy = c1.radius * c1.cos_mu
        if c1.left and c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta, 0)
        elif c1.left and not c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif not c1.left and c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        else:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TeST_tangent_circles(self, c1, c2):
        dx = c1.radius * c1.sin_mu
        dy = c1.radius * c1.cos_mu
        theta = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        elif c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta, 0)
        else:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TiST_path(self, c1, c2):
        q1, q2 = self.TiST_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TeST_path(self, c1, c2):
        q1, q2 = self.TeST_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TST_path(self, c1, c2):
        if self.TiST_exists(c1, c2):
            return self.TiST_path(c1, c2)
        if self.TeST_exists(c1, c2):
            return self.TeST_path(c1, c2)
        return None

    # ---- TSTcT (TiSTcT / TeSTcT) ---------------------------------------
    def TiSTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (2 * c1.radius * c1.sin_mu + 2 * abs(c1.kappa_inv)) ** 2 +
            (2 * c1.radius * c1.cos_mu) ** 2)

    def TeSTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= 2 * (abs(c1.kappa_inv) + c1.radius * c1.sin_mu)

    def TSTcT_exists(self, c1, c2):
        return self.TiSTcT_exists(c1, c2) or self.TeSTcT_exists(c1, c2)

    def TiSTcT_path(self, c1, c2):
        theta = self.angle
        delta_y = (4 * c2.radius * c2.cos_mu) / (abs(c2.kappa) * self.distance)
        delta_x = math.sqrt(max(0.0, (2 * c2.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1, q2 = self.TiST_tangent_circles(c1, tgt1)
        q3 = self.TcT_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  ci.hc_turn_length(q3) + cend.hc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3, ci

    def TeSTcT_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c2.kappa_inv)
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1, q2 = self.TeST_tangent_circles(c1, tgt1)
        q3 = self.TcT_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        ci = HC_CC_Circle(q2, c1.left, c1.forward, True, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  ci.hc_turn_length(q3) + cend.hc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3, ci

    def TSTcT_path(self, c1, c2):
        if self.TiSTcT_exists(c1, c2):
            return self.TiSTcT_path(c1, c2)
        if self.TeSTcT_exists(c1, c2):
            return self.TeSTcT_path(c1, c2)
        return None

    # ---- TcTST (TcTiST / TcTeST) ---------------------------------------
    def TcTiST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (2 * c1.radius * c1.sin_mu + 2 * abs(c1.kappa_inv)) ** 2 +
            (2 * c1.radius * c1.cos_mu) ** 2)

    def TcTeST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= 2 * (abs(c1.kappa_inv) + c1.radius * c1.sin_mu)

    def TcTST_exists(self, c1, c2):
        return self.TcTiST_exists(c1, c2) or self.TcTeST_exists(c1, c2)

    def TcTiST_path(self, c1, c2):
        theta = self.angle
        delta_y = (4 * c2.radius * c2.cos_mu) / (abs(c2.kappa) * self.distance)
        delta_x = math.sqrt(max(0.0, (2 * c2.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, not c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TiST_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        length = (cstart.hc_turn_length(q1) + ci.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + cend.cc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3, ci

    def TcTeST_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c2.kappa_inv)
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, c2.left, not c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TeST_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        length = (cstart.hc_turn_length(q1) + ci.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + cend.cc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3, ci

    def TcTST_path(self, c1, c2):
        if self.TcTiST_exists(c1, c2):
            return self.TcTiST_path(c1, c2)
        if self.TcTeST_exists(c1, c2):
            return self.TcTeST_path(c1, c2)
        return None

    # ---- TcTSTcT (TcTiSTcT / TcTeSTcT) ---------------------------------
    def TcTiSTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= math.sqrt(
            (2 * c1.radius) ** 2 +
            16 * c1.radius * c1.sin_mu * abs(c1.kappa_inv) +
            (4 * c1.kappa_inv) ** 2)

    def TcTeSTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 4 * abs(c1.kappa_inv) + 2 * c1.radius * c1.sin_mu

    def TcTSTcT_exists(self, c1, c2):
        return self.TcTiSTcT_exists(c1, c2) or self.TcTeSTcT_exists(c1, c2)

    def TcTiSTcT_path(self, c1, c2):
        theta = self.angle
        delta_y = (4 * c1.radius * c1.cos_mu) / (self.distance * abs(c1.kappa))
        delta_x = math.sqrt(max(0.0, (2 * c1.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TiST_tangent_circles(tgt1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        ci1 = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        ci2 = HC_CC_Circle(q3, not c2.left, c2.forward, True, p)
        length = (cstart.hc_turn_length(q1) + ci1.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + ci2.hc_turn_length(q4) +
                  cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci1, ci2

    def TcTeSTcT_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c1.kappa_inv)
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TeST_tangent_circles(tgt1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        ci1 = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        ci2 = HC_CC_Circle(q3, not c2.left, c2.forward, True, p)
        length = (cstart.hc_turn_length(q1) + ci1.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + ci2.hc_turn_length(q4) +
                  cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci1, ci2

    def TcTSTcT_path(self, c1, c2):
        if self.TcTiSTcT_exists(c1, c2):
            return self.TcTiSTcT_path(c1, c2)
        if self.TcTeSTcT_exists(c1, c2):
            return self.TcTeSTcT_path(c1, c2)
        return None

    # ---- TTcTT ----------------------------------------------------------
    def TTcTT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance <= 4 * c1.radius + 2 * abs(c1.kappa_inv)

    def TTcTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * abs(c1.kappa_inv)
        r2 = 2 * c1.radius
        p = self._p()
        if self.distance < 4 * c1.radius - 2 * abs(c1.kappa_inv):
            delta_x = (self.distance + r1) / 2
        else:
            delta_x = (self.distance - r1) / 2
        delta_y = math.sqrt(max(0.0, r2 ** 2 - delta_x ** 2))
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, not c2.forward, c2.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt3 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt4 = HC_CC_Circle(x, y, not c2.left, not c2.forward, c2.regular, p)
        q1 = self.TT_tangent_circles(c1, tgt1)
        q2 = self.TcT_tangent_circles(tgt1, tgt2)
        q3 = self.TT_tangent_circles(tgt2, c2)
        q4 = self.TT_tangent_circles(c1, tgt3)
        q5 = self.TcT_tangent_circles(tgt3, tgt4)
        q6 = self.TT_tangent_circles(tgt4, c2)
        return q1, q2, q3, q4, q5, q6

    def TTcTT_path(self, c1, c2):
        qa, qb, qc, qd, qe, qf = self.TTcTT_tangent_circles(c1, c2)
        p = self._p()
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qc, not c2.left, c2.forward, True, p)
        middle3 = HC_CC_Circle(qd, not c1.left, c1.forward, True, p)
        middle4 = HC_CC_Circle(qf, not c2.left, c2.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        l1 = (cstart.cc_turn_length(qa) + middle1.hc_turn_length(qb) +
              middle2.hc_turn_length(qb) + cend.cc_turn_length(qc))
        l2 = (cstart.cc_turn_length(qd) + middle3.hc_turn_length(qe) +
              middle4.hc_turn_length(qe) + cend.cc_turn_length(qf))
        if l1 < l2:
            return l1, cstart, cend, qa, qb, qc, middle1, middle2
        return l2, cstart, cend, qd, qe, qf, middle3, middle4

    # ---- TcTTcT ---------------------------------------------------------
    def TcTTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return ((self.distance <= 4 * abs(c1.kappa_inv) + 2 * c1.radius) and
                (self.distance >= 4 * abs(c1.kappa_inv) - 2 * c1.radius))

    def TcTTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * abs(c1.kappa_inv)
        r2 = c1.radius
        delta_x = (r1 ** 2 + (self.distance / 2) ** 2 - r2 ** 2) / self.distance
        delta_y = math.sqrt(max(0.0, r1 ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt3 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt4 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2 = self.TT_tangent_circles(tgt1, tgt2)
        q3 = self.TcT_tangent_circles(tgt2, c2)
        q4 = self.TcT_tangent_circles(c1, tgt3)
        q5 = self.TT_tangent_circles(tgt3, tgt4)
        q6 = self.TcT_tangent_circles(tgt4, c2)
        return q1, q2, q3, q4, q5, q6

    def TcTTcT_path(self, c1, c2):
        qa, qb, qc, qd, qe, qf = self.TcTTcT_tangent_circles(c1, c2)
        p = self._p()
        middle1 = HC_CC_Circle(qb, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qb, c1.left, not c1.forward, True, p)
        middle3 = HC_CC_Circle(qe, not c1.left, c1.forward, True, p)
        middle4 = HC_CC_Circle(qe, c1.left, not c1.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        l1 = (cstart.hc_turn_length(qa) + middle1.hc_turn_length(qa) +
              middle2.hc_turn_length(qc) + cend.hc_turn_length(qc))
        l2 = (cstart.hc_turn_length(qd) + middle3.hc_turn_length(qd) +
              middle4.hc_turn_length(qf) + cend.hc_turn_length(qf))
        if l1 < l2:
            return l1, cstart, cend, qa, qb, qc, middle1, middle2
        return l2, cstart, cend, qd, qe, qf, middle3, middle4

    # ---- TTT ------------------------------------------------------------
    def TTT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance <= 4 * c1.radius

    def TTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r = 2 * c1.radius
        delta_x = 0.5 * self.distance
        delta_y = math.sqrt(max(0.0, r ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        q1 = self.TT_tangent_circles(c1, tgt1)
        q2 = self.TT_tangent_circles(tgt1, c2)
        q3 = self.TT_tangent_circles(c1, tgt2)
        q4 = self.TT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTT_tangent_circles(c1, c2)
        p = self._p()
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        l1 = (cstart.cc_turn_length(qa) + middle1.cc_turn_length(qb) +
              cend.cc_turn_length(qb))
        l2 = (cstart.cc_turn_length(qc) + middle2.cc_turn_length(qd) +
              cend.cc_turn_length(qd))
        if l1 < l2:
            return l1, cstart, cend, qa, qb, middle1
        return l2, cstart, cend, qc, qd, middle2

    # ---- TcST (TciST / TceST) ------------------------------------------
    def TciST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (c1.radius * c1.sin_mu) ** 2 +
            (c1.radius * c1.cos_mu + abs(c1.kappa_inv)) ** 2)

    def TceST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (c1.radius * c1.sin_mu) ** 2 +
            (c1.radius * c1.cos_mu - abs(c1.kappa_inv)) ** 2)

    def TcST_exists(self, c1, c2):
        return self.TciST_exists(c1, c2) or self.TceST_exists(c1, c2)

    def TciST_path(self, c1, c2):
        alpha = math.asin((c1.radius * c1.cos_mu + abs(c1.kappa_inv)) / self.distance)
        dx1, dy1 = 0.0, abs(c1.kappa_inv)
        dx2 = c1.radius * c1.sin_mu
        dy2 = c1.radius * c1.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.hc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TceST_path(self, c1, c2):
        alpha = math.asin((c1.radius * c1.cos_mu - abs(c1.kappa_inv)) / self.distance)
        dx1, dy1 = 0.0, abs(c1.kappa_inv)
        dx2 = c1.radius * c1.sin_mu
        dy2 = c1.radius * c1.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.hc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TcST_path(self, c1, c2):
        if self.TciST_exists(c1, c2):
            return self.TciST_path(c1, c2)
        if self.TceST_exists(c1, c2):
            return self.TceST_path(c1, c2)
        return None

    # ---- TScT (TiScT / TeScT) ------------------------------------------
    def TiScT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (c1.radius * c1.sin_mu) ** 2 +
            (c1.radius * c1.cos_mu + abs(c1.kappa_inv)) ** 2)

    def TeScT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (c1.radius * c1.sin_mu) ** 2 +
            (c1.radius * c1.cos_mu - abs(c1.kappa_inv)) ** 2)

    def TScT_exists(self, c1, c2):
        return self.TiScT_exists(c1, c2) or self.TeScT_exists(c1, c2)

    def TiScT_path(self, c1, c2):
        alpha = math.asin((c1.radius * c1.cos_mu + abs(c1.kappa_inv)) / self.distance)
        dx1 = c1.radius * c1.sin_mu
        dy1 = c1.radius * c1.cos_mu
        dx2, dy2 = 0.0, abs(c1.kappa_inv)
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q2 = Configuration(x, y, theta, c2.kappa)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.hc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TeScT_path(self, c1, c2):
        alpha = math.asin((c1.radius * c1.cos_mu - abs(c1.kappa_inv)) / self.distance)
        dx1 = c1.radius * c1.sin_mu
        dy1 = c1.radius * c1.cos_mu
        dx2, dy2 = 0.0, abs(c1.kappa_inv)
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q2 = Configuration(x, y, theta, c2.kappa)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.hc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TScT_path(self, c1, c2):
        if self.TiScT_exists(c1, c2):
            return self.TiScT_path(c1, c2)
        if self.TeScT_exists(c1, c2):
            return self.TeScT_path(c1, c2)
        return None

    # ---- TcScT (TciScT / TceScT) ---------------------------------------
    def TciScT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance > 2 * abs(c1.kappa_inv)

    def TceScT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= get_epsilon()

    def TcScT_exists(self, c1, c2):
        return self.TciScT_exists(c1, c2) or self.TceScT_exists(c1, c2)

    def TciScT_path(self, c1, c2):
        alpha = math.asin(2 / (abs(c1.kappa) * self.distance))
        dx, dy = 0.0, abs(c1.kappa_inv)
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        else:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        length = (cstart.hc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.hc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TceScT_path(self, c1, c2):
        theta = self.angle
        dx, dy = 0.0, abs(c1.kappa_inv)
        p = self._p()
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        else:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        length = (cstart.hc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.hc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TcScT_path(self, c1, c2):
        if self.TciScT_exists(c1, c2):
            return self.TciScT_path(c1, c2)
        if self.TceScT_exists(c1, c2):
            return self.TceScT_path(c1, c2)
        return None


class HC00_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """HC Reeds-Shepp with zero curvature at start and end."""

    _CONTROL_TABLE = {
        hc_cc_rs_path_type.E: [('empty',)],
        hc_cc_rs_path_type.S: [('straight', 'start', 'end')],
        hc_cc_rs_path_type.T: [('cc', 'cstart', 'end', True)],
        hc_cc_rs_path_type.TT: [('cc', 'cstart', 'qi1', True), ('cc', 'cend', 'qi1', False)],
        hc_cc_rs_path_type.TcT: [('hc', 'cstart', 'qi1', True), ('hc', 'cend', 'qi1', False)],
        hc_cc_rs_path_type.TcTcT: [('hc', 'cstart', 'qi1', True), ('rs', 'ci1', 'qi2', True), ('hc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TcTT: [('hc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('cc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TTcT: [('cc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi2', True), ('hc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TST: [('cc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('cc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TSTcT: [('cc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('hc', 'ci1', 'qi3', True), ('hc', 'cend', 'qi3', False)],
        hc_cc_rs_path_type.TcTST: [('hc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('straight', 'qi2', 'qi3'), ('cc', 'cend', 'qi3', False)],
        hc_cc_rs_path_type.TcTSTcT: [('hc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'ci2', 'qi4', True), ('hc', 'cend', 'qi4', False)],
        hc_cc_rs_path_type.TTcTT: [('cc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi2', True), ('hc', 'ci2', 'qi2', False), ('cc', 'cend', 'qi3', False)],
        hc_cc_rs_path_type.TcTTcT: [('hc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('hc', 'ci2', 'qi3', True), ('hc', 'cend', 'qi3', False)],
        hc_cc_rs_path_type.TTT: [('cc', 'cstart', 'qi1', True), ('cc', 'ci1', 'qi2', True), ('cc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TcST: [('hc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('cc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TScT: [('cc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('hc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TcScT: [('hc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('hc', 'cend', 'qi2', False)],
    }

    _FAMILY_REGISTRY = [
        ('TT',      hc_cc_rs_path_type.TT,      ('length', 'cstart', 'cend', 'qi1'), False),
        ('TcT',     hc_cc_rs_path_type.TcT,      ('length', 'cstart', 'cend', 'qi1'), False),
        ('TcTcT',   hc_cc_rs_path_type.TcTcT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TcTT',    hc_cc_rs_path_type.TcTT,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TTcT',    hc_cc_rs_path_type.TTcT,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TST',     hc_cc_rs_path_type.TST,      ('length', 'cstart', 'cend', 'qi1', 'qi2'), True),
        ('TSTcT',   hc_cc_rs_path_type.TSTcT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1'), True),
        ('TcTST',   hc_cc_rs_path_type.TcTST,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1'), True),
        ('TcTSTcT', hc_cc_rs_path_type.TcTSTcT,  ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'qi4', 'ci1', 'ci2'), True),
        ('TTcTT',   hc_cc_rs_path_type.TTcTT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1', 'ci2'), False),
        ('TcTTcT',  hc_cc_rs_path_type.TcTTcT,   ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1', 'ci2'), False),
        ('TTT',     hc_cc_rs_path_type.TTT,      ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TcST',    hc_cc_rs_path_type.TcST,     ('length', 'cstart', 'cend', 'qi1', 'qi2'), True),
        ('TScT',    hc_cc_rs_path_type.TScT,     ('length', 'cstart', 'cend', 'qi1', 'qi2'), True),
        ('TcScT',   hc_cc_rs_path_type.TcScT,    ('length', 'cstart', 'cend', 'qi1', 'qi2'), True),
    ]

    def __init__(self, kappa, sigma, discretization=0.1):
        super().__init__(kappa, sigma, discretization)
        self._helper = _HC00_Reeds_Shepp(self)
        self.rs_circle_param_ = HC_CC_Circle_Param()
        import sys
        self.rs_circle_param_.set_param(
            self.kappa_, sys.float_info.max, 1.0 / self.kappa_, 0.0, 0.0, 1.0, 0.0)

    def hc00_circles_rs_path(self, c1, c2):
        """Evaluate all families and return the shortest HC_CC_RS_Path."""
        tp = hc_cc_rs_path_type
        h = self._helper
        h.distance = center_distance(c1, c2)
        h.angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        p = self.hc_cc_circle_param_

        N = nb_hc_cc_rs_paths
        length = [_INF] * N
        qi1 = [None] * N
        qi2 = [None] * N
        qi3 = [None] * N
        qi4 = [None] * N
        cstart = [None] * N
        cend = [None] * N
        ci1 = [None] * N
        ci2 = [None] * N

        skip = False

        # case E
        if configuration_equal(c1.start, c2.start):
            length[tp.E] = 0
            skip = True

        # case S
        if not skip and configuration_aligned(c1.start, c2.start):
            length[tp.S] = configuration_distance(c1.start, c2.start)
            skip = True
        if not skip and configuration_aligned(c2.start, c1.start):
            length[tp.S] = configuration_distance(c2.start, c1.start)
            skip = True

        # case T
        if not skip and configuration_on_hc_cc_circle(c1, c2.start):
            cstart[tp.T] = HC_CC_Circle(
                c1.start, c1.left, c1.forward, _CC_REGULAR, p)
            length[tp.T] = cstart[tp.T].cc_turn_length(c2.start)
            skip = True

        if not skip:
            length, qi1, qi2, qi3, qi4, cstart, cend, ci1, ci2 = \
                _evaluate_rs_families(h, c1, c2, self._FAMILY_REGISTRY)

        # select shortest
        best = min(range(N), key=lambda i: length[i])
        return HC_CC_RS_Path(
            c1.start, c2.start, hc_cc_rs_path_type(best),
            self.kappa_, self.sigma_,
            qi1[best], qi2[best], qi3[best], qi4[best],
            cstart[best], cend[best], ci1[best], ci2[best],
            length[best])

    def _find_path(self, state1, state2):
        """Compute shortest path between two states (4×4 circle combos)."""
        start = Configuration(state1.x, state1.y, state1.theta, 0.0)
        end = Configuration(state2.x, state2.y, state2.theta, 0.0)
        p = self.hc_cc_circle_param_
        start_circles = [
            HC_CC_Circle(start, True, True, True, p),
            HC_CC_Circle(start, False, True, True, p),
            HC_CC_Circle(start, True, False, True, p),
            HC_CC_Circle(start, False, False, True, p),
        ]
        end_circles = [
            HC_CC_Circle(end, True, True, True, p),
            HC_CC_Circle(end, False, True, True, p),
            HC_CC_Circle(end, True, False, True, p),
            HC_CC_Circle(end, False, False, True, p),
        ]
        best_path = None
        for sc in start_circles:
            for ec in end_circles:
                path = self.hc00_circles_rs_path(sc, ec)
                if path is not None and (best_path is None or path.length < best_path.length):
                    best_path = path
        return best_path


# ===================================================================
# HC0pm – zero start, ±max end curvature
# ===================================================================
_HC_REGULAR = False


class _HC0pm_Reeds_Shepp(_HC00_Reeds_Shepp):
    """Inner helper for HC0pm path families (zero start, ±max end curvature).

    Inherits all ``_exists`` and tangent-circle methods from
    :class:`_HC00_Reeds_Shepp`.  Only ``_path`` methods and
    ``TcTcT_tangent_circles`` are overridden to account for RS end circles.
    """

    # ---- TT -------------------------------------------------------------
    def TT_path(self, c1, c2):
        q1 = self.TT_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(q1, c2.left, not c2.forward, _HC_REGULAR, p)
        q2 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = cstart.cc_turn_length(q1) + cend.hc_turn_length(q2)
        return length, cstart, cend, q1, q2

    # ---- TcT ------------------------------------------------------------
    def TcT_path(self, c1, c2):
        q = self.TcT_tangent_circles(c1, c2)
        p = self._p()
        rsp = self._rsp()
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        length = cstart.hc_turn_length(q) + cend.rs_turn_length(q)
        return length, cstart, cend, q

    # ---- TcTcT ----------------------------------------------------------
    def TcTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r = 2 * abs(c1.kappa_inv)
        delta_x = 0.5 * self.distance
        delta_y = math.sqrt(max(0.0, r ** 2 - delta_x ** 2))
        rsp = self._rsp()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, True, rsp)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, True, rsp)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2 = self.TcT_tangent_circles(tgt1, c2)
        q3 = self.TcT_tangent_circles(c1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TcTcT_path(self, c1, c2):
        qa, qb, qc, qd = self.TcTcT_tangent_circles(c1, c2)
        p = self._p()
        rsp = self._rsp()
        middle1 = HC_CC_Circle(qa, not c2.left, c2.forward, True, rsp)
        middle2 = HC_CC_Circle(qc, not c2.left, c2.forward, True, rsp)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        l1 = cstart.hc_turn_length(qa) + middle1.rs_turn_length(qb) + cend.rs_turn_length(qb)
        l2 = cstart.hc_turn_length(qc) + middle2.rs_turn_length(qd) + cend.rs_turn_length(qd)
        if l1 < l2:
            return l1, cstart, cend, qa, qb, middle1
        return l2, cstart, cend, qc, qd, middle2

    # ---- TcTT -----------------------------------------------------------
    def TcTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TcTT_tangent_circles(c1, c2)
        p = self._p()
        middle1 = HC_CC_Circle(qb, not c1.left, c1.forward, True, p)
        end1 = HC_CC_Circle(qb, c2.left, not c2.forward, _HC_REGULAR, p)
        middle2 = HC_CC_Circle(qd, not c1.left, c1.forward, True, p)
        end2 = HC_CC_Circle(qd, c2.left, not c2.forward, _HC_REGULAR, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        q2 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        l1 = cstart.hc_turn_length(qa) + middle1.hc_turn_length(qa) + end1.hc_turn_length(q2)
        l2 = cstart.hc_turn_length(qc) + middle2.hc_turn_length(qc) + end2.hc_turn_length(q2)
        if l1 < l2:
            return l1, cstart, end1, qa, q2, middle1
        return l2, cstart, end2, qc, q2, middle2

    # ---- TTcT -----------------------------------------------------------
    def TTcT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTcT_tangent_circles(c1, c2)
        p = self._p()
        rsp = self._rsp()
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        l1 = cstart.cc_turn_length(qa) + middle1.hc_turn_length(qb) + cend.rs_turn_length(qb)
        l2 = cstart.cc_turn_length(qc) + middle2.hc_turn_length(qd) + cend.rs_turn_length(qd)
        if l1 < l2:
            return l1, cstart, cend, qa, qb, middle1
        return l2, cstart, cend, qc, qd, middle2

    # ---- TST (TiST / TeST) ---------------------------------------------
    def TiST_path(self, c1, c2):
        q1, q2 = self.TiST_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(q2, c2.left, not c2.forward, _HC_REGULAR, p)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.hc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    def TeST_path(self, c1, c2):
        q1, q2 = self.TeST_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(q2, c2.left, not c2.forward, _HC_REGULAR, p)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.hc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    # ---- TSTcT (TiSTcT / TeSTcT) ---------------------------------------
    def TiSTcT_path(self, c1, c2):
        theta = self.angle
        delta_y = (4 * c1.radius * c1.cos_mu) / (abs(c1.kappa) * self.distance)
        delta_x = math.sqrt(max(0.0, (2 * c1.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        rsp = self._rsp()
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1, q2 = self.TiST_tangent_circles(c1, tgt1)
        q3 = self.TcT_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  ci.hc_turn_length(q3) + cend.rs_turn_length(q3))
        return length, cstart, cend, q1, q2, q3, ci

    def TeSTcT_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c2.kappa_inv)
        delta_y = 0.0
        p = self._p()
        rsp = self._rsp()
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1, q2 = self.TeST_tangent_circles(c1, tgt1)
        q3 = self.TcT_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        ci = HC_CC_Circle(q2, c1.left, c1.forward, True, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  ci.hc_turn_length(q3) + cend.rs_turn_length(q3))
        return length, cstart, cend, q1, q2, q3, ci

    # ---- TcTST (TcTiST / TcTeST) ---------------------------------------
    def TcTiST_path(self, c1, c2):
        theta = self.angle
        delta_y = (4 * c1.radius * c1.cos_mu) / (abs(c1.kappa) * self.distance)
        delta_x = math.sqrt(max(0.0, (2 * c1.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TiST_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(q3, c2.left, not c2.forward, _HC_REGULAR, p)
        q4 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        length = (cstart.hc_turn_length(q1) + ci.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci

    def TcTeST_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c2.kappa_inv)
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TeST_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(q3, c2.left, not c2.forward, _HC_REGULAR, p)
        q4 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        length = (cstart.hc_turn_length(q1) + ci.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci

    # ---- TcTSTcT (TcTiSTcT / TcTeSTcT) ---------------------------------
    def TcTiSTcT_path(self, c1, c2):
        theta = self.angle
        delta_y = (4 * c1.radius * c1.cos_mu) / (self.distance * abs(c1.kappa))
        delta_x = math.sqrt(max(0.0, (2 * c1.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        rsp = self._rsp()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TiST_tangent_circles(tgt1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        ci1 = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        ci2 = HC_CC_Circle(q3, not c2.left, c2.forward, True, p)
        length = (cstart.hc_turn_length(q1) + ci1.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + ci2.hc_turn_length(q4) +
                  cend.rs_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci1, ci2

    def TcTeSTcT_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c1.kappa_inv)
        delta_y = 0.0
        p = self._p()
        rsp = self._rsp()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TeST_tangent_circles(tgt1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        ci1 = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        ci2 = HC_CC_Circle(q3, not c2.left, c2.forward, True, p)
        length = (cstart.hc_turn_length(q1) + ci1.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + ci2.hc_turn_length(q4) +
                  cend.rs_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci1, ci2

    # ---- TTcTT ----------------------------------------------------------
    def TTcTT_path(self, c1, c2):
        qa, qb, qc, qd, qe, qf = self.TTcTT_tangent_circles(c1, c2)
        p = self._p()
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qc, not c2.left, c2.forward, True, p)
        end1 = HC_CC_Circle(qc, c2.left, not c2.forward, _HC_REGULAR, p)
        middle3 = HC_CC_Circle(qd, not c1.left, c1.forward, True, p)
        middle4 = HC_CC_Circle(qf, not c2.left, c2.forward, True, p)
        end2 = HC_CC_Circle(qf, c2.left, not c2.forward, _HC_REGULAR, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        l1 = (cstart.cc_turn_length(qa) + middle1.hc_turn_length(qb) +
              middle2.hc_turn_length(qb) + end1.hc_turn_length(q3))
        l2 = (cstart.cc_turn_length(qd) + middle3.hc_turn_length(qe) +
              middle4.hc_turn_length(qe) + end2.hc_turn_length(q3))
        if l1 < l2:
            return l1, cstart, end1, qa, qb, q3, middle1, middle2
        return l2, cstart, end2, qd, qe, q3, middle3, middle4

    # ---- TcTTcT ---------------------------------------------------------
    def TcTTcT_path(self, c1, c2):
        qa, qb, qc, qd, qe, qf = self.TcTTcT_tangent_circles(c1, c2)
        p = self._p()
        rsp = self._rsp()
        middle1 = HC_CC_Circle(qb, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qb, c1.left, not c1.forward, True, p)
        middle3 = HC_CC_Circle(qe, not c1.left, c1.forward, True, p)
        middle4 = HC_CC_Circle(qe, c1.left, not c1.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        l1 = (cstart.hc_turn_length(qa) + middle1.hc_turn_length(qa) +
              middle2.hc_turn_length(qc) + cend.rs_turn_length(qc))
        l2 = (cstart.hc_turn_length(qd) + middle3.hc_turn_length(qd) +
              middle4.hc_turn_length(qf) + cend.rs_turn_length(qf))
        if l1 < l2:
            return l1, cstart, cend, qa, qc, middle1, middle2
        return l2, cstart, cend, qd, qf, middle3, middle4

    # ---- TTT ------------------------------------------------------------
    def TTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTT_tangent_circles(c1, c2)
        p = self._p()
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, _CC_REGULAR, p)
        end1 = HC_CC_Circle(qb, c2.left, not c2.forward, _HC_REGULAR, p)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, _CC_REGULAR, p)
        end2 = HC_CC_Circle(qd, c2.left, not c2.forward, _HC_REGULAR, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        l1 = cstart.cc_turn_length(qa) + middle1.cc_turn_length(qb) + end1.hc_turn_length(q3)
        l2 = cstart.cc_turn_length(qc) + middle2.cc_turn_length(qd) + end2.hc_turn_length(q3)
        if l1 < l2:
            return l1, cstart, end1, qa, qb, q3, middle1
        return l2, cstart, end2, qc, qd, q3, middle2

    # ---- TcST (TciST / TceST) ------------------------------------------
    def TciST_path(self, c1, c2):
        alpha = math.asin((c1.radius * c1.cos_mu + abs(c1.kappa_inv)) / self.distance)
        dx1, dy1 = 0.0, abs(c1.kappa_inv)
        dx2 = c1.radius * c1.sin_mu
        dy2 = c1.radius * c1.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(q2, c2.left, not c2.forward, _HC_REGULAR, p)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = (cstart.hc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.hc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    def TceST_path(self, c1, c2):
        alpha = math.asin((c1.radius * c1.cos_mu - abs(c1.kappa_inv)) / self.distance)
        dx1, dy1 = 0.0, abs(c1.kappa_inv)
        dx2 = c1.radius * c1.sin_mu
        dy2 = c1.radius * c1.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(q2, c2.left, not c2.forward, _HC_REGULAR, p)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = (cstart.hc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.hc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    # ---- TScT (TiScT / TeScT) ------------------------------------------
    def TiScT_path(self, c1, c2):
        alpha = math.asin((c1.radius * c1.cos_mu + abs(c1.kappa_inv)) / self.distance)
        dx1 = c1.radius * c1.sin_mu
        dy1 = c1.radius * c1.cos_mu
        dx2, dy2 = 0.0, abs(c1.kappa_inv)
        p = self._p()
        rsp = self._rsp()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q2 = Configuration(x, y, theta, c2.kappa)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.rs_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TeScT_path(self, c1, c2):
        alpha = math.asin((c1.radius * c1.cos_mu - abs(c1.kappa_inv)) / self.distance)
        dx1 = c1.radius * c1.sin_mu
        dy1 = c1.radius * c1.cos_mu
        dx2, dy2 = 0.0, abs(c1.kappa_inv)
        p = self._p()
        rsp = self._rsp()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q2 = Configuration(x, y, theta, c2.kappa)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.rs_turn_length(q2))
        return length, cstart, cend, q1, q2

    # ---- TcScT (TciScT / TceScT) ---------------------------------------
    def TciScT_path(self, c1, c2):
        alpha = math.asin(2 / (abs(c1.kappa) * self.distance))
        dx, dy = 0.0, abs(c1.kappa_inv)
        p = self._p()
        rsp = self._rsp()
        if c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        else:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        length = (cstart.hc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.rs_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TceScT_path(self, c1, c2):
        theta = self.angle
        dx, dy = 0.0, abs(c1.kappa_inv)
        p = self._p()
        rsp = self._rsp()
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        else:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        length = (cstart.hc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.rs_turn_length(q2))
        return length, cstart, cend, q1, q2


class HC0pm_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """HC Reeds-Shepp with zero curvature at start, ±max at end."""

    _CONTROL_TABLE = {
        hc_cc_rs_path_type.E: [('empty',)],
        hc_cc_rs_path_type.T: [('hc', 'cstart', 'end', True)],
        hc_cc_rs_path_type.TT: [('cc', 'cstart', 'qi1', True), ('hc', 'cend', 'qi2', True)],
        hc_cc_rs_path_type.TcT: [('hc', 'cstart', 'qi1', True), ('rs', 'cend', 'qi1', False)],
        hc_cc_rs_path_type.TcTcT: [('hc', 'cstart', 'qi1', True), ('rs', 'ci1', 'qi2', True), ('rs', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TcTT: [('hc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('hc', 'cend', 'qi2', True)],
        hc_cc_rs_path_type.TTcT: [('cc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi2', True), ('rs', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TST: [('cc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('hc', 'cend', 'qi3', True)],
        hc_cc_rs_path_type.TSTcT: [('cc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('hc', 'ci1', 'qi3', True), ('rs', 'cend', 'qi3', False)],
        hc_cc_rs_path_type.TcTST: [('hc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'cend', 'qi4', True)],
        hc_cc_rs_path_type.TcTSTcT: [('hc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'ci2', 'qi4', True), ('rs', 'cend', 'qi4', False)],
        hc_cc_rs_path_type.TTcTT: [('cc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi2', True), ('hc', 'ci2', 'qi2', False), ('hc', 'cend', 'qi3', True)],
        hc_cc_rs_path_type.TcTTcT: [('hc', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('hc', 'ci2', 'qi2', True), ('rs', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TTT: [('cc', 'cstart', 'qi1', True), ('cc', 'ci1', 'qi2', True), ('hc', 'cend', 'qi3', True)],
        hc_cc_rs_path_type.TcST: [('hc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('hc', 'cend', 'qi3', True)],
        hc_cc_rs_path_type.TScT: [('cc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('rs', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TcScT: [('hc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('rs', 'cend', 'qi2', False)],
    }

    _FAMILY_REGISTRY = [
        ('TT',      hc_cc_rs_path_type.TT,      ('length', 'cstart', 'cend', 'qi1', 'qi2'), False),
        ('TcT',     hc_cc_rs_path_type.TcT,      ('length', 'cstart', 'cend', 'qi1'), False),
        ('TcTcT',   hc_cc_rs_path_type.TcTcT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TcTT',    hc_cc_rs_path_type.TcTT,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TTcT',    hc_cc_rs_path_type.TTcT,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TST',     hc_cc_rs_path_type.TST,      ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3'), True),
        ('TSTcT',   hc_cc_rs_path_type.TSTcT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1'), True),
        ('TcTST',   hc_cc_rs_path_type.TcTST,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'qi4', 'ci1'), True),
        ('TcTSTcT', hc_cc_rs_path_type.TcTSTcT,  ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'qi4', 'ci1', 'ci2'), True),
        ('TTcTT',   hc_cc_rs_path_type.TTcTT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1', 'ci2'), False),
        ('TcTTcT',  hc_cc_rs_path_type.TcTTcT,   ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1', 'ci2'), False),
        ('TTT',     hc_cc_rs_path_type.TTT,      ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1'), False),
        ('TcST',    hc_cc_rs_path_type.TcST,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3'), True),
        ('TScT',    hc_cc_rs_path_type.TScT,     ('length', 'cstart', 'cend', 'qi1', 'qi2'), True),
        ('TcScT',   hc_cc_rs_path_type.TcScT,    ('length', 'cstart', 'cend', 'qi1', 'qi2'), True),
    ]

    def __init__(self, kappa, sigma, discretization=0.1):
        import sys
        super().__init__(kappa, sigma, discretization)
        self._helper = _HC0pm_Reeds_Shepp(self)
        self.rs_circle_param_ = HC_CC_Circle_Param()
        self.rs_circle_param_.set_param(
            self.kappa_, sys.float_info.max,
            1.0 / self.kappa_, 0.0, 0.0, 1.0, 0.0)
        self.radius_ = self.hc_cc_circle_param_.radius
        self.mu_ = self.hc_cc_circle_param_.mu

    def hc0pm_circles_rs_path(self, c1, c2):
        """Evaluate all families and return the shortest HC_CC_RS_Path."""
        tp = hc_cc_rs_path_type
        h = self._helper
        h.distance = center_distance(c1, c2)
        h.angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        p = self.hc_cc_circle_param_

        N = nb_hc_cc_rs_paths
        length = [_INF] * N
        qi1 = [None] * N
        qi2 = [None] * N
        qi3 = [None] * N
        qi4 = [None] * N
        cstart = [None] * N
        cend = [None] * N
        ci1 = [None] * N
        ci2 = [None] * N

        skip = False

        # case E
        if configuration_equal(c1.start, c2.start):
            length[tp.E] = 0
            skip = True

        # case T (distance < epsilon)
        if not skip and h.distance < get_epsilon():
            cstart[tp.T] = HC_CC_Circle(
                c1.start, c1.left, c1.forward, _HC_REGULAR, p)
            length[tp.T] = cstart[tp.T].hc_turn_length(c2.start)
            skip = True

        if not skip:
            length, qi1, qi2, qi3, qi4, cstart, cend, ci1, ci2 = \
                _evaluate_rs_families(h, c1, c2, self._FAMILY_REGISTRY)

        # select shortest
        best = min(range(N), key=lambda i: length[i])
        return HC_CC_RS_Path(
            c1.start, c2.start, hc_cc_rs_path_type(best),
            self.kappa_, self.sigma_,
            qi1[best], qi2[best], qi3[best], qi4[best],
            cstart[best], cend[best], ci1[best], ci2[best],
            length[best])

    def _find_path(self, state1, state2):
        """Compute shortest path between two states (4×4 circle combos)."""
        start = Configuration(state1.x, state1.y, state1.theta, 0.0)
        end1 = Configuration(state2.x, state2.y, state2.theta, self.kappa_)
        end2 = Configuration(state2.x, state2.y, state2.theta, -self.kappa_)
        p = self.hc_cc_circle_param_
        rsp = self.rs_circle_param_
        start_circles = [
            HC_CC_Circle(start, True, True, True, p),
            HC_CC_Circle(start, False, True, True, p),
            HC_CC_Circle(start, True, False, True, p),
            HC_CC_Circle(start, False, False, True, p),
        ]
        end_circles = [
            HC_CC_Circle(end1, True, True, True, rsp),
            HC_CC_Circle(end2, False, True, True, rsp),
            HC_CC_Circle(end1, True, False, True, rsp),
            HC_CC_Circle(end2, False, False, True, rsp),
        ]
        best_path = None
        for sc in start_circles:
            for j, ec in enumerate(end_circles):
                if j == 0 and state2.kappa < 0:
                    continue
                if j == 1 and state2.kappa > 0:
                    continue
                if j == 2 and state2.kappa < 0:
                    continue
                if j == 3 and state2.kappa > 0:
                    continue
                path = self.hc0pm_circles_rs_path(sc, ec)
                if path is not None and (best_path is None or path.length < best_path.length):
                    best_path = path
        return best_path


# ===================================================================
# HCpm0 – ±max start, zero end curvature
# ===================================================================
class _HCpm0_Reeds_Shepp(_HC00_Reeds_Shepp):
    """Inner helper for HCpm0 path families (±max start, zero end curvature).

    In HCpm0, c1 is an RS circle (kappa=±kappa_max) and c2 is an HC circle
    (kappa=0).  Geometric properties (radius, mu, sin_mu, cos_mu) come from
    c2, not c1.  All ``_exists``, tangent-circle, and ``_path`` methods that
    reference those properties are overridden accordingly.
    """

    # ---- TT -------------------------------------------------------------
    def TT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return abs(self.distance - 2 * c2.radius) < get_epsilon()

    def TT_tangent_circles(self, c1, c2):
        x = (c1.xc + c2.xc) / 2
        y = (c1.yc + c2.yc) / 2
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        if c1.left:
            theta = angle + HALF_PI - c2.mu if c1.forward else angle + HALF_PI + c2.mu
        else:
            theta = angle - HALF_PI + c2.mu if c1.forward else angle - HALF_PI - c2.mu
        return Configuration(x, y, theta, 0)

    def TT_path(self, c1, c2):
        q2 = self.TT_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(q2, not c2.left, c2.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        length = cstart.hc_turn_length(q1) + cend.cc_turn_length(q2)
        return length, cstart, cend, q1, q2

    # ---- TcT ------------------------------------------------------------
    def TcT_path(self, c1, c2):
        q = self.TcT_tangent_circles(c1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._p())
        length = cstart.rs_turn_length(q) + cend.hc_turn_length(q)
        return length, cstart, cend, q

    # ---- TcTcT ----------------------------------------------------------
    def TcTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r = 2 * abs(c1.kappa_inv)
        delta_x = 0.5 * self.distance
        delta_y = math.sqrt(max(0.0, r ** 2 - delta_x ** 2))
        rsp = self._rsp()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, True, rsp)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, True, rsp)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2 = self.TcT_tangent_circles(tgt1, c2)
        q3 = self.TcT_tangent_circles(c1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TcTcT_path(self, c1, c2):
        qa, qb, qc, qd = self.TcTcT_tangent_circles(c1, c2)
        rsp = self._rsp()
        middle1 = HC_CC_Circle(qa, not c2.left, c2.forward, True, rsp)
        middle2 = HC_CC_Circle(qc, not c2.left, c2.forward, True, rsp)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, rsp)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._p())
        l1 = (cstart.rs_turn_length(qa) + middle1.rs_turn_length(qb) +
              cend.hc_turn_length(qb))
        l2 = (cstart.rs_turn_length(qc) + middle2.rs_turn_length(qd) +
              cend.hc_turn_length(qd))
        if l1 < l2:
            return l1, cstart, cend, qa, qb, middle1
        return l2, cstart, cend, qc, qd, middle2

    # ---- TcTT -----------------------------------------------------------
    def TcTT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return ((self.distance <= 2 * c2.radius + 2 * abs(c2.kappa_inv)) and
                (self.distance >= 2 * c2.radius - 2 * abs(c2.kappa_inv)))

    def TcTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * abs(c2.kappa_inv)
        r2 = 2 * c2.radius
        delta_x = (r1 ** 2 + self.distance ** 2 - r2 ** 2) / (2 * self.distance)
        delta_y = math.sqrt(max(0.0, r1 ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2 = self.TT_tangent_circles(tgt1, c2)
        q3 = self.TcT_tangent_circles(c1, tgt2)
        q4 = self.TT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TcTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TcTT_tangent_circles(c1, c2)
        p = self._p()
        rsp = self._rsp()
        middle1 = HC_CC_Circle(qb, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qd, not c1.left, c1.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, rsp)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        l1 = (cstart.rs_turn_length(qa) + middle1.hc_turn_length(qa) +
              cend.cc_turn_length(qb))
        l2 = (cstart.rs_turn_length(qc) + middle2.hc_turn_length(qc) +
              cend.cc_turn_length(qd))
        if l1 < l2:
            return l1, cstart, cend, qa, qb, middle1
        return l2, cstart, cend, qc, qd, middle2

    # ---- TTcT -----------------------------------------------------------
    def TTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return ((self.distance <= 2 * c2.radius + 2 * abs(c2.kappa_inv)) and
                (self.distance >= 2 * c2.radius - 2 * abs(c2.kappa_inv)))

    def TTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * c2.radius
        r2 = 2 * abs(c2.kappa_inv)
        delta_x = (r1 ** 2 + self.distance ** 2 - r2 ** 2) / (2 * self.distance)
        delta_y = math.sqrt(max(0.0, r1 ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        q1 = self.TT_tangent_circles(c1, tgt1)
        q2 = self.TcT_tangent_circles(tgt1, c2)
        q3 = self.TT_tangent_circles(c1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TTcT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTcT_tangent_circles(c1, c2)
        p = self._p()
        start1 = HC_CC_Circle(
            qa, c1.left, not c1.forward, _HC_REGULAR, p)
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True, p)
        start2 = HC_CC_Circle(
            qc, c1.left, not c1.forward, _HC_REGULAR, p)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, True, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._p())
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        l1 = (start1.hc_turn_length(q1) + middle1.hc_turn_length(qb) +
              cend.hc_turn_length(qb))
        l2 = (start2.hc_turn_length(q1) + middle2.hc_turn_length(qd) +
              cend.hc_turn_length(qd))
        if l1 < l2:
            return l1, start1, cend, q1, qb, middle1
        return l2, start2, cend, q1, qd, middle2

    # ---- TST (TiST / TeST) ---------------------------------------------
    def TiST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 2 * c2.radius

    def TiST_tangent_circles(self, c1, c2):
        dist = center_distance(c1, c2)
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        alpha = math.asin(2 * c2.radius * c2.cos_mu / dist)
        dx = c2.radius * c2.sin_mu
        dy = c2.radius * c2.cos_mu
        if c1.left and c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta, 0)
        elif c1.left and not c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif not c1.left and c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        else:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TeST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 2 * c2.radius * c2.sin_mu

    def TeST_tangent_circles(self, c1, c2):
        dx = c2.radius * c2.sin_mu
        dy = c2.radius * c2.cos_mu
        theta = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        elif c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta, 0)
        else:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TiST_path(self, c1, c2):
        q2, q3 = self.TiST_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(
            c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  cend.cc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    def TeST_path(self, c1, c2):
        q2, q3 = self.TeST_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(
            c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  cend.cc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    # ---- TSTcT (TiSTcT / TeSTcT) ---------------------------------------
    def TiSTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (2 * c2.radius * c2.sin_mu + 2 * abs(c2.kappa_inv)) ** 2 +
            (2 * c2.radius * c2.cos_mu) ** 2)

    def TeSTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= 2 * (abs(c2.kappa_inv) + c2.radius * c2.sin_mu)

    def TiSTcT_path(self, c1, c2):
        theta = self.angle
        delta_y = (4 * c2.radius * c2.cos_mu) / (abs(c2.kappa) * self.distance)
        delta_x = math.sqrt(max(0.0, (2 * c2.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q2, q3 = self.TiST_tangent_circles(c1, tgt1)
        q4 = self.TcT_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        ci = HC_CC_Circle(q3, not c1.left, c1.forward, True, p)
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  ci.hc_turn_length(q4) + cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci

    def TeSTcT_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c2.kappa_inv)
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q2, q3 = self.TeST_tangent_circles(c1, tgt1)
        q4 = self.TcT_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        ci = HC_CC_Circle(q3, not c2.left, c2.forward, True, p)
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  ci.hc_turn_length(q4) + cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci

    # ---- TcTST (TcTiST / TcTeST) ---------------------------------------
    def TcTiST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (2 * c2.radius * c2.sin_mu + 2 * abs(c2.kappa_inv)) ** 2 +
            (2 * c2.radius * c2.cos_mu) ** 2)

    def TcTeST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= 2 * (abs(c2.kappa_inv) + c2.radius * c2.sin_mu)

    def TcTiST_path(self, c1, c2):
        theta = self.angle
        delta_y = (4 * c2.radius * c2.cos_mu) / (abs(c2.kappa) * self.distance)
        delta_x = math.sqrt(max(0.0, (2 * c2.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt1 = HC_CC_Circle(
            x, y, not c1.left, not c1.forward, c1.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TiST_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        length = (cstart.rs_turn_length(q1) + ci.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + cend.cc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3, ci

    def TcTeST_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c2.kappa_inv)
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(
            x, y, not c1.left, not c1.forward, c1.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TeST_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        length = (cstart.rs_turn_length(q1) + ci.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + cend.cc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3, ci

    # ---- TcTSTcT (TcTiSTcT / TcTeSTcT) ---------------------------------
    def TcTiSTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= math.sqrt(
            (2 * c2.radius) ** 2 +
            16 * c2.radius * c2.sin_mu * abs(c2.kappa_inv) +
            (4 * c2.kappa_inv) ** 2)

    def TcTeSTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 4 * abs(c2.kappa_inv) + 2 * c2.radius * c2.sin_mu

    def TcTiSTcT_path(self, c1, c2):
        theta = self.angle
        delta_y = (4 * c2.radius * c2.cos_mu) / (self.distance * abs(c2.kappa))
        delta_x = math.sqrt(max(0.0, (2 * c2.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(
            x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TiST_tangent_circles(tgt1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        ci1 = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        ci2 = HC_CC_Circle(q3, not c2.left, c2.forward, True, p)
        length = (cstart.rs_turn_length(q1) + ci1.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + ci2.hc_turn_length(q4) +
                  cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci1, ci2

    def TcTeSTcT_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c1.kappa_inv)
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(
            x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TeST_tangent_circles(tgt1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        ci1 = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        ci2 = HC_CC_Circle(q3, not c2.left, c2.forward, True, p)
        length = (cstart.rs_turn_length(q1) + ci1.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + ci2.hc_turn_length(q4) +
                  cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci1, ci2

    # ---- TTcTT ----------------------------------------------------------
    def TTcTT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance <= 4 * c2.radius + 2 * abs(c2.kappa_inv)

    def TTcTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * abs(c2.kappa_inv)
        r2 = 2 * c2.radius
        p = self._p()
        if self.distance < 4 * c2.radius - 2 * abs(c2.kappa_inv):
            delta_x = (self.distance + r1) / 2
        else:
            delta_x = (self.distance - r1) / 2
        delta_y = math.sqrt(max(0.0, r2 ** 2 - delta_x ** 2))
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt2 = HC_CC_Circle(
            x, y, not c2.left, not c2.forward, c2.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt3 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt4 = HC_CC_Circle(
            x, y, not c2.left, not c2.forward, c2.regular, p)
        q1 = self.TT_tangent_circles(c1, tgt1)
        q2 = self.TcT_tangent_circles(tgt1, tgt2)
        q3 = self.TT_tangent_circles(tgt2, c2)
        q4 = self.TT_tangent_circles(c1, tgt3)
        q5 = self.TcT_tangent_circles(tgt3, tgt4)
        q6 = self.TT_tangent_circles(tgt4, c2)
        return q1, q2, q3, q4, q5, q6

    def TTcTT_path(self, c1, c2):
        qa, qb, qc, qd, qe, qf = self.TTcTT_tangent_circles(c1, c2)
        p = self._p()
        start1 = HC_CC_Circle(
            qa, c1.left, not c1.forward, _HC_REGULAR, p)
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qc, not c2.left, c2.forward, True, p)
        start2 = HC_CC_Circle(
            qd, c1.left, not c1.forward, _HC_REGULAR, p)
        middle3 = HC_CC_Circle(qd, not c1.left, c1.forward, True, p)
        middle4 = HC_CC_Circle(qf, not c2.left, c2.forward, True, p)
        cend = HC_CC_Circle(
            c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        l1 = (start1.hc_turn_length(q1) + middle1.hc_turn_length(qb) +
              middle2.hc_turn_length(qb) + cend.cc_turn_length(qc))
        l2 = (start2.hc_turn_length(q1) + middle3.hc_turn_length(qe) +
              middle4.hc_turn_length(qe) + cend.cc_turn_length(qf))
        if l1 < l2:
            return l1, start1, cend, q1, qb, qc, middle1, middle2
        return l2, start2, cend, q1, qe, qf, middle3, middle4

    # ---- TcTTcT ---------------------------------------------------------
    def TcTTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return ((self.distance <= 4 * abs(c2.kappa_inv) + 2 * c2.radius) and
                (self.distance >= 4 * abs(c2.kappa_inv) - 2 * c2.radius))

    def TcTTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * abs(c2.kappa_inv)
        r2 = c2.radius
        delta_x = (r1 ** 2 + (self.distance / 2) ** 2 - r2 ** 2) / self.distance
        delta_y = math.sqrt(max(0.0, r1 ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(
            x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt3 = HC_CC_Circle(
            x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt4 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2 = self.TT_tangent_circles(tgt1, tgt2)
        q3 = self.TcT_tangent_circles(tgt2, c2)
        q4 = self.TcT_tangent_circles(c1, tgt3)
        q5 = self.TT_tangent_circles(tgt3, tgt4)
        q6 = self.TcT_tangent_circles(tgt4, c2)
        return q1, q2, q3, q4, q5, q6

    def TcTTcT_path(self, c1, c2):
        qa, qb, qc, qd, qe, qf = self.TcTTcT_tangent_circles(c1, c2)
        p = self._p()
        rsp = self._rsp()
        middle1 = HC_CC_Circle(qb, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qb, c1.left, not c1.forward, True, p)
        middle3 = HC_CC_Circle(qe, not c1.left, c1.forward, True, p)
        middle4 = HC_CC_Circle(qe, c1.left, not c1.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, rsp)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        l1 = (cstart.rs_turn_length(qa) + middle1.hc_turn_length(qa) +
              middle2.hc_turn_length(qc) + cend.hc_turn_length(qc))
        l2 = (cstart.rs_turn_length(qd) + middle3.hc_turn_length(qd) +
              middle4.hc_turn_length(qf) + cend.hc_turn_length(qf))
        if l1 < l2:
            return l1, cstart, cend, qa, qc, middle1, middle2
        return l2, cstart, cend, qd, qf, middle3, middle4

    # ---- TTT ------------------------------------------------------------
    def TTT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance <= 4 * c2.radius

    def TTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r = 2 * c2.radius
        delta_x = 0.5 * self.distance
        delta_y = math.sqrt(max(0.0, r ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        q1 = self.TT_tangent_circles(c1, tgt1)
        q2 = self.TT_tangent_circles(tgt1, c2)
        q3 = self.TT_tangent_circles(c1, tgt2)
        q4 = self.TT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTT_tangent_circles(c1, c2)
        p = self._p()
        start1 = HC_CC_Circle(
            qa, c1.left, not c1.forward, _HC_REGULAR, p)
        middle1 = HC_CC_Circle(
            qa, not c1.left, c1.forward, _CC_REGULAR, p)
        start2 = HC_CC_Circle(
            qc, c1.left, not c1.forward, _HC_REGULAR, p)
        middle2 = HC_CC_Circle(
            qc, not c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(
            c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        l1 = (start1.hc_turn_length(q1) + middle1.cc_turn_length(qb) +
              cend.cc_turn_length(qb))
        l2 = (start2.hc_turn_length(q1) + middle2.cc_turn_length(qd) +
              cend.cc_turn_length(qd))
        if l1 < l2:
            return l1, start1, cend, q1, qb, middle1
        return l2, start2, cend, q1, qd, middle2

    # ---- TcST (TciST / TceST) ------------------------------------------
    def TciST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (c2.radius * c2.sin_mu) ** 2 +
            (c2.radius * c2.cos_mu + abs(c2.kappa_inv)) ** 2)

    def TceST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (c2.radius * c2.sin_mu) ** 2 +
            (c2.radius * c2.cos_mu - abs(c2.kappa_inv)) ** 2)

    def TciST_path(self, c1, c2):
        alpha = math.asin(
            (c2.radius * c2.cos_mu + abs(c2.kappa_inv)) / self.distance)
        dx1, dy1 = 0.0, abs(c2.kappa_inv)
        dx2 = c2.radius * c2.sin_mu
        dy2 = c2.radius * c2.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(
            c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.rs_turn_length(q1) +
                  configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TceST_path(self, c1, c2):
        alpha = math.asin(
            (c2.radius * c2.cos_mu - abs(c2.kappa_inv)) / self.distance)
        dx1, dy1 = 0.0, abs(c2.kappa_inv)
        dx2 = c2.radius * c2.sin_mu
        dy2 = c2.radius * c2.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(
            c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.rs_turn_length(q1) +
                  configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
        return length, cstart, cend, q1, q2

    # ---- TScT (TiScT / TeScT) ------------------------------------------
    def TiScT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (c2.radius * c2.sin_mu) ** 2 +
            (c2.radius * c2.cos_mu + abs(c2.kappa_inv)) ** 2)

    def TeScT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= math.sqrt(
            (c2.radius * c2.sin_mu) ** 2 +
            (c2.radius * c2.cos_mu - abs(c2.kappa_inv)) ** 2)

    def TiScT_path(self, c1, c2):
        alpha = math.asin(
            (c2.radius * c2.cos_mu + abs(c2.kappa_inv)) / self.distance)
        dx1 = c2.radius * c2.sin_mu
        dy1 = c2.radius * c2.cos_mu
        dx2, dy2 = 0.0, abs(c2.kappa_inv)
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q2 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q3 = Configuration(x, y, theta, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q2 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q3 = Configuration(x, y, theta + PI, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q2 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q3 = Configuration(x, y, theta, c2.kappa)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q2 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q3 = Configuration(x, y, theta + PI, c2.kappa)
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  cend.hc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    def TeScT_path(self, c1, c2):
        alpha = math.asin(
            (c2.radius * c2.cos_mu - abs(c2.kappa_inv)) / self.distance)
        dx1 = c2.radius * c2.sin_mu
        dy1 = c2.radius * c2.cos_mu
        dx2, dy2 = 0.0, abs(c2.kappa_inv)
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q2 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q3 = Configuration(x, y, theta, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q2 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q3 = Configuration(x, y, theta + PI, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q2 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q3 = Configuration(x, y, theta, c2.kappa)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q2 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q3 = Configuration(x, y, theta + PI, c2.kappa)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  cend.hc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    # ---- TcScT (TciScT / TceScT) ---------------------------------------
    def TciScT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance > 2 * abs(c2.kappa_inv)

    def TciScT_path(self, c1, c2):
        alpha = math.asin(2 / (abs(c2.kappa) * self.distance))
        dx, dy = 0.0, abs(c2.kappa_inv)
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        else:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        length = (cstart.rs_turn_length(q1) +
                  configuration_distance(q1, q2) +
                  cend.hc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TceScT_path(self, c1, c2):
        theta = self.angle
        dx, dy = 0.0, abs(c2.kappa_inv)
        p = self._p()
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        else:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, p)
        length = (cstart.rs_turn_length(q1) +
                  configuration_distance(q1, q2) +
                  cend.hc_turn_length(q2))
        return length, cstart, cend, q1, q2


class HCpm0_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """HC Reeds-Shepp with ±max curvature at start, zero at end."""

    _CONTROL_TABLE = {
        hc_cc_rs_path_type.E: [('empty',)],
        hc_cc_rs_path_type.T: [('hc', 'cend', 'start', False)],
        hc_cc_rs_path_type.TT: [('hc', 'cstart', 'qi1', False), ('cc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TcT: [('rs', 'cstart', 'qi1', True), ('hc', 'cend', 'qi1', False)],
        hc_cc_rs_path_type.TcTcT: [('rs', 'cstart', 'qi1', True), ('rs', 'ci1', 'qi2', True), ('hc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TcTT: [('rs', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('cc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TTcT: [('hc', 'cstart', 'qi1', False), ('hc', 'ci1', 'qi2', True), ('hc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TST: [('hc', 'cstart', 'qi1', False), ('straight', 'qi2', 'qi3'), ('cc', 'cend', 'qi3', False)],
        hc_cc_rs_path_type.TSTcT: [('hc', 'cstart', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'ci1', 'qi4', True), ('hc', 'cend', 'qi4', False)],
        hc_cc_rs_path_type.TcTST: [('rs', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('straight', 'qi2', 'qi3'), ('cc', 'cend', 'qi3', False)],
        hc_cc_rs_path_type.TcTSTcT: [('rs', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'ci2', 'qi4', True), ('hc', 'cend', 'qi4', False)],
        hc_cc_rs_path_type.TTcTT: [('hc', 'cstart', 'qi1', False), ('hc', 'ci1', 'qi2', True), ('hc', 'ci2', 'qi2', False), ('cc', 'cend', 'qi3', False)],
        hc_cc_rs_path_type.TcTTcT: [('rs', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('hc', 'ci2', 'qi2', True), ('hc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TTT: [('hc', 'cstart', 'qi1', False), ('cc', 'ci1', 'qi2', True), ('cc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TcST: [('rs', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('cc', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TScT: [('hc', 'cstart', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'cend', 'qi3', False)],
        hc_cc_rs_path_type.TcScT: [('rs', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('hc', 'cend', 'qi2', False)],
    }

    _FAMILY_REGISTRY = [
        ('TT',      hc_cc_rs_path_type.TT,      ('length', 'cstart', 'cend', 'qi1', 'qi2'), False),
        ('TcT',     hc_cc_rs_path_type.TcT,      ('length', 'cstart', 'cend', 'qi1'), False),
        ('TcTcT',   hc_cc_rs_path_type.TcTcT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TcTT',    hc_cc_rs_path_type.TcTT,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TTcT',    hc_cc_rs_path_type.TTcT,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TST',     hc_cc_rs_path_type.TST,      ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3'), True),
        ('TSTcT',   hc_cc_rs_path_type.TSTcT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'qi4', 'ci1'), True),
        ('TcTST',   hc_cc_rs_path_type.TcTST,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1'), True),
        ('TcTSTcT', hc_cc_rs_path_type.TcTSTcT,  ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'qi4', 'ci1', 'ci2'), True),
        ('TTcTT',   hc_cc_rs_path_type.TTcTT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1', 'ci2'), False),
        ('TcTTcT',  hc_cc_rs_path_type.TcTTcT,   ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1', 'ci2'), False),
        ('TTT',     hc_cc_rs_path_type.TTT,      ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TcST',    hc_cc_rs_path_type.TcST,     ('length', 'cstart', 'cend', 'qi1', 'qi2'), True),
        ('TScT',    hc_cc_rs_path_type.TScT,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3'), True),
        ('TcScT',   hc_cc_rs_path_type.TcScT,    ('length', 'cstart', 'cend', 'qi1', 'qi2'), True),
    ]

    def __init__(self, kappa, sigma, discretization=0.1):
        import sys
        super().__init__(kappa, sigma, discretization)
        self._helper = _HCpm0_Reeds_Shepp(self)
        self.rs_circle_param_ = HC_CC_Circle_Param()
        self.rs_circle_param_.set_param(
            self.kappa_, sys.float_info.max,
            1.0 / self.kappa_, 0.0, 0.0, 1.0, 0.0)
        self.radius_ = self.hc_cc_circle_param_.radius
        self.mu_ = self.hc_cc_circle_param_.mu

    def hcpm0_circles_rs_path(self, c1, c2):
        """Evaluate all families and return the shortest HC_CC_RS_Path."""
        tp = hc_cc_rs_path_type
        h = self._helper
        h.distance = center_distance(c1, c2)
        h.angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        p = self.hc_cc_circle_param_

        N = nb_hc_cc_rs_paths
        length = [_INF] * N
        qi1 = [None] * N
        qi2 = [None] * N
        qi3 = [None] * N
        qi4 = [None] * N
        cstart = [None] * N
        cend = [None] * N
        ci1 = [None] * N
        ci2 = [None] * N

        skip = False

        # case E
        if configuration_equal(c1.start, c2.start):
            length[tp.E] = 0
            skip = True

        # case T (distance < epsilon)
        if not skip and h.distance < get_epsilon():
            cend[tp.T] = HC_CC_Circle(
                c2.start, c2.left, c2.forward, _HC_REGULAR, p)
            length[tp.T] = cend[tp.T].hc_turn_length(c1.start)
            skip = True

        if not skip:
            length, qi1, qi2, qi3, qi4, cstart, cend, ci1, ci2 = \
                _evaluate_rs_families(h, c1, c2, self._FAMILY_REGISTRY)

        # select shortest
        best = min(range(N), key=lambda i: length[i])
        return HC_CC_RS_Path(
            c1.start, c2.start, hc_cc_rs_path_type(best),
            self.kappa_, self.sigma_,
            qi1[best], qi2[best], qi3[best], qi4[best],
            cstart[best], cend[best], ci1[best], ci2[best],
            length[best])

    def _find_path(self, state1, state2):
        """Compute shortest path between two states (4×4 circle combos)."""
        start1 = Configuration(state1.x, state1.y, state1.theta, self.kappa_)
        start2 = Configuration(state1.x, state1.y, state1.theta, -self.kappa_)
        end = Configuration(state2.x, state2.y, state2.theta, 0.0)
        rsp = self.rs_circle_param_
        p = self.hc_cc_circle_param_
        start_circles = [
            HC_CC_Circle(start1, True, True, True, rsp),
            HC_CC_Circle(start2, False, True, True, rsp),
            HC_CC_Circle(start1, True, False, True, rsp),
            HC_CC_Circle(start2, False, False, True, rsp),
        ]
        end_circles = [
            HC_CC_Circle(end, True, True, True, p),
            HC_CC_Circle(end, False, True, True, p),
            HC_CC_Circle(end, True, False, True, p),
            HC_CC_Circle(end, False, False, True, p),
        ]
        best_path = None
        for i, sc in enumerate(start_circles):
            if i == 0 and state1.kappa < 0:
                continue
            if i == 1 and state1.kappa > 0:
                continue
            if i == 2 and state1.kappa < 0:
                continue
            if i == 3 and state1.kappa > 0:
                continue
            for ec in end_circles:
                path = self.hcpm0_circles_rs_path(sc, ec)
                if path is not None and (best_path is None
                                         or path.length < best_path.length):
                    best_path = path
        return best_path


# ===================================================================
# HCpmpm – ±max start, ±max end curvature
# ===================================================================
class _HCpmpm_Reeds_Shepp(_HC00_Reeds_Shepp):
    """Inner helper for HCpmpm path families (±max start, ±max end curvature).

    In HCpmpm, both c1 and c2 are RS circles (kappa=±kappa_max).  Geometric
    properties (radius, mu, sin_mu, cos_mu) come from the parent's HC circle
    param, not from c1/c2.
    """

    # ---- TT -------------------------------------------------------------
    def TT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return abs(self.distance - 2 * self._parent.radius_) < get_epsilon()

    def TT_tangent_circles(self, c1, c2):
        x = (c1.xc + c2.xc) / 2
        y = (c1.yc + c2.yc) / 2
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        mu = self._parent.mu_
        if c1.left:
            theta = angle + HALF_PI - mu if c1.forward else angle + HALF_PI + mu
        else:
            theta = angle - HALF_PI + mu if c1.forward else angle - HALF_PI - mu
        return Configuration(x, y, theta, 0)

    def TT_path(self, c1, c2):
        q2 = self.TT_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(q2, c2.left, not c2.forward, _HC_REGULAR, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = cstart.hc_turn_length(q1) + cend.hc_turn_length(q3)
        return length, cstart, cend, q1, q2, q3

    # ---- TcT ------------------------------------------------------------
    def TcT_path(self, c1, c2):
        q = self.TcT_tangent_circles(c1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._rsp())
        length = cstart.rs_turn_length(q) + cend.rs_turn_length(q)
        return length, cstart, cend, q

    # ---- TcTcT ----------------------------------------------------------
    def TcTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r = 2 * abs(c1.kappa_inv)
        delta_x = 0.5 * self.distance
        delta_y = math.sqrt(max(0.0, r ** 2 - delta_x ** 2))
        rsp = self._rsp()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, True, rsp)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, not c1.forward, True, rsp)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2 = self.TcT_tangent_circles(tgt1, c2)
        q3 = self.TcT_tangent_circles(c1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TcTcT_path(self, c1, c2):
        qa, qb, qc, qd = self.TcTcT_tangent_circles(c1, c2)
        rsp = self._rsp()
        middle1 = HC_CC_Circle(qa, not c1.left, not c1.forward, True, rsp)
        middle2 = HC_CC_Circle(qc, not c1.left, not c1.forward, True, rsp)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, rsp)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        l1 = (cstart.rs_turn_length(qa) + middle1.rs_turn_length(qb) +
              cend.rs_turn_length(qb))
        l2 = (cstart.rs_turn_length(qc) + middle2.rs_turn_length(qd) +
              cend.rs_turn_length(qd))
        if l1 < l2:
            return l1, cstart, cend, qa, qb, middle1
        return l2, cstart, cend, qc, qd, middle2

    # ---- TcTT -----------------------------------------------------------
    def TcTT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return ((self.distance <= 2 * self._parent.radius_ + 2 * abs(c1.kappa_inv)) and
                (self.distance >= 2 * self._parent.radius_ - 2 * abs(c1.kappa_inv)))

    def TcTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * abs(c1.kappa_inv)
        r2 = 2 * self._parent.radius_
        delta_x = (r1 ** 2 + self.distance ** 2 - r2 ** 2) / (2 * self.distance)
        delta_y = math.sqrt(max(0.0, r1 ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2 = self.TT_tangent_circles(tgt1, c2)
        q3 = self.TcT_tangent_circles(c1, tgt2)
        q4 = self.TT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TcTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TcTT_tangent_circles(c1, c2)
        p = self._p()
        end1 = HC_CC_Circle(
            qb, c2.left, not c2.forward, _HC_REGULAR, p)
        end2 = HC_CC_Circle(
            qd, c2.left, not c2.forward, _HC_REGULAR, p)
        middle1 = HC_CC_Circle(qb, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qd, not c1.left, c1.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        q2 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        l1 = (cstart.rs_turn_length(qa) + middle1.hc_turn_length(qa) +
              end1.hc_turn_length(q2))
        l2 = (cstart.rs_turn_length(qc) + middle2.hc_turn_length(qc) +
              end2.hc_turn_length(q2))
        if l1 < l2:
            return l1, cstart, end1, qa, q2, middle1
        return l2, cstart, end2, qc, q2, middle2

    # ---- TTcT -----------------------------------------------------------
    def TTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return ((self.distance <= 2 * self._parent.radius_ + 2 * abs(c1.kappa_inv)) and
                (self.distance >= 2 * self._parent.radius_ - 2 * abs(c1.kappa_inv)))

    def TTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * self._parent.radius_
        r2 = 2 * abs(c1.kappa_inv)
        delta_x = (r1 ** 2 + self.distance ** 2 - r2 ** 2) / (2 * self.distance)
        delta_y = math.sqrt(max(0.0, r1 ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        q1 = self.TT_tangent_circles(c1, tgt1)
        q2 = self.TcT_tangent_circles(tgt1, c2)
        q3 = self.TT_tangent_circles(c1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TTcT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTcT_tangent_circles(c1, c2)
        p = self._p()
        start1 = HC_CC_Circle(
            qa, c1.left, not c1.forward, _HC_REGULAR, p)
        start2 = HC_CC_Circle(
            qc, c1.left, not c1.forward, _HC_REGULAR, p)
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, True, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._rsp())
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        l1 = (start1.hc_turn_length(q1) + middle1.hc_turn_length(qb) +
              cend.rs_turn_length(qb))
        l2 = (start2.hc_turn_length(q1) + middle2.hc_turn_length(qd) +
              cend.rs_turn_length(qd))
        if l1 < l2:
            return l1, start1, cend, q1, qb, middle1
        return l2, start2, cend, q1, qd, middle2

    # ---- TST (TiST / TeST) ---------------------------------------------
    def TiST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 2 * self._parent.radius_

    def TiST_tangent_circles(self, c1, c2):
        dist = center_distance(c1, c2)
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        cos_mu = self._parent.cos_mu_
        sin_mu = self._parent.sin_mu_
        radius = self._parent.radius_
        alpha = math.asin(2 * radius * cos_mu / dist)
        dx = radius * sin_mu
        dy = radius * cos_mu
        if c1.left and c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta, 0)
        elif c1.left and not c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif not c1.left and c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        else:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TeST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        sin_mu = self._parent.sin_mu_
        return self.distance >= 2 * self._parent.radius_ * sin_mu

    def TeST_tangent_circles(self, c1, c2):
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        dx = radius * sin_mu
        dy = radius * cos_mu
        theta = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        elif c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta, 0)
        else:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TiST_path(self, c1, c2):
        q2, q3 = self.TiST_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(
            q3, c2.left, not c2.forward, _HC_REGULAR, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        q4 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4

    def TeST_path(self, c1, c2):
        q2, q3 = self.TeST_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(
            q3, c2.left, not c2.forward, _HC_REGULAR, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        q4 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4

    # ---- TSTcT (TiSTcT / TeSTcT) ---------------------------------------
    def TiSTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        return self.distance >= math.sqrt(
            (2 * radius * sin_mu + 2 * abs(c1.kappa_inv)) ** 2 +
            (2 * radius * cos_mu) ** 2)

    def TeSTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        return self.distance >= 2 * (abs(c1.kappa_inv) + radius * sin_mu)

    def TiSTcT_path(self, c1, c2):
        theta = self.angle
        radius = self._parent.radius_
        cos_mu = self._parent.cos_mu_
        delta_y = (4 * radius * cos_mu) / (abs(c2.kappa) * self.distance)
        delta_x = math.sqrt(max(0.0, (2 * c2.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q2, q3 = self.TiST_tangent_circles(c1, tgt1)
        q4 = self.TcT_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._rsp())
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        ci = HC_CC_Circle(q3, not c1.left, c1.forward, True, p)
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  ci.hc_turn_length(q4) + cend.rs_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci

    def TeSTcT_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c2.kappa_inv)
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q2, q3 = self.TeST_tangent_circles(c1, tgt1)
        q4 = self.TcT_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._rsp())
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        ci = HC_CC_Circle(q3, c1.left, c1.forward, True, p)
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  ci.hc_turn_length(q4) + cend.rs_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci

    # ---- TcTST (TcTiST / TcTeST) ---------------------------------------
    def TcTiST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        return self.distance >= math.sqrt(
            (2 * radius * sin_mu + 2 * abs(c1.kappa_inv)) ** 2 +
            (2 * radius * cos_mu) ** 2)

    def TcTeST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        return self.distance >= 2 * (abs(c1.kappa_inv) + radius * sin_mu)

    def TcTiST_path(self, c1, c2):
        theta = self.angle
        radius = self._parent.radius_
        cos_mu = self._parent.cos_mu_
        delta_y = (4 * radius * cos_mu) / (abs(c2.kappa) * self.distance)
        delta_x = math.sqrt(max(0.0, (2 * c2.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt1 = HC_CC_Circle(
            x, y, not c2.left, not c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TiST_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(
            q3, c2.left, not c2.forward, _HC_REGULAR, p)
        q4 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        length = (cstart.rs_turn_length(q1) + ci.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci

    def TcTeST_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c2.kappa_inv)
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(
            x, y, c2.left, not c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TeST_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(
            q3, c2.left, not c2.forward, _HC_REGULAR, p)
        q4 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        length = (cstart.rs_turn_length(q1) + ci.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + cend.hc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci

    # ---- TcTSTcT (TcTiSTcT / TcTeSTcT) ---------------------------------
    def TcTiSTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        return self.distance >= math.sqrt(
            (2 * radius) ** 2 +
            16 * radius * sin_mu * abs(c1.kappa_inv) +
            (4 * c1.kappa_inv) ** 2)

    def TcTeSTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        return self.distance >= 4 * abs(c1.kappa_inv) + 2 * radius * sin_mu

    def TcTiSTcT_path(self, c1, c2):
        theta = self.angle
        radius = self._parent.radius_
        cos_mu = self._parent.cos_mu_
        delta_y = (4 * radius * cos_mu) / (self.distance * abs(c1.kappa))
        delta_x = math.sqrt(max(0.0, (2 * c1.kappa_inv) ** 2 - delta_y ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(
            x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TiST_tangent_circles(tgt1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._rsp())
        ci1 = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        ci2 = HC_CC_Circle(q3, not c2.left, c2.forward, True, p)
        length = (cstart.rs_turn_length(q1) + ci1.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + ci2.hc_turn_length(q4) +
                  cend.rs_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci1, ci2

    def TcTeSTcT_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * abs(c1.kappa_inv)
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(
            x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TeST_tangent_circles(tgt1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._rsp())
        ci1 = HC_CC_Circle(q2, not c1.left, c1.forward, True, p)
        ci2 = HC_CC_Circle(q3, not c2.left, c2.forward, True, p)
        length = (cstart.rs_turn_length(q1) + ci1.hc_turn_length(q1) +
                  configuration_distance(q2, q3) + ci2.hc_turn_length(q4) +
                  cend.rs_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci1, ci2

    # ---- TTcTT ----------------------------------------------------------
    def TTcTT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance <= 4 * self._parent.radius_ + 2 * abs(c1.kappa_inv)

    def TTcTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * abs(c1.kappa_inv)
        r2 = 2 * self._parent.radius_
        p = self._p()
        if self.distance < 4 * self._parent.radius_ - 2 * abs(c1.kappa_inv):
            delta_x = (self.distance + r1) / 2
        else:
            delta_x = (self.distance - r1) / 2
        delta_y = math.sqrt(max(0.0, r2 ** 2 - delta_x ** 2))
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt2 = HC_CC_Circle(
            x, y, not c2.left, not c2.forward, c2.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt3 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt4 = HC_CC_Circle(
            x, y, not c2.left, not c2.forward, c2.regular, p)
        q1 = self.TT_tangent_circles(c1, tgt1)
        q2 = self.TcT_tangent_circles(tgt1, tgt2)
        q3 = self.TT_tangent_circles(tgt2, c2)
        q4 = self.TT_tangent_circles(c1, tgt3)
        q5 = self.TcT_tangent_circles(tgt3, tgt4)
        q6 = self.TT_tangent_circles(tgt4, c2)
        return q1, q2, q3, q4, q5, q6

    def TTcTT_path(self, c1, c2):
        qa, qb, qc, qd, qe, qf = self.TTcTT_tangent_circles(c1, c2)
        p = self._p()
        start1 = HC_CC_Circle(
            qa, c1.left, not c1.forward, _HC_REGULAR, p)
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qc, not c2.left, c2.forward, True, p)
        end1 = HC_CC_Circle(
            qc, c2.left, not c2.forward, _HC_REGULAR, p)
        start2 = HC_CC_Circle(
            qd, c1.left, not c1.forward, _HC_REGULAR, p)
        middle3 = HC_CC_Circle(qd, not c1.left, c1.forward, True, p)
        middle4 = HC_CC_Circle(qf, not c2.left, c2.forward, True, p)
        end2 = HC_CC_Circle(
            qf, c2.left, not c2.forward, _HC_REGULAR, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        l1 = (start1.hc_turn_length(q1) + middle1.hc_turn_length(qb) +
              middle2.hc_turn_length(qb) + end1.hc_turn_length(q3))
        l2 = (start2.hc_turn_length(q1) + middle3.hc_turn_length(qe) +
              middle4.hc_turn_length(qe) + end2.hc_turn_length(q3))
        if l1 < l2:
            return l1, start1, end1, q1, qb, q3, middle1, middle2
        return l2, start2, end2, q1, qe, q3, middle3, middle4

    # ---- TcTTcT ---------------------------------------------------------
    def TcTTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return ((self.distance <= 4 * abs(c1.kappa_inv) + 2 * self._parent.radius_) and
                (self.distance >= 4 * abs(c1.kappa_inv) - 2 * self._parent.radius_))

    def TcTTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * abs(c1.kappa_inv)
        r2 = self._parent.radius_
        delta_x = (r1 ** 2 + (self.distance / 2) ** 2 - r2 ** 2) / self.distance
        delta_y = math.sqrt(max(0.0, r1 ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(
            x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt3 = HC_CC_Circle(
            x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt4 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2 = self.TT_tangent_circles(tgt1, tgt2)
        q3 = self.TcT_tangent_circles(tgt2, c2)
        q4 = self.TcT_tangent_circles(c1, tgt3)
        q5 = self.TT_tangent_circles(tgt3, tgt4)
        q6 = self.TcT_tangent_circles(tgt4, c2)
        return q1, q2, q3, q4, q5, q6

    def TcTTcT_path(self, c1, c2):
        qa, qb, qc, qd, qe, qf = self.TcTTcT_tangent_circles(c1, c2)
        p = self._p()
        rsp = self._rsp()
        middle1 = HC_CC_Circle(qb, not c1.left, c1.forward, True, p)
        middle2 = HC_CC_Circle(qb, c1.left, not c1.forward, True, p)
        middle3 = HC_CC_Circle(qe, not c1.left, c1.forward, True, p)
        middle4 = HC_CC_Circle(qe, c1.left, not c1.forward, True, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular, rsp)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular, rsp)
        l1 = (cstart.rs_turn_length(qa) + middle1.hc_turn_length(qa) +
              middle2.hc_turn_length(qc) + cend.rs_turn_length(qc))
        l2 = (cstart.rs_turn_length(qd) + middle3.hc_turn_length(qd) +
              middle4.hc_turn_length(qf) + cend.rs_turn_length(qf))
        if l1 < l2:
            return l1, cstart, cend, qa, qc, middle1, middle2
        return l2, cstart, cend, qd, qf, middle3, middle4

    # ---- TTT ------------------------------------------------------------
    def TTT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance <= 4 * self._parent.radius_

    def TTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r = 2 * self._parent.radius_
        delta_x = 0.5 * self.distance
        delta_y = math.sqrt(max(0.0, r ** 2 - delta_x ** 2))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular, p)
        q1 = self.TT_tangent_circles(c1, tgt1)
        q2 = self.TT_tangent_circles(tgt1, c2)
        q3 = self.TT_tangent_circles(c1, tgt2)
        q4 = self.TT_tangent_circles(tgt2, c2)
        return q1, q2, q3, q4

    def TTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTT_tangent_circles(c1, c2)
        p = self._p()
        start1 = HC_CC_Circle(
            qa, c1.left, not c1.forward, _HC_REGULAR, p)
        middle1 = HC_CC_Circle(
            qa, not c1.left, c1.forward, _CC_REGULAR, p)
        end1 = HC_CC_Circle(
            qb, c2.left, not c2.forward, _HC_REGULAR, p)
        start2 = HC_CC_Circle(
            qc, c1.left, not c1.forward, _HC_REGULAR, p)
        middle2 = HC_CC_Circle(
            qc, not c1.left, c1.forward, _CC_REGULAR, p)
        end2 = HC_CC_Circle(
            qd, c2.left, not c2.forward, _HC_REGULAR, p)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        l1 = (start1.hc_turn_length(q1) + middle1.cc_turn_length(qb) +
              end1.hc_turn_length(q3))
        l2 = (start2.hc_turn_length(q1) + middle2.cc_turn_length(qd) +
              end2.hc_turn_length(q3))
        if l1 < l2:
            return l1, start1, end1, q1, qb, q3, middle1
        return l2, start2, end2, q1, qd, q3, middle2

    # ---- TcST (TciST / TceST) ------------------------------------------
    def TciST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        return self.distance >= math.sqrt(
            (radius * sin_mu) ** 2 +
            (radius * cos_mu + abs(c1.kappa_inv)) ** 2)

    def TceST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        return self.distance >= math.sqrt(
            (radius * sin_mu) ** 2 +
            (radius * cos_mu - abs(c1.kappa_inv)) ** 2)

    def TciST_path(self, c1, c2):
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        alpha = math.asin(
            (radius * cos_mu + abs(c1.kappa_inv)) / self.distance)
        dx1, dy1 = 0.0, abs(c1.kappa_inv)
        dx2 = radius * sin_mu
        dy2 = radius * cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta, 0)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(
            q2, c2.left, not c2.forward, _HC_REGULAR, p)
        length = (cstart.rs_turn_length(q1) +
                  configuration_distance(q1, q2) +
                  cend.hc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    def TceST_path(self, c1, c2):
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        alpha = math.asin(
            (radius * cos_mu - abs(c1.kappa_inv)) / self.distance)
        dx1, dy1 = 0.0, abs(c1.kappa_inv)
        dx2 = radius * sin_mu
        dy2 = radius * cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2)
            q2 = Configuration(x, y, theta, 0)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(
            q2, c2.left, not c2.forward, _HC_REGULAR, p)
        length = (cstart.rs_turn_length(q1) +
                  configuration_distance(q1, q2) +
                  cend.hc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    # ---- TScT (TiScT / TeScT) ------------------------------------------
    def TiScT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        return self.distance >= math.sqrt(
            (radius * sin_mu) ** 2 +
            (radius * cos_mu + abs(c1.kappa_inv)) ** 2)

    def TeScT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        return self.distance >= math.sqrt(
            (radius * sin_mu) ** 2 +
            (radius * cos_mu - abs(c1.kappa_inv)) ** 2)

    def TiScT_path(self, c1, c2):
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        alpha = math.asin(
            (radius * cos_mu + abs(c1.kappa_inv)) / self.distance)
        dx1 = radius * sin_mu
        dy1 = radius * cos_mu
        dx2, dy2 = 0.0, abs(c1.kappa_inv)
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q2 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q3 = Configuration(x, y, theta, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q2 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q3 = Configuration(x, y, theta + PI, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q2 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q3 = Configuration(x, y, theta, c2.kappa)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q2 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q3 = Configuration(x, y, theta + PI, c2.kappa)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._rsp())
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  cend.rs_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    def TeScT_path(self, c1, c2):
        radius = self._parent.radius_
        sin_mu = self._parent.sin_mu_
        cos_mu = self._parent.cos_mu_
        alpha = math.asin(
            (radius * cos_mu - abs(c1.kappa_inv)) / self.distance)
        dx1 = radius * sin_mu
        dy1 = radius * cos_mu
        dx2, dy2 = 0.0, abs(c1.kappa_inv)
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q2 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q3 = Configuration(x, y, theta, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q2 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q3 = Configuration(x, y, theta + PI, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1)
            q2 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2)
            q3 = Configuration(x, y, theta, c2.kappa)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1)
            q2 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2)
            q3 = Configuration(x, y, theta + PI, c2.kappa)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        cstart = HC_CC_Circle(
            q2, c1.left, not c1.forward, _HC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._rsp())
        length = (cstart.hc_turn_length(q1) +
                  configuration_distance(q2, q3) +
                  cend.rs_turn_length(q3))
        return length, cstart, cend, q1, q2, q3

    # ---- TcScT (TciScT / TceScT) ---------------------------------------
    def TciScT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance > 2 * abs(c1.kappa_inv)

    def TceScT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= get_epsilon()

    def TciScT_path(self, c1, c2):
        alpha = math.asin(2 / (abs(c1.kappa) * self.distance))
        dx, dy = 0.0, abs(c1.kappa_inv)
        if c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif c1.left and not c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif not c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        else:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._rsp())
        length = (cstart.rs_turn_length(q1) +
                  configuration_distance(q1, q2) +
                  cend.rs_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TceScT_path(self, c1, c2):
        theta = self.angle
        dx, dy = 0.0, abs(c1.kappa_inv)
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        elif c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        elif not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, c2.kappa)
        else:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, c1.kappa)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, c2.kappa)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                              self._rsp())
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                            self._rsp())
        length = (cstart.rs_turn_length(q1) +
                  configuration_distance(q1, q2) +
                  cend.rs_turn_length(q2))
        return length, cstart, cend, q1, q2


class HCpmpm_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """HC Reeds-Shepp with ±max curvature at start and end."""

    _CONTROL_TABLE = {
        hc_cc_rs_path_type.E: [('empty',)],
        hc_cc_rs_path_type.T: [('rs', 'cstart', 'end', True)],
        hc_cc_rs_path_type.TT: [('hc', 'cstart', 'qi1', False), ('hc', 'cend', 'qi3', True)],
        hc_cc_rs_path_type.TcT: [('rs', 'cstart', 'qi1', True), ('rs', 'cend', 'qi1', False)],
        hc_cc_rs_path_type.TcTcT: [('rs', 'cstart', 'qi1', True), ('rs', 'ci1', 'qi2', True), ('rs', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TcTT: [('rs', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('hc', 'cend', 'qi2', True)],
        hc_cc_rs_path_type.TTcT: [('hc', 'cstart', 'qi1', False), ('hc', 'ci1', 'qi2', True), ('rs', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TST: [('hc', 'cstart', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'cend', 'qi4', True)],
        hc_cc_rs_path_type.TSTcT: [('hc', 'cstart', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'ci1', 'qi4', True), ('rs', 'cend', 'qi4', False)],
        hc_cc_rs_path_type.TcTST: [('rs', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'cend', 'qi4', True)],
        hc_cc_rs_path_type.TcTSTcT: [('rs', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'ci2', 'qi4', True), ('rs', 'cend', 'qi4', False)],
        hc_cc_rs_path_type.TTcTT: [('hc', 'cstart', 'qi1', False), ('hc', 'ci1', 'qi2', True), ('hc', 'ci2', 'qi2', False), ('hc', 'cend', 'qi3', True)],
        hc_cc_rs_path_type.TcTTcT: [('rs', 'cstart', 'qi1', True), ('hc', 'ci1', 'qi1', False), ('hc', 'ci2', 'qi2', True), ('rs', 'cend', 'qi2', False)],
        hc_cc_rs_path_type.TTT: [('hc', 'cstart', 'qi1', False), ('cc', 'ci1', 'qi2', True), ('hc', 'cend', 'qi3', True)],
        hc_cc_rs_path_type.TcST: [('rs', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('hc', 'cend', 'qi3', True)],
        hc_cc_rs_path_type.TScT: [('hc', 'cstart', 'qi1', False), ('straight', 'qi2', 'qi3'), ('rs', 'cend', 'qi3', False)],
        hc_cc_rs_path_type.TcScT: [('rs', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('rs', 'cend', 'qi2', False)],
    }

    _FAMILY_REGISTRY = [
        ('TT',      hc_cc_rs_path_type.TT,      ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3'), False),
        ('TcT',     hc_cc_rs_path_type.TcT,      ('length', 'cstart', 'cend', 'qi1'), False),
        ('TcTcT',   hc_cc_rs_path_type.TcTcT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TcTT',    hc_cc_rs_path_type.TcTT,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TTcT',    hc_cc_rs_path_type.TTcT,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1'), False),
        ('TST',     hc_cc_rs_path_type.TST,      ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'qi4'), True),
        ('TSTcT',   hc_cc_rs_path_type.TSTcT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'qi4', 'ci1'), True),
        ('TcTST',   hc_cc_rs_path_type.TcTST,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'qi4', 'ci1'), True),
        ('TcTSTcT', hc_cc_rs_path_type.TcTSTcT,  ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'qi4', 'ci1', 'ci2'), True),
        ('TTcTT',   hc_cc_rs_path_type.TTcTT,    ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1', 'ci2'), False),
        ('TcTTcT',  hc_cc_rs_path_type.TcTTcT,   ('length', 'cstart', 'cend', 'qi1', 'qi2', 'ci1', 'ci2'), False),
        ('TTT',     hc_cc_rs_path_type.TTT,      ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3', 'ci1'), False),
        ('TcST',    hc_cc_rs_path_type.TcST,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3'), True),
        ('TScT',    hc_cc_rs_path_type.TScT,     ('length', 'cstart', 'cend', 'qi1', 'qi2', 'qi3'), True),
        ('TcScT',   hc_cc_rs_path_type.TcScT,    ('length', 'cstart', 'cend', 'qi1', 'qi2'), True),
    ]

    def __init__(self, kappa, sigma, discretization=0.1):
        import sys
        super().__init__(kappa, sigma, discretization)
        self._helper = _HCpmpm_Reeds_Shepp(self)
        self.rs_circle_param_ = HC_CC_Circle_Param()
        self.rs_circle_param_.set_param(
            self.kappa_, sys.float_info.max,
            1.0 / self.kappa_, 0.0, 0.0, 1.0, 0.0)
        self.radius_ = self.hc_cc_circle_param_.radius
        self.mu_ = self.hc_cc_circle_param_.mu
        self.sin_mu_ = self.hc_cc_circle_param_.sin_mu
        self.cos_mu_ = self.hc_cc_circle_param_.cos_mu

    def hcpmpm_circles_rs_path(self, c1, c2):
        """Evaluate all families and return the shortest HC_CC_RS_Path."""
        tp = hc_cc_rs_path_type
        h = self._helper
        h.distance = center_distance(c1, c2)
        h.angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)

        N = nb_hc_cc_rs_paths
        length = [_INF] * N
        qi1 = [None] * N
        qi2 = [None] * N
        qi3 = [None] * N
        qi4 = [None] * N
        cstart = [None] * N
        cend = [None] * N
        ci1 = [None] * N
        ci2 = [None] * N

        skip = False

        # case E
        if configuration_equal(c1.start, c2.start):
            length[tp.E] = 0
            skip = True

        # case T
        if not skip and configuration_on_hc_cc_circle(c1, c2.start):
            cstart[tp.T] = HC_CC_Circle(
                c1.start, c1.left, c1.forward, True, self.rs_circle_param_)
            length[tp.T] = cstart[tp.T].rs_turn_length(c2.start)
            skip = True

        if not skip:
            length, qi1, qi2, qi3, qi4, cstart, cend, ci1, ci2 = \
                _evaluate_rs_families(h, c1, c2, self._FAMILY_REGISTRY)

        # select shortest
        best = min(range(N), key=lambda i: length[i])
        return HC_CC_RS_Path(
            c1.start, c2.start, hc_cc_rs_path_type(best),
            self.kappa_, self.sigma_,
            qi1[best], qi2[best], qi3[best], qi4[best],
            cstart[best], cend[best], ci1[best], ci2[best],
            length[best])

    def _find_path(self, state1, state2):
        """Compute shortest path between two states (4×4 circle combos)."""
        start1 = Configuration(state1.x, state1.y, state1.theta, self.kappa_)
        start2 = Configuration(state1.x, state1.y, state1.theta, -self.kappa_)
        end1 = Configuration(state2.x, state2.y, state2.theta, self.kappa_)
        end2 = Configuration(state2.x, state2.y, state2.theta, -self.kappa_)
        rsp = self.rs_circle_param_
        start_circles = [
            HC_CC_Circle(start1, True, True, True, rsp),
            HC_CC_Circle(start2, False, True, True, rsp),
            HC_CC_Circle(start1, True, False, True, rsp),
            HC_CC_Circle(start2, False, False, True, rsp),
        ]
        end_circles = [
            HC_CC_Circle(end1, True, True, True, rsp),
            HC_CC_Circle(end2, False, True, True, rsp),
            HC_CC_Circle(end1, True, False, True, rsp),
            HC_CC_Circle(end2, False, False, True, rsp),
        ]
        best_path = None
        for i, sc in enumerate(start_circles):
            if i == 0 and state1.kappa < 0:
                continue
            if i == 1 and state1.kappa > 0:
                continue
            if i == 2 and state1.kappa < 0:
                continue
            if i == 3 and state1.kappa > 0:
                continue
            for j, ec in enumerate(end_circles):
                if j == 0 and state2.kappa < 0:
                    continue
                if j == 1 and state2.kappa > 0:
                    continue
                if j == 2 and state2.kappa < 0:
                    continue
                if j == 3 and state2.kappa > 0:
                    continue
                path = self.hcpmpm_circles_rs_path(sc, ec)
                if path is not None and (best_path is None
                                         or path.length < best_path.length):
                    best_path = path
        return best_path


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

    def _select_sub_space(self, start_state, end_state):
        """Select appropriate sub-state-space based on endpoint curvatures."""
        eps = get_epsilon()
        if abs(start_state.kappa) < eps:
            return self._hc0pm if abs(end_state.kappa) >= eps else self._hc00
        else:
            return self._hcpmpm if abs(end_state.kappa) >= eps else self._hcpm0

    def get_distance(self, state1, state2):
        start_preds = self.predict_state(state1)
        end_preds = self.predict_state(state2)
        best = _INF
        eps = get_epsilon()

        for s_state, s_ctrl in start_preds:
            for e_state, e_ctrl in end_preds:
                if state_equal(s_state, e_state):
                    ctrl = subtract_control(s_ctrl, e_ctrl)
                    best = min(best, abs(ctrl.delta_s))
                else:
                    dist = self._select_sub_space(s_state, e_state).get_distance(s_state, e_state)
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
        eps = get_epsilon()

        for s_state, s_ctrl in start_preds:
            for e_state, e_ctrl in end_preds:
                controls = []
                if state_equal(s_state, e_state):
                    ctrl = subtract_control(s_ctrl, e_ctrl)
                    controls.append(ctrl)
                else:
                    controls = self._select_sub_space(s_state, e_state).get_controls(s_state, e_state)
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
