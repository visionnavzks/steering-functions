"""CC00 Reeds-Shepp state space (Python port).

Ported from cc00_reeds_shepp_state_space.cpp.  CC00 has zero curvature at
both start **and** end, so every turn is a CC (clothoid-circle-clothoid)
turn and all circles are created with ``CC_REGULAR = False``.
"""

import math
from steering_functions.hc_cc_state_space import HC_CC_StateSpace
from steering_functions.state import State, Control
from steering_functions.configuration import (
    Configuration, configuration_distance, configuration_equal,
    configuration_aligned,
)
from steering_functions.hc_cc_circle import (
    HC_CC_Circle, center_distance, configuration_on_hc_cc_circle,
)
from steering_functions.paths import (
    HC_CC_RS_Path, hc_cc_rs_path_type,
    empty_controls, straight_controls, cc_turn_controls,
)
from steering_functions.utilities import (
    get_epsilon, global_frame_change,
    PI, HALF_PI,
)
from steering_functions.hc_reeds_shepp_state_space import (
    _evaluate_rs_families, _select_rs_best, _HC00_LAYOUTS,
)

_INF = float("inf")
_CC_REGULAR = False


# ===================================================================
# Inner helper – evaluates all CC00 path families between two circles
# ===================================================================
class _CC00_Reeds_Shepp:

    def __init__(self, parent):
        self._parent = parent
        self.distance = 0.0
        self.angle = 0.0

    def _p(self):
        return self._parent.hc_cc_circle_param_

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
        return abs(self.distance - 2 * c1.radius * c1.cos_mu) < get_epsilon()

    def TcT_tangent_circles(self, c1, c2):
        dist = center_distance(c1, c2)
        delta_x = 0.5 * dist
        delta_y = math.sqrt(max(0.0, c1.radius ** 2 - delta_x ** 2))
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        if c1.left:
            theta = angle + HALF_PI
            if c1.forward:
                x, y = global_frame_change(c1.xc, c1.yc, angle, delta_x, delta_y)
            else:
                x, y = global_frame_change(c1.xc, c1.yc, angle, delta_x, -delta_y)
        else:
            theta = angle - HALF_PI
            if c1.forward:
                x, y = global_frame_change(c1.xc, c1.yc, angle, delta_x, -delta_y)
            else:
                x, y = global_frame_change(c1.xc, c1.yc, angle, delta_x, delta_y)
        return Configuration(x, y, theta, 0)

    def TcT_path(self, c1, c2):
        q = self.TcT_tangent_circles(c1, c2)
        p = self._p()
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = cstart.cc_turn_length(q) + cend.cc_turn_length(q)
        return length, cstart, cend, q

    # ---- TcTcT ----------------------------------------------------------
    def TcTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance <= 4 * c1.radius * c1.cos_mu

    def TcTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r = 2 * c1.radius * c1.cos_mu
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
        middle1 = HC_CC_Circle(qa, not c1.left, not c1.forward, c1.regular, p)
        middle2 = HC_CC_Circle(qc, not c1.left, not c1.forward, c1.regular, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        l1 = (cstart.cc_turn_length(qa) + middle1.cc_turn_length(qb) +
              cend.cc_turn_length(qb))
        l2 = (cstart.cc_turn_length(qc) + middle2.cc_turn_length(qd) +
              cend.cc_turn_length(qd))
        if l1 < l2:
            return l1, cstart, cend, qa, qb, middle1
        return l2, cstart, cend, qc, qd, middle2

    # ---- TcTT -----------------------------------------------------------
    def TcTT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return ((self.distance <= 2 * c1.radius * (1 + c1.cos_mu)) and
                (self.distance >= 2 * c1.radius * (1 - c1.cos_mu)))

    def TcTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * c1.radius * c1.cos_mu
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
        middle1 = HC_CC_Circle(qa, not c1.left, not c1.forward, c1.regular, p)
        middle2 = HC_CC_Circle(qc, not c1.left, not c1.forward, c1.regular, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        l1 = (cstart.cc_turn_length(qa) + middle1.cc_turn_length(qb) +
              cend.cc_turn_length(qb))
        l2 = (cstart.cc_turn_length(qc) + middle2.cc_turn_length(qd) +
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
        return ((self.distance <= 2 * c1.radius * (1 + c1.cos_mu)) and
                (self.distance >= 2 * c1.radius * (1 - c1.cos_mu)))

    def TTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * c1.radius
        r2 = 2 * c1.radius * c1.cos_mu
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
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, c1.regular, p)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, c1.regular, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        l1 = (cstart.cc_turn_length(qa) + middle1.cc_turn_length(qb) +
              cend.cc_turn_length(qb))
        l2 = (cstart.cc_turn_length(qc) + middle2.cc_turn_length(qd) +
              cend.cc_turn_length(qd))
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
        return self.distance >= (
            2 * c1.radius *
            math.sqrt(1 + 2 * c1.sin_mu * c1.cos_mu + c1.cos_mu ** 2))

    def TeSTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= 2 * c1.radius * (c1.cos_mu + c1.sin_mu)

    def TSTcT_exists(self, c1, c2):
        return self.TiSTcT_exists(c1, c2) or self.TeSTcT_exists(c1, c2)

    def TiSTcT_path(self, c1, c2):
        theta = self.angle
        r = c2.radius * c2.cos_mu
        delta_y = (2 * r) ** 2 / self.distance
        delta_x = 2 * r * math.sqrt(max(0.0, 1 - delta_y / self.distance))
        p = self._p()
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1, q2 = self.TiST_tangent_circles(c1, tgt1)
        q3 = self.TcT_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        ci = HC_CC_Circle(q2, not c1.left, c1.forward, c1.regular, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  ci.cc_turn_length(q3) + cend.cc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3, ci

    def TeSTcT_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * c2.radius * c2.cos_mu
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1, q2 = self.TeST_tangent_circles(c1, tgt1)
        q3 = self.TcT_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        ci = HC_CC_Circle(q2, c1.left, c1.forward, c1.regular, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  ci.cc_turn_length(q3) + cend.cc_turn_length(q3))
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
        return self.distance >= (
            2 * c1.radius *
            math.sqrt(1 + 2 * c1.sin_mu * c1.cos_mu + c1.cos_mu ** 2))

    def TcTeST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= 2 * c1.radius * (c1.cos_mu + c1.sin_mu)

    def TcTST_exists(self, c1, c2):
        return self.TcTiST_exists(c1, c2) or self.TcTeST_exists(c1, c2)

    def TcTiST_path(self, c1, c2):
        theta = self.angle
        r = c1.radius * c1.cos_mu
        delta_y = (2 * r) ** 2 / self.distance
        delta_x = 2 * r * math.sqrt(max(0.0, 1 - delta_y / self.distance))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt1 = HC_CC_Circle(x, y, not c2.left, not c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TiST_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        ci = HC_CC_Circle(q1, not c1.left, not c1.forward, c1.regular, p)
        length = (cstart.cc_turn_length(q1) + ci.cc_turn_length(q2) +
                  configuration_distance(q2, q3) + cend.cc_turn_length(q3))
        return length, cstart, cend, q1, q2, q3, ci

    def TcTeST_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * c2.radius * c2.cos_mu
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, c2.left, not c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TeST_tangent_circles(tgt1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        ci = HC_CC_Circle(q1, not c1.left, not c1.forward, c1.regular, p)
        length = (cstart.cc_turn_length(q1) + ci.cc_turn_length(q2) +
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
        return self.distance >= (
            2 * c1.radius *
            math.sqrt(1 + 4 * c1.cos_mu * c1.sin_mu + 4 * c1.cos_mu ** 2))

    def TcTeSTcT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 2 * c1.radius * (2 * c1.cos_mu + c1.sin_mu)

    def TcTSTcT_exists(self, c1, c2):
        return self.TcTiSTcT_exists(c1, c2) or self.TcTeSTcT_exists(c1, c2)

    def TcTiSTcT_path(self, c1, c2):
        theta = self.angle
        r = c1.radius * c1.cos_mu
        delta_y = (2 * r) ** 2 / self.distance
        delta_x = 2 * r * math.sqrt(max(0.0, 1 - delta_y / self.distance))
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TiST_tangent_circles(tgt1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        ci1 = HC_CC_Circle(q1, not c1.left, not c1.forward, c1.regular, p)
        ci2 = HC_CC_Circle(q3, not c2.left, c2.forward, c2.regular, p)
        length = (cstart.cc_turn_length(q1) + ci1.cc_turn_length(q2) +
                  configuration_distance(q2, q3) + ci2.cc_turn_length(q4) +
                  cend.cc_turn_length(q4))
        return length, cstart, cend, q1, q2, q3, q4, ci1, ci2

    def TcTeSTcT_path(self, c1, c2):
        theta = self.angle
        delta_x = 2 * c1.radius * c1.cos_mu
        delta_y = 0.0
        p = self._p()
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, not c1.forward, c1.regular, p)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, c2.forward, c2.regular, p)
        q1 = self.TcT_tangent_circles(c1, tgt1)
        q2, q3 = self.TeST_tangent_circles(tgt1, tgt2)
        q4 = self.TcT_tangent_circles(tgt2, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        ci1 = HC_CC_Circle(q1, not c1.left, not c1.forward, c1.regular, p)
        ci2 = HC_CC_Circle(q3, not c2.left, c2.forward, c2.regular, p)
        length = (cstart.cc_turn_length(q1) + ci1.cc_turn_length(q2) +
                  configuration_distance(q2, q3) + ci2.cc_turn_length(q4) +
                  cend.cc_turn_length(q4))
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
        return self.distance <= 2 * c1.radius * (c1.cos_mu + 2)

    def TTcTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * c1.radius * c1.cos_mu
        r2 = 2 * c1.radius
        p = self._p()
        if self.distance < 2 * c1.radius * (-c1.cos_mu + 2):
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
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, c1.regular, p)
        middle2 = HC_CC_Circle(qb, not c2.left, not c2.forward, c2.regular, p)
        middle3 = HC_CC_Circle(qd, not c1.left, c1.forward, c1.regular, p)
        middle4 = HC_CC_Circle(qe, not c2.left, not c2.forward, c2.regular, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        l1 = (cstart.cc_turn_length(qa) + middle1.cc_turn_length(qb) +
              middle2.cc_turn_length(qc) + cend.cc_turn_length(qc))
        l2 = (cstart.cc_turn_length(qd) + middle3.cc_turn_length(qe) +
              middle4.cc_turn_length(qf) + cend.cc_turn_length(qf))
        if l1 < l2:
            return l1, cstart, cend, qa, qb, qc, middle1, middle2
        return l2, cstart, cend, qd, qe, qf, middle3, middle4

    # ---- TcTTcT ---------------------------------------------------------
    def TcTTcT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return ((self.distance <= 2 * c1.radius * (2 * c1.cos_mu + 1)) and
                (self.distance >= 2 * c1.radius * (2 * c1.cos_mu - 1)))

    def TcTTcT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * c1.radius * c1.cos_mu
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
        middle1 = HC_CC_Circle(qa, not c1.left, not c1.forward, c1.regular, p)
        middle2 = HC_CC_Circle(qb, c1.left, not c1.forward, c1.regular, p)
        middle3 = HC_CC_Circle(qd, not c1.left, not c1.forward, c1.regular, p)
        middle4 = HC_CC_Circle(qe, c1.left, not c1.forward, c1.regular, p)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        l1 = (cstart.cc_turn_length(qa) + middle1.cc_turn_length(qb) +
              middle2.cc_turn_length(qc) + cend.cc_turn_length(qc))
        l2 = (cstart.cc_turn_length(qd) + middle3.cc_turn_length(qe) +
              middle4.cc_turn_length(qf) + cend.cc_turn_length(qf))
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
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, c1.regular, p)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, c1.regular, p)
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
        return self.distance >= 2 * c1.radius * c1.cos_mu

    def TceST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= get_epsilon()

    def TcST_exists(self, c1, c2):
        return self.TciST_exists(c1, c2) or self.TceST_exists(c1, c2)

    def TciST_path(self, c1, c2):
        alpha = math.asin(2 * c1.radius * c1.cos_mu / self.distance)
        dx = c1.radius * c1.sin_mu
        dy = c1.radius * c1.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TceST_path(self, c1, c2):
        theta = self.angle
        dx = c1.radius * c1.sin_mu
        dy = c1.radius * c1.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -dx, dy)
            q2 = Configuration(x, y, theta, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
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
        return self.distance >= 2 * c1.radius * c1.cos_mu

    def TeScT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward != c2.forward:
            return False
        return self.distance >= get_epsilon()

    def TScT_exists(self, c1, c2):
        return self.TiScT_exists(c1, c2) or self.TeScT_exists(c1, c2)

    def TiScT_path(self, c1, c2):
        alpha = math.asin(2 * c1.radius * c1.cos_mu / self.distance)
        dx = c1.radius * c1.sin_mu
        dy = c1.radius * c1.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, 0)
        elif c1.left and not c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif not c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        else:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TeScT_path(self, c1, c2):
        theta = self.angle
        dx = c1.radius * c1.sin_mu
        dy = c1.radius * c1.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        elif c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, 0)
        else:
            x, y = global_frame_change(c1.xc, c1.yc, theta, dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
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
        return self.distance >= 2 * c1.radius * c1.cos_mu

    def TceScT_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= get_epsilon()

    def TcScT_exists(self, c1, c2):
        return self.TciScT_exists(c1, c2) or self.TceScT_exists(c1, c2)

    def TciScT_path(self, c1, c2):
        alpha = math.asin(2 * c1.radius * c1.cos_mu / self.distance)
        dx = c1.radius * c1.sin_mu
        dy = c1.radius * c1.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            theta = self.angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            theta = self.angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TceScT_path(self, c1, c2):
        theta = self.angle
        dx = c1.radius * c1.sin_mu
        dy = c1.radius * c1.cos_mu
        p = self._p()
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta + PI, 0)
        elif c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta, 0)
        elif not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, -dy)
            q2 = Configuration(x, y, theta + PI, 0)
        else:
            x, y = global_frame_change(c1.xc, c1.yc, theta, -dx, dy)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, dx, dy)
            q2 = Configuration(x, y, theta, 0)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, _CC_REGULAR, p)
        length = (cstart.cc_turn_length(q1) + configuration_distance(q1, q2) +
                  cend.cc_turn_length(q2))
        return length, cstart, cend, q1, q2

    def TcScT_path(self, c1, c2):
        if self.TciScT_exists(c1, c2):
            return self.TciScT_path(c1, c2)
        if self.TceScT_exists(c1, c2):
            return self.TceScT_path(c1, c2)
        return None


# ===================================================================
# Public class
# ===================================================================
class CC00_Reeds_Shepp_State_Space(HC_CC_StateSpace):
    """CC Reeds-Shepp with zero curvature at start and end (G² continuous)."""

    def __init__(self, kappa, sigma, discretization=0.1):
        super().__init__(kappa, sigma, discretization)
        self._helper = _CC00_Reeds_Shepp(self)

    # ------------------------------------------------------------------
    def cc00_circles_rs_path(self, c1, c2):
        """Evaluate all families for one (start-circle, end-circle) pair."""
        tp = hc_cc_rs_path_type
        h = self._helper
        h.distance = center_distance(c1, c2)
        h.angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        p = self.hc_cc_circle_param_

        # --- special cases ---
        if configuration_equal(c1.start, c2.start):
            return HC_CC_RS_Path(c1.start, c2.start, tp.E, self.kappa_, self.sigma_,
                                 None, None, None, None, None, None, None, None, 0)
        if configuration_aligned(c1.start, c2.start):
            d = configuration_distance(c1.start, c2.start)
            return HC_CC_RS_Path(c1.start, c2.start, tp.S, self.kappa_, self.sigma_,
                                 None, None, None, None, None, None, None, None, d)
        if configuration_aligned(c2.start, c1.start):
            d = configuration_distance(c2.start, c1.start)
            return HC_CC_RS_Path(c1.start, c2.start, tp.S, self.kappa_, self.sigma_,
                                 None, None, None, None, None, None, None, None, d)
        if configuration_on_hc_cc_circle(c1, c2.start):
            cs = HC_CC_Circle(c1.start, c1.left, c1.forward, _CC_REGULAR, p)
            return HC_CC_RS_Path(c1.start, c2.start, tp.T, self.kappa_, self.sigma_,
                                 None, None, None, None, cs, None, None, None,
                                 cs.cc_turn_length(c2.start))

        # --- evaluate all families ---
        length, arrays = _evaluate_rs_families(h, c1, c2, _HC00_LAYOUTS)
        return _select_rs_best(c1.start, c2.start, self.kappa_, self.sigma_,
                               length, arrays)

    # ------------------------------------------------------------------
    def cc00_reeds_shepp(self, state1, state2):
        """Compute shortest path between two states (4×4 circle combos)."""
        start = Configuration(state1.x, state1.y, state1.theta, 0.0)
        end = Configuration(state2.x, state2.y, state2.theta, 0.0)
        p = self.hc_cc_circle_param_
        start_circles = [
            HC_CC_Circle(start, True, True, _CC_REGULAR, p),
            HC_CC_Circle(start, False, True, _CC_REGULAR, p),
            HC_CC_Circle(start, True, False, _CC_REGULAR, p),
            HC_CC_Circle(start, False, False, _CC_REGULAR, p),
        ]
        end_circles = [
            HC_CC_Circle(end, True, True, _CC_REGULAR, p),
            HC_CC_Circle(end, False, True, _CC_REGULAR, p),
            HC_CC_Circle(end, True, False, _CC_REGULAR, p),
            HC_CC_Circle(end, False, False, _CC_REGULAR, p),
        ]
        best_path = None
        for sc in start_circles:
            for ec in end_circles:
                path = self.cc00_circles_rs_path(sc, ec)
                if path is not None and (best_path is None or path.length < best_path.length):
                    best_path = path
        return best_path

    # ------------------------------------------------------------------
    def get_distance(self, state1, state2):
        path = self.cc00_reeds_shepp(state1, state2)
        return path.length

    # ------------------------------------------------------------------
    def get_controls(self, state1, state2):
        tp = hc_cc_rs_path_type
        path = self.cc00_reeds_shepp(state1, state2)
        controls = []
        t = path.type

        if t == tp.E:
            empty_controls(controls)
        elif t == tp.S:
            straight_controls(path.start, path.end, controls)
        elif t == tp.T:
            cc_turn_controls(path.cstart, path.end, True, controls)
        elif t in (tp.TT, tp.TcT):
            cc_turn_controls(path.cstart, path.qi1, True, controls)
            cc_turn_controls(path.cend, path.qi1, False, controls)
        elif t in (tp.TcTcT, tp.TcTT, tp.TTcT):
            cc_turn_controls(path.cstart, path.qi1, True, controls)
            cc_turn_controls(path.ci1, path.qi2, True, controls)
            cc_turn_controls(path.cend, path.qi2, False, controls)
        elif t == tp.TST:
            cc_turn_controls(path.cstart, path.qi1, True, controls)
            straight_controls(path.qi1, path.qi2, controls)
            cc_turn_controls(path.cend, path.qi2, False, controls)
        elif t == tp.TSTcT:
            cc_turn_controls(path.cstart, path.qi1, True, controls)
            straight_controls(path.qi1, path.qi2, controls)
            cc_turn_controls(path.ci1, path.qi3, True, controls)
            cc_turn_controls(path.cend, path.qi3, False, controls)
        elif t == tp.TcTST:
            cc_turn_controls(path.cstart, path.qi1, True, controls)
            cc_turn_controls(path.ci1, path.qi2, True, controls)
            straight_controls(path.qi2, path.qi3, controls)
            cc_turn_controls(path.cend, path.qi3, False, controls)
        elif t == tp.TcTSTcT:
            cc_turn_controls(path.cstart, path.qi1, True, controls)
            cc_turn_controls(path.ci1, path.qi2, True, controls)
            straight_controls(path.qi2, path.qi3, controls)
            cc_turn_controls(path.ci2, path.qi4, True, controls)
            cc_turn_controls(path.cend, path.qi4, False, controls)
        elif t in (tp.TTcTT, tp.TcTTcT):
            cc_turn_controls(path.cstart, path.qi1, True, controls)
            cc_turn_controls(path.ci1, path.qi2, True, controls)
            cc_turn_controls(path.ci2, path.qi3, True, controls)
            cc_turn_controls(path.cend, path.qi3, False, controls)
        elif t == tp.TTT:
            cc_turn_controls(path.cstart, path.qi1, True, controls)
            cc_turn_controls(path.ci1, path.qi2, True, controls)
            cc_turn_controls(path.cend, path.qi2, False, controls)
        elif t in (tp.TcST, tp.TScT, tp.TcScT):
            cc_turn_controls(path.cstart, path.qi1, True, controls)
            straight_controls(path.qi1, path.qi2, controls)
            cc_turn_controls(path.cend, path.qi2, False, controls)

        return controls
