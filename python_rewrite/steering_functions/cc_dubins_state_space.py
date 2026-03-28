"""
CC Dubins state space module.

Ports of the C++ CC Dubins state space classes for continuous curvature
Dubins path planning. Contains five classes:
  - CC00_Dubins_State_Space: zero curvature at start and end
  - CC0pm_Dubins_State_Space: zero start, ±max end curvature
  - CCpm0_Dubins_State_Space: ±max start, zero end curvature
  - CCpmpm_Dubins_State_Space: ±max curvature at both start and end
  - CC_Dubins_State_Space: general wrapper selecting best among all 4

Derived from Continuous Curvature (CC) Steer.
Copyright (c) 2016, Thierry Fraichard and Inria.
Copyright (c) 2017, see NOTICE file.
"""

import math

from steering_functions.hc_cc_state_space import HC_CC_StateSpace
from steering_functions.state import State, Control
from steering_functions.configuration import (
    Configuration,
    configuration_distance,
    configuration_equal,
    configuration_aligned,
)
from steering_functions.hc_cc_circle import (
    HC_CC_Circle,
    HC_CC_Circle_Param,
    center_distance,
    configuration_on_hc_cc_circle,
)
from steering_functions.paths import (
    CC_Dubins_Path,
    cc_dubins_path_type,
    nb_cc_dubins_paths,
    empty_controls,
    straight_controls,
    cc_turn_controls,
    rs_turn_controls,
    hc_turn_controls,
    state_equal,
    subtract_control,
    reverse_control,
    _build_controls,
)
from steering_functions.utilities import (
    get_epsilon,
    sgn,
    double_array_init,
    array_index_min,
    global_frame_change,
    end_of_clothoid,
    PI,
    HALF_PI,
)

_INF = float("inf")


# =============================================================================
# CC00 Helper
# =============================================================================
class _CC00_Dubins:
    """Helper for computing CC Dubins paths between two CC circles
    with zero curvature at both start and end."""

    def __init__(self, parent):
        self._parent = parent
        self.distance = 0.0
        self.angle = 0.0

    # ----- TT ----------------------------------------------------------------
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
        length = c1.cc_turn_length(q) + c2.cc_turn_length(q)
        return length, q

    # ----- TST ---------------------------------------------------------------
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
        distance = center_distance(c1, c2)
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        alpha = math.asin(2 * c1.radius * c1.cos_mu / distance)
        delta_x = c1.radius * c1.sin_mu
        delta_y = c1.radius * c1.cos_mu
        q1 = q2 = None
        if c1.left and c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta, 0)
        if c1.left and not c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        if not c1.left and c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta, 0)
        if not c1.left and not c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TeST_tangent_circles(self, c1, c2):
        delta_x = c1.radius * c1.sin_mu
        delta_y = c1.radius * c1.cos_mu
        theta = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        q1 = q2 = None
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta, 0)
        if c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        if not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta, 0)
        if not c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TiST_path(self, c1, c2):
        q1, q2 = self.TiST_tangent_circles(c1, c2)
        length = (
            c1.cc_turn_length(q1)
            + configuration_distance(q1, q2)
            + c2.cc_turn_length(q2)
        )
        return length, q1, q2

    def TeST_path(self, c1, c2):
        q1, q2 = self.TeST_tangent_circles(c1, c2)
        length = (
            c1.cc_turn_length(q1)
            + configuration_distance(q1, q2)
            + c2.cc_turn_length(q2)
        )
        return length, q1, q2

    def TST_path(self, c1, c2):
        if self.TiST_exists(c1, c2):
            return self.TiST_path(c1, c2)
        if self.TeST_exists(c1, c2):
            return self.TeST_path(c1, c2)
        return _INF, None, None

    # ----- TTT ---------------------------------------------------------------
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
        delta_y = math.sqrt(r ** 2 - delta_x ** 2)

        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular,
                            self._parent.hc_cc_circle_param_)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular,
                            self._parent.hc_cc_circle_param_)

        qa = self.TT_tangent_circles(c1, tgt1)
        qb = self.TT_tangent_circles(tgt1, c2)
        qc = self.TT_tangent_circles(c1, tgt2)
        qd = self.TT_tangent_circles(tgt2, c2)
        return qa, qb, qc, qd

    def TTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTT_tangent_circles(c1, c2)
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, c1.regular,
                               self._parent.hc_cc_circle_param_)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, c1.regular,
                               self._parent.hc_cc_circle_param_)

        length1 = (
            c1.cc_turn_length(qa)
            + middle1.cc_turn_length(qb)
            + c2.cc_turn_length(qb)
        )
        length2 = (
            c1.cc_turn_length(qc)
            + middle2.cc_turn_length(qd)
            + c2.cc_turn_length(qd)
        )
        if length1 < length2:
            return length1, qa, qb, middle1
        else:
            return length2, qc, qd, middle2


# =============================================================================
# CC00_Dubins_State_Space
# =============================================================================
class CC00_Dubins_State_Space(HC_CC_StateSpace):
    """CC Dubins state space with zero curvature at start and end."""

    _CONTROL_TABLE = {
        cc_dubins_path_type.E: [('empty',)],
        cc_dubins_path_type.S: [('straight', 'start', 'end')],
        cc_dubins_path_type.T: [('cc', 'cstart', 'end', True)],
        cc_dubins_path_type.TT: [('cc', 'cstart', 'qi1', True), ('cc', 'cend', 'qi1', False)],
        cc_dubins_path_type.TST: [('cc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('cc', 'cend', 'qi2', False)],
        cc_dubins_path_type.TTT: [('cc', 'cstart', 'qi1', True), ('cc', 'ci1', 'qi2', True), ('cc', 'cend', 'qi2', False)],
    }

    def __init__(self, kappa, sigma, discretization=0.1, forwards=True):
        super().__init__(kappa, sigma, discretization)
        self.forwards_ = forwards
        self._cc00_dubins = _CC00_Dubins(self)

    def cc00_circles_dubins_path(self, c1, c2):
        length = double_array_init(nb_cc_dubins_paths, _INF)
        qi1 = [None] * nb_cc_dubins_paths
        qi2 = [None] * nb_cc_dubins_paths
        cstart = [None] * nb_cc_dubins_paths
        ci1 = [None] * nb_cc_dubins_paths
        cend = [None] * nb_cc_dubins_paths

        self._cc00_dubins.distance = center_distance(c1, c2)
        self._cc00_dubins.angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)

        E = cc_dubins_path_type.E
        S = cc_dubins_path_type.S
        T = cc_dubins_path_type.T
        TT = cc_dubins_path_type.TT
        TST = cc_dubins_path_type.TST
        TTT = cc_dubins_path_type.TTT

        # case E
        if configuration_equal(c1.start, c2.start):
            length[E] = 0
        # case S forwards
        elif self.forwards_ and configuration_aligned(c1.start, c2.start):
            length[S] = configuration_distance(c1.start, c2.start)
        # case S backwards
        elif not self.forwards_ and configuration_aligned(c2.start, c1.start):
            length[S] = configuration_distance(c2.start, c1.start)
        # case T
        elif configuration_on_hc_cc_circle(c1, c2.start):
            cstart[T] = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                                     self.hc_cc_circle_param_)
            length[T] = cstart[T].cc_turn_length(c2.start)
        else:
            # case TT
            if self._cc00_dubins.TT_exists(c1, c2):
                cstart[TT] = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                                          self.hc_cc_circle_param_)
                cend[TT] = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                                        self.hc_cc_circle_param_)
                l, q = self._cc00_dubins.TT_path(cstart[TT], cend[TT])
                length[TT] = l
                qi1[TT] = q
            # case TST
            if self._cc00_dubins.TST_exists(c1, c2):
                cstart[TST] = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                                           self.hc_cc_circle_param_)
                cend[TST] = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                                         self.hc_cc_circle_param_)
                l, q1, q2 = self._cc00_dubins.TST_path(cstart[TST], cend[TST])
                length[TST] = l
                qi1[TST] = q1
                qi2[TST] = q2
            # case TTT
            if self._cc00_dubins.TTT_exists(c1, c2):
                cstart[TTT] = HC_CC_Circle(c1.start, c1.left, c1.forward, c1.regular,
                                           self.hc_cc_circle_param_)
                cend[TTT] = HC_CC_Circle(c2.start, c2.left, c2.forward, c2.regular,
                                         self.hc_cc_circle_param_)
                l, q1, q2, mid = self._cc00_dubins.TTT_path(cstart[TTT], cend[TTT])
                length[TTT] = l
                qi1[TTT] = q1
                qi2[TTT] = q2
                ci1[TTT] = mid

        best = array_index_min(length)
        best_path = cc_dubins_path_type(best)
        return CC_Dubins_Path(
            c1.start, c2.start, best_path, self.kappa_, self.sigma_,
            qi1[best], qi2[best], None, None,
            cstart[best], cend[best], ci1[best], None,
            length[best],
        )

    def _find_path(self, state1, state2):
        start = Configuration(state1.x, state1.y, state1.theta, 0.0)
        end = Configuration(state2.x, state2.y, state2.theta, 0.0)

        if self.forwards_:
            start_circles = [
                HC_CC_Circle(start, True, True, True, self.hc_cc_circle_param_),
                HC_CC_Circle(start, False, True, True, self.hc_cc_circle_param_),
            ]
            end_circles = [
                HC_CC_Circle(end, True, False, True, self.hc_cc_circle_param_),
                HC_CC_Circle(end, False, False, True, self.hc_cc_circle_param_),
            ]
        else:
            start_circles = [
                HC_CC_Circle(start, True, False, True, self.hc_cc_circle_param_),
                HC_CC_Circle(start, False, False, True, self.hc_cc_circle_param_),
            ]
            end_circles = [
                HC_CC_Circle(end, True, True, True, self.hc_cc_circle_param_),
                HC_CC_Circle(end, False, True, True, self.hc_cc_circle_param_),
            ]

        paths = [None] * 4
        lg = [_INF] * 4
        for i in range(2):
            for j in range(2):
                idx = 2 * i + j
                paths[idx] = self.cc00_circles_dubins_path(start_circles[i], end_circles[j])
                if paths[idx] is not None:
                    lg[idx] = paths[idx].length

        best = array_index_min(lg)
        return paths[best]


# =============================================================================
# CC0pm Helper
# =============================================================================
class _CC0pm_Dubins:
    """Helper for computing CC Dubins paths between two circles
    with zero curvature at start and ±max curvature at end."""

    def __init__(self, parent):
        self._parent = parent
        self.distance = 0.0
        self.angle = 0.0

    # ----- TT ----------------------------------------------------------------
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
        q1 = self.TT_tangent_circles(c1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        cend = HC_CC_Circle(q1, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        q2 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = cstart.cc_turn_length(q1) + cend.hc_turn_length(q2)
        return length, cstart, cend, q1, q2

    # ----- TST ---------------------------------------------------------------
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
        distance = center_distance(c1, c2)
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        alpha = math.asin(2 * c1.radius * c1.cos_mu / distance)
        delta_x = c1.radius * c1.sin_mu
        delta_y = c1.radius * c1.cos_mu
        q1 = q2 = None
        if c1.left and c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta, 0)
        if c1.left and not c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        if not c1.left and c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta, 0)
        if not c1.left and not c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TeST_tangent_circles(self, c1, c2):
        delta_x = c1.radius * c1.sin_mu
        delta_y = c1.radius * c1.cos_mu
        theta = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        q1 = q2 = None
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta, 0)
        if c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        if not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta, 0)
        if not c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TiST_path(self, c1, c2):
        q1, q2 = self.TiST_tangent_circles(c1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        cend = HC_CC_Circle(q2, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = (
            cstart.cc_turn_length(q1)
            + configuration_distance(q1, q2)
            + cend.hc_turn_length(q3)
        )
        return length, cstart, cend, q1, q2, q3

    def TeST_path(self, c1, c2):
        q1, q2 = self.TeST_tangent_circles(c1, c2)
        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        cend = HC_CC_Circle(q2, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = (
            cstart.cc_turn_length(q1)
            + configuration_distance(q1, q2)
            + cend.hc_turn_length(q3)
        )
        return length, cstart, cend, q1, q2, q3

    def TST_path(self, c1, c2):
        if self.TiST_exists(c1, c2):
            return self.TiST_path(c1, c2)
        if self.TeST_exists(c1, c2):
            return self.TeST_path(c1, c2)
        return _INF, None, None, None, None, None

    # ----- TTT ---------------------------------------------------------------
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
        delta_y = math.sqrt(r ** 2 - delta_x ** 2)

        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular,
                            self._parent.hc_cc_circle_param_)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular,
                            self._parent.hc_cc_circle_param_)

        qa = self.TT_tangent_circles(c1, tgt1)
        qb = self.TT_tangent_circles(tgt1, c2)
        qc = self.TT_tangent_circles(c1, tgt2)
        qd = self.TT_tangent_circles(tgt2, c2)
        return qa, qb, qc, qd

    def TTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTT_tangent_circles(c1, c2)

        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True,
                               self._parent.hc_cc_circle_param_)
        end1 = HC_CC_Circle(qb, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, True,
                               self._parent.hc_cc_circle_param_)
        end2 = HC_CC_Circle(qd, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)

        cstart = HC_CC_Circle(c1.start, c1.left, c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)

        length1 = (
            cstart.cc_turn_length(qa)
            + middle1.cc_turn_length(qb)
            + end1.hc_turn_length(q3)
        )
        length2 = (
            cstart.cc_turn_length(qc)
            + middle2.cc_turn_length(qd)
            + end2.hc_turn_length(q3)
        )
        if length1 < length2:
            return length1, cstart, end1, qa, qb, q3, middle1
        else:
            return length2, cstart, end2, qc, qd, q3, middle2


# =============================================================================
# CC0pm_Dubins_State_Space
# =============================================================================
class CC0pm_Dubins_State_Space(HC_CC_StateSpace):
    """CC Dubins state space with zero curvature at start and ±max at end."""

    _CONTROL_TABLE = {
        cc_dubins_path_type.E: [('empty',)],
        cc_dubins_path_type.T: [('hc', 'cstart', 'end', True)],
        cc_dubins_path_type.TT: [('cc', 'cstart', 'qi1', True), ('hc', 'cend', 'qi2', True)],
        cc_dubins_path_type.TST: [('cc', 'cstart', 'qi1', True), ('straight', 'qi1', 'qi2'), ('hc', 'cend', 'qi3', True)],
        cc_dubins_path_type.TTT: [('cc', 'cstart', 'qi1', True), ('cc', 'ci1', 'qi2', True), ('hc', 'cend', 'qi3', True)],
    }

    def __init__(self, kappa, sigma, discretization=0.1, forwards=True):
        super().__init__(kappa, sigma, discretization)
        self.forwards_ = forwards
        self._cc0pm_dubins = _CC0pm_Dubins(self)
        self.rs_circle_param_ = HC_CC_Circle_Param()
        self.rs_circle_param_.set_param(kappa, _INF, 1.0 / kappa, 0.0, 0.0, 1.0, 0.0)
        self.radius_ = self.hc_cc_circle_param_.radius
        self.mu_ = self.hc_cc_circle_param_.mu

    def cc0pm_circles_dubins_path(self, c1, c2):
        length = double_array_init(nb_cc_dubins_paths, _INF)
        qi1 = [None] * nb_cc_dubins_paths
        qi2 = [None] * nb_cc_dubins_paths
        qi3 = [None] * nb_cc_dubins_paths
        cstart = [None] * nb_cc_dubins_paths
        ci1 = [None] * nb_cc_dubins_paths
        cend = [None] * nb_cc_dubins_paths

        self._cc0pm_dubins.distance = center_distance(c1, c2)
        self._cc0pm_dubins.angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)

        E = cc_dubins_path_type.E
        T = cc_dubins_path_type.T
        TT = cc_dubins_path_type.TT
        TST = cc_dubins_path_type.TST
        TTT = cc_dubins_path_type.TTT

        # case E
        if configuration_equal(c1.start, c2.start):
            length[E] = 0
        # case T
        elif self._cc0pm_dubins.distance < get_epsilon():
            cstart[T] = HC_CC_Circle(c1.start, c1.left, c1.forward, True,
                                     self.hc_cc_circle_param_)
            length[T] = cstart[T].hc_turn_length(c2.start)
        else:
            # case TT
            if self._cc0pm_dubins.TT_exists(c1, c2):
                l, cs, ce, q1, q2 = self._cc0pm_dubins.TT_path(c1, c2)
                length[TT] = l
                cstart[TT] = cs
                cend[TT] = ce
                qi1[TT] = q1
                qi2[TT] = q2
            # case TST
            if self._cc0pm_dubins.TST_exists(c1, c2):
                l, cs, ce, q1, q2, q3 = self._cc0pm_dubins.TST_path(c1, c2)
                length[TST] = l
                cstart[TST] = cs
                cend[TST] = ce
                qi1[TST] = q1
                qi2[TST] = q2
                qi3[TST] = q3
            # case TTT
            if self._cc0pm_dubins.TTT_exists(c1, c2):
                l, cs, ce, q1, q2, q3, mid = self._cc0pm_dubins.TTT_path(c1, c2)
                length[TTT] = l
                cstart[TTT] = cs
                cend[TTT] = ce
                qi1[TTT] = q1
                qi2[TTT] = q2
                qi3[TTT] = q3
                ci1[TTT] = mid

        best = array_index_min(length)
        best_path = cc_dubins_path_type(best)
        return CC_Dubins_Path(
            c1.start, c2.start, best_path, self.kappa_, self.sigma_,
            qi1[best], qi2[best], qi3[best], None,
            cstart[best], cend[best], ci1[best], None,
            length[best],
        )

    def _find_path(self, state1, state2):
        start = Configuration(state1.x, state1.y, state1.theta, 0.0)
        end1 = Configuration(state2.x, state2.y, state2.theta, self.kappa_)
        end2 = Configuration(state2.x, state2.y, state2.theta, -self.kappa_)

        if self.forwards_:
            start_circles = [
                HC_CC_Circle(start, True, True, True, self.hc_cc_circle_param_),
                HC_CC_Circle(start, False, True, True, self.hc_cc_circle_param_),
            ]
            end_circles = [
                HC_CC_Circle(end1, True, False, True, self.rs_circle_param_),
                HC_CC_Circle(end2, False, False, True, self.rs_circle_param_),
            ]
        else:
            start_circles = [
                HC_CC_Circle(start, True, False, True, self.hc_cc_circle_param_),
                HC_CC_Circle(start, False, False, True, self.hc_cc_circle_param_),
            ]
            end_circles = [
                HC_CC_Circle(end1, True, True, True, self.rs_circle_param_),
                HC_CC_Circle(end2, False, True, True, self.rs_circle_param_),
            ]

        paths = [None] * 4
        lg = [_INF] * 4
        for i in range(2):
            for j in range(2):
                # skip circle at the end for curvature continuity
                if j == 0 and state2.kappa < 0:
                    continue
                elif j == 1 and state2.kappa > 0:
                    continue
                idx = 2 * i + j
                paths[idx] = self.cc0pm_circles_dubins_path(
                    start_circles[i], end_circles[j]
                )
                if paths[idx] is not None:
                    lg[idx] = paths[idx].length

        best = array_index_min(lg)
        return paths[best]


# =============================================================================
# CCpm0 Helper
# =============================================================================
class _CCpm0_Dubins:
    """Helper for computing CC Dubins paths between two circles
    with ±max curvature at start and zero curvature at end."""

    def __init__(self, parent):
        self._parent = parent
        self.distance = 0.0
        self.angle = 0.0

    # ----- TT ----------------------------------------------------------------
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
        cstart = HC_CC_Circle(q2, not c2.left, c2.forward, True,
                              self._parent.hc_cc_circle_param_)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        length = cstart.hc_turn_length(q1) + cend.cc_turn_length(q2)
        return length, cstart, cend, q1, q2

    # ----- TST ---------------------------------------------------------------
    def TiST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 2 * c2.radius

    def TeST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 2 * c2.radius * c2.sin_mu

    def TST_exists(self, c1, c2):
        return self.TiST_exists(c1, c2) or self.TeST_exists(c1, c2)

    def TiST_tangent_circles(self, c1, c2):
        distance = center_distance(c1, c2)
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        alpha = math.asin(2 * c2.radius * c2.cos_mu / distance)
        delta_x = c2.radius * c2.sin_mu
        delta_y = c2.radius * c2.cos_mu
        q1 = q2 = None
        if c1.left and c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta, 0)
        if c1.left and not c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        if not c1.left and c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta, 0)
        if not c1.left and not c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TeST_tangent_circles(self, c1, c2):
        delta_x = c2.radius * c2.sin_mu
        delta_y = c2.radius * c2.cos_mu
        theta = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        q1 = q2 = None
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta, 0)
        if c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        if not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta, 0)
        if not c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TiST_path(self, c1, c2):
        q2, q3 = self.TiST_tangent_circles(c1, c2)
        cstart = HC_CC_Circle(q2, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        length = (
            cstart.hc_turn_length(q1)
            + configuration_distance(q2, q3)
            + cend.cc_turn_length(q3)
        )
        return length, cstart, cend, q1, q2, q3

    def TeST_path(self, c1, c2):
        q2, q3 = self.TeST_tangent_circles(c1, c2)
        cstart = HC_CC_Circle(q2, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        length = (
            cstart.hc_turn_length(q1)
            + configuration_distance(q2, q3)
            + cend.cc_turn_length(q3)
        )
        return length, cstart, cend, q1, q2, q3

    def TST_path(self, c1, c2):
        if self.TiST_exists(c1, c2):
            return self.TiST_path(c1, c2)
        if self.TeST_exists(c1, c2):
            return self.TeST_path(c1, c2)
        return _INF, None, None, None, None, None

    # ----- TTT ---------------------------------------------------------------
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
        delta_y = math.sqrt(r ** 2 - delta_x ** 2)

        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular,
                            self._parent.hc_cc_circle_param_)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular,
                            self._parent.hc_cc_circle_param_)

        qa = self.TT_tangent_circles(c1, tgt1)
        qb = self.TT_tangent_circles(tgt1, c2)
        qc = self.TT_tangent_circles(c1, tgt2)
        qd = self.TT_tangent_circles(tgt2, c2)
        return qa, qb, qc, qd

    def TTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTT_tangent_circles(c1, c2)

        start1 = HC_CC_Circle(qa, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True,
                               self._parent.hc_cc_circle_param_)
        start2 = HC_CC_Circle(qc, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, True,
                               self._parent.hc_cc_circle_param_)

        cend = HC_CC_Circle(c2.start, c2.left, c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)

        length1 = (
            start1.hc_turn_length(q1)
            + middle1.cc_turn_length(qb)
            + cend.cc_turn_length(qb)
        )
        length2 = (
            start2.hc_turn_length(q1)
            + middle2.cc_turn_length(qd)
            + cend.cc_turn_length(qd)
        )
        if length1 < length2:
            return length1, start1, cend, q1, qb, middle1
        else:
            return length2, start2, cend, q1, qd, middle2


# =============================================================================
# CCpm0_Dubins_State_Space
# =============================================================================
class CCpm0_Dubins_State_Space(HC_CC_StateSpace):
    """CC Dubins state space with ±max curvature at start and zero at end."""

    _CONTROL_TABLE = {
        cc_dubins_path_type.E: [('empty',)],
        cc_dubins_path_type.T: [('hc', 'cend', 'start', False)],
        cc_dubins_path_type.TT: [('hc', 'cstart', 'qi1', False), ('cc', 'cend', 'qi2', False)],
        cc_dubins_path_type.TST: [('hc', 'cstart', 'qi1', False), ('straight', 'qi2', 'qi3'), ('cc', 'cend', 'qi3', False)],
        cc_dubins_path_type.TTT: [('hc', 'cstart', 'qi1', False), ('cc', 'ci1', 'qi2', True), ('cc', 'cend', 'qi2', False)],
    }

    def __init__(self, kappa, sigma, discretization=0.1, forwards=True):
        super().__init__(kappa, sigma, discretization)
        self.forwards_ = forwards
        self._ccpm0_dubins = _CCpm0_Dubins(self)
        self.rs_circle_param_ = HC_CC_Circle_Param()
        self.rs_circle_param_.set_param(kappa, _INF, 1.0 / kappa, 0.0, 0.0, 1.0, 0.0)
        self.radius_ = self.hc_cc_circle_param_.radius
        self.mu_ = self.hc_cc_circle_param_.mu

    def ccpm0_circles_dubins_path(self, c1, c2):
        length = double_array_init(nb_cc_dubins_paths, _INF)
        qi1 = [None] * nb_cc_dubins_paths
        qi2 = [None] * nb_cc_dubins_paths
        qi3 = [None] * nb_cc_dubins_paths
        cstart = [None] * nb_cc_dubins_paths
        ci1 = [None] * nb_cc_dubins_paths
        cend = [None] * nb_cc_dubins_paths

        self._ccpm0_dubins.distance = center_distance(c1, c2)
        self._ccpm0_dubins.angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)

        E = cc_dubins_path_type.E
        T = cc_dubins_path_type.T
        TT = cc_dubins_path_type.TT
        TST = cc_dubins_path_type.TST
        TTT = cc_dubins_path_type.TTT

        # case E
        if configuration_equal(c1.start, c2.start):
            length[E] = 0
        # case T
        elif self._ccpm0_dubins.distance < get_epsilon():
            cend[T] = HC_CC_Circle(c2.start, c2.left, c2.forward, True,
                                   self.hc_cc_circle_param_)
            length[T] = cend[T].hc_turn_length(c1.start)
        else:
            # case TT
            if self._ccpm0_dubins.TT_exists(c1, c2):
                l, cs, ce, q1, q2 = self._ccpm0_dubins.TT_path(c1, c2)
                length[TT] = l
                cstart[TT] = cs
                cend[TT] = ce
                qi1[TT] = q1
                qi2[TT] = q2
            # case TST
            if self._ccpm0_dubins.TST_exists(c1, c2):
                l, cs, ce, q1, q2, q3 = self._ccpm0_dubins.TST_path(c1, c2)
                length[TST] = l
                cstart[TST] = cs
                cend[TST] = ce
                qi1[TST] = q1
                qi2[TST] = q2
                qi3[TST] = q3
            # case TTT
            if self._ccpm0_dubins.TTT_exists(c1, c2):
                l, cs, ce, q1, q2, mid = self._ccpm0_dubins.TTT_path(c1, c2)
                length[TTT] = l
                cstart[TTT] = cs
                cend[TTT] = ce
                qi1[TTT] = q1
                qi2[TTT] = q2
                ci1[TTT] = mid

        best = array_index_min(length)
        best_path = cc_dubins_path_type(best)
        return CC_Dubins_Path(
            c1.start, c2.start, best_path, self.kappa_, self.sigma_,
            qi1[best], qi2[best], qi3[best], None,
            cstart[best], cend[best], ci1[best], None,
            length[best],
        )

    def _find_path(self, state1, state2):
        start1 = Configuration(state1.x, state1.y, state1.theta, self.kappa_)
        start2 = Configuration(state1.x, state1.y, state1.theta, -self.kappa_)
        end = Configuration(state2.x, state2.y, state2.theta, 0.0)

        if self.forwards_:
            start_circles = [
                HC_CC_Circle(start1, True, True, True, self.rs_circle_param_),
                HC_CC_Circle(start2, False, True, True, self.rs_circle_param_),
            ]
            end_circles = [
                HC_CC_Circle(end, True, False, True, self.hc_cc_circle_param_),
                HC_CC_Circle(end, False, False, True, self.hc_cc_circle_param_),
            ]
        else:
            start_circles = [
                HC_CC_Circle(start1, True, False, True, self.rs_circle_param_),
                HC_CC_Circle(start2, False, False, True, self.rs_circle_param_),
            ]
            end_circles = [
                HC_CC_Circle(end, True, True, True, self.hc_cc_circle_param_),
                HC_CC_Circle(end, False, True, True, self.hc_cc_circle_param_),
            ]

        paths = [None] * 4
        lg = [_INF] * 4
        for i in range(2):
            # skip circle at the beginning for curvature continuity
            if i == 0 and state1.kappa < 0:
                continue
            elif i == 1 and state1.kappa > 0:
                continue
            for j in range(2):
                idx = 2 * i + j
                paths[idx] = self.ccpm0_circles_dubins_path(
                    start_circles[i], end_circles[j]
                )
                if paths[idx] is not None:
                    lg[idx] = paths[idx].length

        best = array_index_min(lg)
        return paths[best]


# =============================================================================
# CCpmpm Helper
# =============================================================================
class _CCpmpm_Dubins:
    """Helper for computing CC Dubins paths between two circles
    with ±max curvature at both start and end."""

    def __init__(self, parent):
        self._parent = parent
        self.distance = 0.0
        self.angle = 0.0

    # ----- TT ----------------------------------------------------------------
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
        if c1.left:
            theta = (
                angle + HALF_PI - self._parent.mu_
                if c1.forward
                else angle + HALF_PI + self._parent.mu_
            )
        else:
            theta = (
                angle - HALF_PI + self._parent.mu_
                if c1.forward
                else angle - HALF_PI - self._parent.mu_
            )
        return Configuration(x, y, theta, 0)

    def TT_path(self, c1, c2):
        q2 = self.TT_tangent_circles(c1, c2)
        cstart = HC_CC_Circle(q2, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        cend = HC_CC_Circle(q2, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = cstart.hc_turn_length(q1) + cend.hc_turn_length(q3)
        return length, cstart, cend, q1, q2, q3

    # ----- TST ---------------------------------------------------------------
    def TiST_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 2 * self._parent.radius_

    def TeST_exists(self, c1, c2):
        if c1.left != c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance >= 2 * self._parent.radius_ * self._parent.sin_mu_

    def TST_exists(self, c1, c2):
        return self.TiST_exists(c1, c2) or self.TeST_exists(c1, c2)

    def TiST_tangent_circles(self, c1, c2):
        distance = center_distance(c1, c2)
        angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        alpha = math.asin(2 * self._parent.radius_ * self._parent.cos_mu_ / distance)
        delta_x = self._parent.radius_ * self._parent.sin_mu_
        delta_y = self._parent.radius_ * self._parent.cos_mu_
        q1 = q2 = None
        if c1.left and c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta, 0)
        if c1.left and not c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        if not c1.left and c1.forward:
            theta = angle - alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta, 0)
        if not c1.left and not c1.forward:
            theta = angle + alpha
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TeST_tangent_circles(self, c1, c2):
        delta_x = self._parent.radius_ * self._parent.sin_mu_
        delta_y = self._parent.radius_ * self._parent.cos_mu_
        theta = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)
        q1 = q2 = None
        if c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta, 0)
        if c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        if not c1.left and c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
            q1 = Configuration(x, y, theta, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
            q2 = Configuration(x, y, theta, 0)
        if not c1.left and not c1.forward:
            x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
            q1 = Configuration(x, y, theta + PI, 0)
            x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
            q2 = Configuration(x, y, theta + PI, 0)
        return q1, q2

    def TiST_path(self, c1, c2):
        q2, q3 = self.TiST_tangent_circles(c1, c2)
        cstart = HC_CC_Circle(q2, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        cend = HC_CC_Circle(q3, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        q4 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = (
            cstart.hc_turn_length(q1)
            + configuration_distance(q2, q3)
            + cend.hc_turn_length(q4)
        )
        return length, cstart, cend, q1, q2, q3, q4

    def TeST_path(self, c1, c2):
        q2, q3 = self.TeST_tangent_circles(c1, c2)
        cstart = HC_CC_Circle(q2, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        cend = HC_CC_Circle(q3, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        q4 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)
        length = (
            cstart.hc_turn_length(q1)
            + configuration_distance(q2, q3)
            + cend.hc_turn_length(q4)
        )
        return length, cstart, cend, q1, q2, q3, q4

    def TST_path(self, c1, c2):
        if self.TiST_exists(c1, c2):
            return self.TiST_path(c1, c2)
        if self.TeST_exists(c1, c2):
            return self.TeST_path(c1, c2)
        return _INF, None, None, None, None, None, None

    # ----- TTT ---------------------------------------------------------------
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
        delta_y = math.sqrt(r ** 2 - delta_x ** 2)

        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular,
                            self._parent.hc_cc_circle_param_)
        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt2 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular,
                            self._parent.hc_cc_circle_param_)

        qa = self.TT_tangent_circles(c1, tgt1)
        qb = self.TT_tangent_circles(tgt1, c2)
        qc = self.TT_tangent_circles(c1, tgt2)
        qd = self.TT_tangent_circles(tgt2, c2)
        return qa, qb, qc, qd

    def TTT_path(self, c1, c2):
        qa, qb, qc, qd = self.TTT_tangent_circles(c1, c2)

        start1 = HC_CC_Circle(qa, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True,
                               self._parent.hc_cc_circle_param_)
        end1 = HC_CC_Circle(qb, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        start2 = HC_CC_Circle(qc, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        middle2 = HC_CC_Circle(qc, not c1.left, c1.forward, True,
                               self._parent.hc_cc_circle_param_)
        end2 = HC_CC_Circle(qd, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)

        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)

        length1 = (
            start1.hc_turn_length(q1)
            + middle1.cc_turn_length(qb)
            + end1.hc_turn_length(q3)
        )
        length2 = (
            start2.hc_turn_length(q1)
            + middle2.cc_turn_length(qd)
            + end2.hc_turn_length(q3)
        )
        if length1 < length2:
            return length1, start1, end1, q1, qb, q3, middle1
        else:
            return length2, start2, end2, q1, qd, q3, middle2

    # ----- TTTT --------------------------------------------------------------
    def TTTT_exists(self, c1, c2):
        if c1.left == c2.left:
            return False
        if c1.forward == c2.forward:
            return False
        return self.distance <= 6 * self._parent.radius_

    def TTTT_tangent_circles(self, c1, c2):
        theta = self.angle
        r1 = 2 * self._parent.radius_
        if self.distance < r1:
            delta_x = (self.distance + r1) / 2
            delta_y = math.sqrt(r1 ** 2 - delta_x ** 2)
        else:
            delta_x = (self.distance - r1) / 2
            delta_y = math.sqrt(r1 ** 2 - delta_x ** 2)

        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y)
        tgt1 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular,
                            self._parent.hc_cc_circle_param_)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y)
        tgt2 = HC_CC_Circle(x, y, not c2.left, not c2.forward, c2.regular,
                            self._parent.hc_cc_circle_param_)

        x, y = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y)
        tgt3 = HC_CC_Circle(x, y, not c1.left, c1.forward, c1.regular,
                            self._parent.hc_cc_circle_param_)
        x, y = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y)
        tgt4 = HC_CC_Circle(x, y, not c2.left, not c2.forward, c2.regular,
                            self._parent.hc_cc_circle_param_)

        qa = self.TT_tangent_circles(c1, tgt1)
        qb = self.TT_tangent_circles(tgt1, tgt2)
        qc = self.TT_tangent_circles(tgt2, c2)

        qd = self.TT_tangent_circles(c1, tgt3)
        qe = self.TT_tangent_circles(tgt3, tgt4)
        qf = self.TT_tangent_circles(tgt4, c2)
        return qa, qb, qc, qd, qe, qf

    def TTTT_path(self, c1, c2):
        qa, qb, qc, qd, qe, qf = self.TTTT_tangent_circles(c1, c2)

        start1 = HC_CC_Circle(qa, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        middle1 = HC_CC_Circle(qa, not c1.left, c1.forward, True,
                               self._parent.hc_cc_circle_param_)
        middle2 = HC_CC_Circle(qc, not c2.left, c2.forward, True,
                               self._parent.hc_cc_circle_param_)
        end1 = HC_CC_Circle(qc, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)
        start2 = HC_CC_Circle(qd, c1.left, not c1.forward, True,
                              self._parent.hc_cc_circle_param_)
        middle3 = HC_CC_Circle(qd, not c1.left, c1.forward, True,
                               self._parent.hc_cc_circle_param_)
        middle4 = HC_CC_Circle(qf, not c2.left, c2.forward, True,
                               self._parent.hc_cc_circle_param_)
        end2 = HC_CC_Circle(qf, c2.left, not c2.forward, True,
                            self._parent.hc_cc_circle_param_)

        q1 = Configuration(c1.start.x, c1.start.y, c1.start.theta, c1.kappa)
        q3 = Configuration(c2.start.x, c2.start.y, c2.start.theta, c2.kappa)

        length1 = (
            start1.hc_turn_length(q1)
            + middle1.cc_turn_length(qb)
            + middle2.cc_turn_length(qb)
            + end1.hc_turn_length(q3)
        )
        length2 = (
            start2.hc_turn_length(q1)
            + middle3.cc_turn_length(qe)
            + middle4.cc_turn_length(qe)
            + end2.hc_turn_length(q3)
        )
        if length1 < length2:
            return length1, start1, end1, q1, qb, q3, middle1, middle2
        else:
            return length2, start2, end2, q1, qe, q3, middle3, middle4


# =============================================================================
# CCpmpm_Dubins_State_Space
# =============================================================================
class CCpmpm_Dubins_State_Space(HC_CC_StateSpace):
    """CC Dubins state space with ±max curvature at both start and end."""

    _CONTROL_TABLE = {
        cc_dubins_path_type.E: [('empty',)],
        cc_dubins_path_type.T: [('rs', 'cstart', 'end', True)],
        cc_dubins_path_type.TT: [('hc', 'cstart', 'qi1', False), ('hc', 'cend', 'qi3', True)],
        cc_dubins_path_type.TST: [('hc', 'cstart', 'qi1', False), ('straight', 'qi2', 'qi3'), ('hc', 'cend', 'qi4', True)],
        cc_dubins_path_type.TTT: [('hc', 'cstart', 'qi1', False), ('cc', 'ci1', 'qi2', True), ('hc', 'cend', 'qi3', True)],
        cc_dubins_path_type.TTTT: [('hc', 'cstart', 'qi1', False), ('cc', 'ci1', 'qi2', True), ('cc', 'ci2', 'qi2', False), ('hc', 'cend', 'qi3', True)],
    }

    def __init__(self, kappa, sigma, discretization=0.1, forwards=True):
        super().__init__(kappa, sigma, discretization)
        self.forwards_ = forwards
        self._ccpmpm_dubins = _CCpmpm_Dubins(self)
        self.rs_circle_param_ = HC_CC_Circle_Param()
        self.rs_circle_param_.set_param(kappa, _INF, 1.0 / kappa, 0.0, 0.0, 1.0, 0.0)
        self.radius_ = self.hc_cc_circle_param_.radius
        self.mu_ = self.hc_cc_circle_param_.mu
        self.sin_mu_ = self.hc_cc_circle_param_.sin_mu
        self.cos_mu_ = self.hc_cc_circle_param_.cos_mu

    def ccpmpm_circles_dubins_path(self, c1, c2):
        length = double_array_init(nb_cc_dubins_paths, _INF)
        qi1 = [None] * nb_cc_dubins_paths
        qi2 = [None] * nb_cc_dubins_paths
        qi3 = [None] * nb_cc_dubins_paths
        qi4 = [None] * nb_cc_dubins_paths
        cstart = [None] * nb_cc_dubins_paths
        ci1 = [None] * nb_cc_dubins_paths
        ci2 = [None] * nb_cc_dubins_paths
        cend = [None] * nb_cc_dubins_paths

        self._ccpmpm_dubins.distance = center_distance(c1, c2)
        self._ccpmpm_dubins.angle = math.atan2(c2.yc - c1.yc, c2.xc - c1.xc)

        E = cc_dubins_path_type.E
        T = cc_dubins_path_type.T
        TT = cc_dubins_path_type.TT
        TST = cc_dubins_path_type.TST
        TTT = cc_dubins_path_type.TTT
        TTTT = cc_dubins_path_type.TTTT

        # case E
        if configuration_equal(c1.start, c2.start):
            length[E] = 0
        # case T
        elif configuration_on_hc_cc_circle(c1, c2.start):
            cstart[T] = HC_CC_Circle(c1.start, c1.left, c1.forward, True,
                                     self.rs_circle_param_)
            length[T] = cstart[T].rs_turn_length(c2.start)
        else:
            # case TT
            if self._ccpmpm_dubins.TT_exists(c1, c2):
                l, cs, ce, q1, q2, q3 = self._ccpmpm_dubins.TT_path(c1, c2)
                length[TT] = l
                cstart[TT] = cs
                cend[TT] = ce
                qi1[TT] = q1
                qi2[TT] = q2
                qi3[TT] = q3
            # case TST
            if self._ccpmpm_dubins.TST_exists(c1, c2):
                l, cs, ce, q1, q2, q3, q4 = self._ccpmpm_dubins.TST_path(c1, c2)
                length[TST] = l
                cstart[TST] = cs
                cend[TST] = ce
                qi1[TST] = q1
                qi2[TST] = q2
                qi3[TST] = q3
                qi4[TST] = q4
            # case TTT
            if self._ccpmpm_dubins.TTT_exists(c1, c2):
                l, cs, ce, q1, q2, q3, mid = self._ccpmpm_dubins.TTT_path(c1, c2)
                length[TTT] = l
                cstart[TTT] = cs
                cend[TTT] = ce
                qi1[TTT] = q1
                qi2[TTT] = q2
                qi3[TTT] = q3
                ci1[TTT] = mid
            # case TTTT
            if self._ccpmpm_dubins.TTTT_exists(c1, c2):
                l, cs, ce, q1, q2, q3, m1, m2 = self._ccpmpm_dubins.TTTT_path(c1, c2)
                length[TTTT] = l
                cstart[TTTT] = cs
                cend[TTTT] = ce
                qi1[TTTT] = q1
                qi2[TTTT] = q2
                qi3[TTTT] = q3
                ci1[TTTT] = m1
                ci2[TTTT] = m2

        best = array_index_min(length)
        best_path = cc_dubins_path_type(best)
        return CC_Dubins_Path(
            c1.start, c2.start, best_path, self.kappa_, self.sigma_,
            qi1[best], qi2[best], qi3[best], qi4[best],
            cstart[best], cend[best], ci1[best], ci2[best],
            length[best],
        )

    def _find_path(self, state1, state2):
        start1 = Configuration(state1.x, state1.y, state1.theta, self.kappa_)
        start2 = Configuration(state1.x, state1.y, state1.theta, -self.kappa_)
        end1 = Configuration(state2.x, state2.y, state2.theta, self.kappa_)
        end2 = Configuration(state2.x, state2.y, state2.theta, -self.kappa_)

        if self.forwards_:
            start_circles = [
                HC_CC_Circle(start1, True, True, True, self.rs_circle_param_),
                HC_CC_Circle(start2, False, True, True, self.rs_circle_param_),
            ]
            end_circles = [
                HC_CC_Circle(end1, True, False, True, self.rs_circle_param_),
                HC_CC_Circle(end2, False, False, True, self.rs_circle_param_),
            ]
        else:
            start_circles = [
                HC_CC_Circle(start1, True, False, True, self.rs_circle_param_),
                HC_CC_Circle(start2, False, False, True, self.rs_circle_param_),
            ]
            end_circles = [
                HC_CC_Circle(end1, True, True, True, self.rs_circle_param_),
                HC_CC_Circle(end2, False, True, True, self.rs_circle_param_),
            ]

        paths = [None] * 4
        lg = [_INF] * 4
        for i in range(2):
            # skip circle at the beginning for curvature continuity
            if i == 0 and state1.kappa < 0:
                continue
            elif i == 1 and state1.kappa > 0:
                continue
            for j in range(2):
                # skip circle at the end for curvature continuity
                if j == 0 and state2.kappa < 0:
                    continue
                elif j == 1 and state2.kappa > 0:
                    continue
                idx = 2 * i + j
                paths[idx] = self.ccpmpm_circles_dubins_path(
                    start_circles[i], end_circles[j]
                )
                if paths[idx] is not None:
                    lg[idx] = paths[idx].length

        best = array_index_min(lg)
        return paths[best]


# =============================================================================
# CC_Dubins_State_Space (General wrapper)
# =============================================================================
class CC_Dubins_State_Space(HC_CC_StateSpace):
    """General CC Dubins state space that selects the shortest path
    among all four curvature boundary variants."""

    def __init__(self, kappa, sigma, discretization=0.1, forwards=True):
        super().__init__(kappa, sigma, discretization)
        self.forwards_ = forwards
        self.cc00_dubins_state_space_ = CC00_Dubins_State_Space(
            kappa, sigma, discretization, forwards
        )
        self.cc0pm_dubins_state_space_ = CC0pm_Dubins_State_Space(
            kappa, sigma, discretization, forwards
        )
        self.ccpm0_dubins_state_space_ = CCpm0_Dubins_State_Space(
            kappa, sigma, discretization, forwards
        )
        self.ccpmpm_dubins_state_space_ = CCpmpm_Dubins_State_Space(
            kappa, sigma, discretization, forwards
        )

    def predict_state(self, state, forwards):
        """Predict states forwards or backwards to zero and max curvature."""
        states_controls = []

        # no prediction required
        if (abs(state.kappa) < get_epsilon()) or (
            (self.kappa_ - abs(state.kappa)) < get_epsilon()
        ):
            sc = (
                State(
                    x=state.x, y=state.y, theta=state.theta,
                    kappa=state.kappa, d=state.d,
                ),
                Control(delta_s=0.0, kappa=state.kappa, sigma=0.0),
            )
            states_controls.append(sc)
            return states_controls

        sgn_kappa = sgn(state.kappa)

        if forwards:
            c1 = Control(
                delta_s=(self.kappa_ - sgn_kappa * state.kappa) / self.sigma_,
                kappa=state.kappa,
                sigma=sgn_kappa * self.sigma_,
            )
            c2 = Control(
                delta_s=sgn_kappa * state.kappa / self.sigma_,
                kappa=state.kappa,
                sigma=-sgn_kappa * self.sigma_,
            )
        else:
            c1 = Control(
                delta_s=-(self.kappa_ - sgn_kappa * state.kappa) / self.sigma_,
                kappa=state.kappa,
                sigma=sgn_kappa * self.sigma_,
            )
            c2 = Control(
                delta_s=-sgn_kappa * state.kappa / self.sigma_,
                kappa=state.kappa,
                sigma=-sgn_kappa * self.sigma_,
            )

        for ctrl in (c1, c2):
            d = sgn(ctrl.delta_s)
            abs_delta_s = abs(ctrl.delta_s)
            xf, yf, thf, kf = end_of_clothoid(
                state.x, state.y, state.theta, state.kappa,
                ctrl.sigma, d, abs_delta_s,
            )
            predicted = State(x=xf, y=yf, theta=thf, kappa=kf)
            states_controls.append((predicted, ctrl))

        return states_controls

    def _select_sub_space(self, start_state, end_state):
        """Select appropriate sub-state-space based on endpoint curvatures."""
        eps = get_epsilon()
        if abs(start_state.kappa) < eps:
            return self.cc0pm_dubins_state_space_ if abs(end_state.kappa) >= eps else self.cc00_dubins_state_space_
        else:
            return self.ccpmpm_dubins_state_space_ if abs(end_state.kappa) >= eps else self.ccpm0_dubins_state_space_

    def get_distance(self, state1, state2):
        start_scs = self.predict_state(state1, self.forwards_)
        end_scs = self.predict_state(state2, not self.forwards_)
        distances = []

        for start_state, start_control in start_scs:
            for end_state, end_control in end_scs:
                if state_equal(start_state, end_state):
                    ctrl = subtract_control(start_control, end_control)
                    distances.append(abs(ctrl.delta_s))
                else:
                    distance = self._select_sub_space(start_state, end_state).get_distance(start_state, end_state)
                    if abs(start_control.delta_s) > get_epsilon():
                        distance += abs(start_control.delta_s)
                    if abs(end_control.delta_s) > get_epsilon():
                        distance += abs(end_control.delta_s)
                    distances.append(distance)

        return min(distances)

    def get_controls(self, state1, state2):
        start_scs = self.predict_state(state1, self.forwards_)
        end_scs = self.predict_state(state2, not self.forwards_)
        controls_distance_pairs = []

        for start_state, start_control in start_scs:
            for end_state, end_control in end_scs:
                cc_dubins_controls = []
                if state_equal(start_state, end_state):
                    ctrl = subtract_control(start_control, end_control)
                    cc_dubins_controls.append(ctrl)
                else:
                    cc_dubins_controls = self._select_sub_space(start_state, end_state).get_controls(start_state, end_state)
                    # adjust controls by initial and final control
                    if abs(start_control.delta_s) > get_epsilon():
                        cc_dubins_controls.insert(0, start_control)
                    if abs(end_control.delta_s) > get_epsilon():
                        ec = Control(
                            delta_s=end_control.delta_s,
                            kappa=end_control.kappa,
                            sigma=end_control.sigma,
                        )
                        reverse_control(ec)
                        cc_dubins_controls.append(ec)

                distance = sum(abs(c.delta_s) for c in cc_dubins_controls)
                controls_distance_pairs.append((cc_dubins_controls, distance))

        controls_distance_pairs.sort(key=lambda p: p[1])
        return controls_distance_pairs[0][0]

    def get_controls_reverse(self, state1, state2):
        """Return controls for pure reverse driving by computing the forward
        path from state2 to state1 and reversing it."""
        cc_dubins_controls = self.get_controls(state2, state1)
        cc_dubins_controls.reverse()
        for ctrl in cc_dubins_controls:
            reverse_control(ctrl)
        return cc_dubins_controls

    def get_all_controls(self, state1, state2):
        """Return all possible control sequences (forward and reverse)."""
        return [
            self.get_controls(state1, state2),
            self.get_controls_reverse(state1, state2),
        ]
