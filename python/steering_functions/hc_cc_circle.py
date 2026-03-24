"""Port of hc_cc_core/hc_cc_circle.hpp and hc_cc_circle.cpp."""

import math

from steering_functions.utilities import (
    PI,
    TWO_PI,
    HALF_PI,
    get_epsilon,
    point_distance,
    twopify,
    fresnel,
    global_frame_change,
)
from steering_functions.configuration import Configuration


class HC_CC_Circle_Param:
    """Parameters shared by all HC/CC circles."""

    def __init__(self):
        self.kappa = 0.0
        self.kappa_inv = 0.0
        self.sigma = 0.0
        self.radius = 0.0
        self.mu = 0.0
        self.sin_mu = 0.0
        self.cos_mu = 0.0
        self.delta_min = 0.0

    def set_param(self, kappa, sigma, radius, mu, sin_mu, cos_mu, delta_min):
        """Set circle parameters."""
        self.kappa = kappa
        self.kappa_inv = 1.0 / kappa
        self.sigma = sigma
        self.radius = radius
        self.mu = mu
        self.sin_mu = sin_mu
        self.cos_mu = cos_mu
        self.delta_min = delta_min


class HC_CC_Circle(HC_CC_Circle_Param):
    """HC/CC circle defined by a start configuration and turning parameters."""

    def __init__(self, start_or_xc, y_or_left=None, left_or_forward=None,
                 forward_or_regular=None, regular_or_param=None, param=None):
        super().__init__()
        if isinstance(start_or_xc, Configuration):
            # Constructor 1: from Configuration
            self._init_from_configuration(
                start_or_xc, y_or_left, left_or_forward,
                forward_or_regular, regular_or_param)
        else:
            # Constructor 2: from center coordinates
            self._init_from_center(
                start_or_xc, y_or_left, left_or_forward,
                forward_or_regular, regular_or_param, param)

    def _init_from_configuration(self, _start, _left, _forward, _regular, _param):
        self.start = _start
        self.left = _left
        self.forward = _forward
        self.regular = _regular
        delta_x = _param.radius * _param.sin_mu
        delta_y = _param.radius * _param.cos_mu
        if _left:
            self.kappa = _param.kappa
            self.kappa_inv = _param.kappa_inv
            self.sigma = _param.sigma
            if _forward:
                self.xc, self.yc = global_frame_change(
                    _start.x, _start.y, _start.theta, delta_x, delta_y)
            else:
                self.xc, self.yc = global_frame_change(
                    _start.x, _start.y, _start.theta, -delta_x, delta_y)
        else:
            self.kappa = -_param.kappa
            self.kappa_inv = -_param.kappa_inv
            self.sigma = -_param.sigma
            if _forward:
                self.xc, self.yc = global_frame_change(
                    _start.x, _start.y, _start.theta, delta_x, -delta_y)
            else:
                self.xc, self.yc = global_frame_change(
                    _start.x, _start.y, _start.theta, -delta_x, -delta_y)
        self.radius = _param.radius
        self.mu = _param.mu
        self.sin_mu = _param.sin_mu
        self.cos_mu = _param.cos_mu
        self.delta_min = _param.delta_min

    def _init_from_center(self, _xc, _yc, _left, _forward, _regular, _param):
        self.start = Configuration(0, 0, 0, 0)
        self.left = _left
        self.forward = _forward
        self.regular = _regular
        if _left:
            self.kappa = _param.kappa
            self.kappa_inv = _param.kappa_inv
            self.sigma = _param.sigma
        else:
            self.kappa = -_param.kappa
            self.kappa_inv = -_param.kappa_inv
            self.sigma = -_param.sigma
        self.xc = _xc
        self.yc = _yc
        self.radius = _param.radius
        self.mu = _param.mu
        self.sin_mu = _param.sin_mu
        self.cos_mu = _param.cos_mu
        self.delta_min = _param.delta_min

    def deflection(self, q):
        """Angle between start configuration of circle and configuration *q*."""
        alpha_c = self.start.theta
        alpha_q = q.theta
        if (self.left and self.forward) or (not self.left and not self.forward):
            return twopify(alpha_q - alpha_c)
        else:
            return twopify(alpha_c - alpha_q)

    def D1(self, alpha):
        """D1 helper for elementary path evaluation."""
        s = math.sqrt(2.0 * alpha / PI)
        fresnel_s, fresnel_c = fresnel(s)
        return math.cos(alpha) * fresnel_c + math.sin(alpha) * fresnel_s

    def rs_circular_deflection(self, delta):
        """Circular deflection of a rs-turn."""
        if self.regular:
            return delta
        else:
            return delta if delta <= PI else delta - TWO_PI

    def rs_turn_length(self, q):
        """Length of a rs-turn."""
        assert (abs(abs(self.kappa) - abs(q.kappa)) < get_epsilon()
                and abs(abs(self.sigma) - float('inf')) < get_epsilon())
        delta = self.deflection(q)
        return abs(self.kappa_inv * self.rs_circular_deflection(delta))

    def hc_circular_deflection(self, delta):
        """Circular deflection of a hc-turn."""
        delta_min_twopified = twopify(self.delta_min)
        if self.regular:
            if delta < delta_min_twopified:
                return TWO_PI + delta - delta_min_twopified
            else:
                return delta - delta_min_twopified
        else:
            if delta < delta_min_twopified:
                delta_arc1 = delta - delta_min_twopified
                delta_arc2 = delta_arc1 + TWO_PI
            else:
                delta_arc1 = delta - delta_min_twopified
                delta_arc2 = delta_arc1 - TWO_PI
            return delta_arc1 if abs(delta_arc1) < abs(delta_arc2) else delta_arc2

    def hc_turn_length(self, q):
        """Length of a hc-turn."""
        assert abs(abs(self.kappa) - abs(q.kappa)) < get_epsilon()
        delta = self.deflection(q)
        return (abs(self.kappa / self.sigma)
                + abs(self.kappa_inv * self.hc_circular_deflection(delta)))

    def cc_elementary_sharpness(self, q, delta):
        """Sharpness of an elementary path. Returns (success, sigma0)."""
        distance = point_distance(self.start.x, self.start.y, q.x, q.y)
        if delta < 4.5948 and distance > get_epsilon():
            sigma0 = 4.0 * PI * pow(self.D1(0.5 * delta), 2) / pow(distance, 2)
            if not self.left:
                sigma0 = -sigma0
            return True, sigma0
        return False, 0.0

    def cc_circular_deflection(self, delta):
        """Circular deflection of a cc-turn."""
        two_delta_min_twopified = twopify(2.0 * self.delta_min)
        if self.regular:
            if delta < two_delta_min_twopified:
                return TWO_PI + delta - two_delta_min_twopified
            else:
                return delta - two_delta_min_twopified
        else:
            if delta < two_delta_min_twopified:
                delta_arc1 = delta - two_delta_min_twopified
                delta_arc2 = delta_arc1 + TWO_PI
            else:
                delta_arc1 = delta - two_delta_min_twopified
                delta_arc2 = delta_arc1 - TWO_PI
            return delta_arc1 if abs(delta_arc1) < abs(delta_arc2) else delta_arc2

    def cc_turn_length(self, q):
        """Length of a cc-turn."""
        assert abs(q.kappa) < get_epsilon()
        delta = self.deflection(q)
        # delta = 0
        if delta < get_epsilon():
            return 2.0 * self.radius * self.sin_mu
        # 0 < delta < 2 * delta_min
        length_min = abs(self.kappa / self.sigma)
        length_default = (
            2.0 * length_min
            + abs(self.kappa_inv * self.cc_circular_deflection(delta))
        )
        if delta < 2.0 * self.delta_min:
            success, sigma0 = self.cc_elementary_sharpness(q, delta)
            if success:
                length_elementary = 2.0 * math.sqrt(delta / abs(sigma0))
                return min(length_elementary, length_default)
            else:
                return length_default
        # delta >= 2 * delta_min
        else:
            return length_default


def center_distance(c1, c2):
    """Cartesian distance between the centres of two circles."""
    return math.sqrt((c2.xc - c1.xc) ** 2 + (c2.yc - c1.yc) ** 2)


def configuration_on_hc_cc_circle(c, q):
    """Check whether configuration *q* lies on circle *c*."""
    distance = point_distance(c.xc, c.yc, q.x, q.y)
    if abs(distance - c.radius) > get_epsilon():
        return False
    angle = math.atan2(q.y - c.yc, q.x - c.xc)
    if c.left and c.forward:
        angle = angle + HALF_PI - c.mu
    if c.left and not c.forward:
        angle = angle + HALF_PI + c.mu
    if not c.left and c.forward:
        angle = angle - HALF_PI + c.mu
    if not c.left and not c.forward:
        angle = angle - HALF_PI - c.mu
    angle = twopify(angle)
    return abs(q.theta - angle) < get_epsilon()
