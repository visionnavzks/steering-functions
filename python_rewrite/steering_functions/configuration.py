"""Port of hc_cc_core/configuration.hpp and configuration.cpp."""

import math

from steering_functions.utilities import twopify, point_distance, get_epsilon


class Configuration:
    """Kinematic configuration (position, orientation, curvature)."""

    def __init__(self, x=0.0, y=0.0, theta=0.0, kappa=0.0):
        self.x = x
        self.y = y
        self.theta = twopify(theta)
        self.kappa = kappa

    def __repr__(self):
        return (
            f"Configuration({self.x}, {self.y}, {self.theta}, {self.kappa})"
        )


def configuration_distance(q1, q2):
    """Cartesian distance between two configurations."""
    return point_distance(q1.x, q1.y, q2.x, q2.y)


def configuration_aligned(q1, q2):
    """Check whether two configurations are aligned."""
    if abs(q2.theta - q1.theta) > get_epsilon():
        return False
    angle = twopify(math.atan2(q2.y - q1.y, q2.x - q1.x))
    return abs(angle - q1.theta) <= get_epsilon()


def configuration_equal(q1, q2):
    """Check whether two configurations are equal."""
    if abs(q2.theta - q1.theta) > get_epsilon():
        return False
    if configuration_distance(q1, q2) > get_epsilon():
        return False
    return True
