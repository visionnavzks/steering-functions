"""Port of hc_cc_core/paths.hpp and paths.cpp."""

import math
from enum import IntEnum

from steering_functions.utilities import (
    get_epsilon,
    point_distance,
    twopify,
    sgn,
)
from steering_functions.state import Control
from steering_functions.configuration import Configuration
from steering_functions.hc_cc_circle import HC_CC_Circle


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------

class cc_dubins_path_type(IntEnum):
    E = 0
    S = 1
    T = 2
    TT = 3
    TST = 4
    TTT = 5
    TTTT = 6


nb_cc_dubins_paths = 7


class hc_cc_rs_path_type(IntEnum):
    E = 0
    S = 1
    T = 2
    TT = 3
    TcT = 4
    TcTcT = 5
    TcTT = 6
    TTcT = 7
    TST = 8
    TSTcT = 9
    TcTST = 10
    TcTSTcT = 11
    TTcTT = 12
    TcTTcT = 13
    TTT = 14
    TcST = 15
    TScT = 16
    TcScT = 17


nb_hc_cc_rs_paths = 18


# ---------------------------------------------------------------------------
# Path classes
# ---------------------------------------------------------------------------

class Path:
    """Base path with start/end configurations and path parameters."""

    def __init__(self, start, end, kappa, sigma, length):
        self.start = start
        self.end = end
        self.kappa = kappa
        self.sigma = sigma
        self.length = length


class CC_Dubins_Path(Path):
    """CC-Dubins path with intermediate configurations and circles."""

    def __init__(self, start, end, path_type, kappa, sigma,
                 qi1, qi2, qi3, qi4,
                 cstart, cend, ci1, ci2, length):
        super().__init__(start, end, kappa, sigma, length)
        self.type = path_type
        self.qi1 = qi1
        self.qi2 = qi2
        self.qi3 = qi3
        self.qi4 = qi4
        self.cstart = cstart
        self.cend = cend
        self.ci1 = ci1
        self.ci2 = ci2


class HC_CC_RS_Path(Path):
    """HC/CC Reeds-Shepp path with intermediate configurations and circles."""

    def __init__(self, start, end, path_type, kappa, sigma,
                 qi1, qi2, qi3, qi4,
                 cstart, cend, ci1, ci2, length):
        super().__init__(start, end, kappa, sigma, length)
        self.type = path_type
        self.qi1 = qi1
        self.qi2 = qi2
        self.qi3 = qi3
        self.qi4 = qi4
        self.cstart = cstart
        self.cend = cend
        self.ci1 = ci1
        self.ci2 = ci2


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def state_equal(state1, state2):
    """Check whether two states are equal."""
    if abs(state2.kappa - state1.kappa) > get_epsilon():
        return False
    if abs(twopify(state2.theta) - twopify(state1.theta)) > get_epsilon():
        return False
    if point_distance(state1.x, state1.y, state2.x, state2.y) > get_epsilon():
        return False
    return True


def reverse_control(control):
    """Reverse a control in-place."""
    control.delta_s = -control.delta_s
    control.kappa = control.kappa + abs(control.delta_s) * control.sigma
    control.sigma = -control.sigma


def subtract_control(control1, control2):
    """Subtract *control2* from *control1* and return the result."""
    assert (sgn(control1.delta_s) * control1.sigma
            == sgn(control2.delta_s) * control2.sigma)
    control = Control()
    control.delta_s = control1.delta_s - control2.delta_s
    control.kappa = control1.kappa
    control.sigma = control1.sigma
    return control


def empty_controls(controls):
    """Append a zero-input control."""
    control = Control()
    control.delta_s = 0.0
    control.kappa = 0.0
    control.sigma = 0.0
    controls.append(control)


def straight_controls(q1, q2, controls):
    """Append controls for a straight line from *q1* to *q2*."""
    length = point_distance(q1.x, q1.y, q2.x, q2.y)
    dot_product = (math.cos(q1.theta) * (q2.x - q1.x)
                   + math.sin(q1.theta) * (q2.y - q1.y))
    d = sgn(dot_product)
    control = Control()
    control.delta_s = d * length
    control.kappa = 0.0
    control.sigma = 0.0
    controls.append(control)


def _direction(forward, order):
    """Return +1 or -1 depending on driving direction and segment order."""
    if (forward and order) or (not forward and not order):
        return 1
    else:
        return -1


def rs_turn_controls(c, q, order, controls):
    """Append controls for a rs-turn."""
    import sys as _sys
    assert (abs(abs(c.kappa) - abs(q.kappa)) < get_epsilon()
            and abs(abs(c.sigma) - _sys.float_info.max) < get_epsilon())
    delta = c.deflection(q)
    length_arc = abs(c.kappa_inv) * c.rs_circular_deflection(delta)
    d = _direction(c.forward, order)

    arc = Control()
    arc.delta_s = d * length_arc
    arc.kappa = c.kappa
    arc.sigma = 0.0
    controls.append(arc)


def hc_turn_controls(c, q, order, controls):
    """Append controls for a hc-turn."""
    assert abs(abs(c.kappa) - abs(q.kappa)) < get_epsilon()
    delta = c.deflection(q)
    length_min = abs(c.kappa / c.sigma)
    length_arc = abs(c.kappa_inv) * c.hc_circular_deflection(delta)
    d = _direction(c.forward, order)

    if order:
        clothoid = Control()
        clothoid.delta_s = d * length_min
        clothoid.kappa = 0.0
        clothoid.sigma = c.sigma
        controls.append(clothoid)

    arc = Control()
    arc.delta_s = d * length_arc
    arc.kappa = c.kappa
    arc.sigma = 0.0
    controls.append(arc)

    if not order:
        clothoid = Control()
        clothoid.delta_s = d * length_min
        clothoid.kappa = c.kappa
        clothoid.sigma = -c.sigma
        controls.append(clothoid)


def cc_elementary_controls(c, q, delta, order, controls):
    """Append controls for an elementary path if one exists. Returns True on success."""
    success, sigma0 = c.cc_elementary_sharpness(q, delta)
    if success:
        length = math.sqrt(delta / abs(sigma0))
        d = _direction(c.forward, order)

        clothoid1 = Control()
        clothoid1.delta_s = d * length
        clothoid1.kappa = 0.0
        clothoid1.sigma = sigma0
        controls.append(clothoid1)

        clothoid2 = Control()
        clothoid2.delta_s = d * length
        clothoid2.kappa = sigma0 * length
        clothoid2.sigma = -sigma0
        controls.append(clothoid2)
        return True
    return False


def cc_default_controls(c, q, delta, order, controls):
    """Append controls for a default cc-turn (two clothoids and a circular arc)."""
    length_min = abs(c.kappa / c.sigma)
    length_arc = abs(c.kappa_inv) * c.cc_circular_deflection(delta)
    d = _direction(c.forward, order)

    clothoid1 = Control()
    clothoid1.delta_s = d * length_min
    clothoid1.kappa = 0.0
    clothoid1.sigma = c.sigma
    controls.append(clothoid1)

    arc = Control()
    arc.delta_s = d * length_arc
    arc.kappa = c.kappa
    arc.sigma = 0.0
    controls.append(arc)

    clothoid2 = Control()
    clothoid2.delta_s = d * length_min
    clothoid2.kappa = c.kappa
    clothoid2.sigma = -c.sigma
    controls.append(clothoid2)


def cc_turn_controls(c, q, order, controls):
    """Append controls for a cc-turn."""
    assert abs(q.kappa) < get_epsilon()
    delta = c.deflection(q)
    # delta ≈ 0
    if delta < get_epsilon():
        if order:
            straight_controls(c.start, q, controls)
        else:
            straight_controls(q, c.start, controls)
        return
    # 0 < delta < 2 * delta_min
    if delta < 2.0 * c.delta_min:
        controls_elementary = []
        if cc_elementary_controls(c, q, delta, order, controls_elementary):
            controls_default = []
            cc_default_controls(c, q, delta, order, controls_default)
            length_elementary = sum(abs(ctrl.delta_s) for ctrl in controls_elementary)
            length_default = sum(abs(ctrl.delta_s) for ctrl in controls_default)
            if length_elementary < length_default:
                controls.extend(controls_elementary)
            else:
                controls.extend(controls_default)
            return
        else:
            cc_default_controls(c, q, delta, order, controls)
            return
    # delta >= 2 * delta_min
    else:
        cc_default_controls(c, q, delta, order, controls)
        return


def _build_controls(path, control_table):
    """Build controls from a path using a declarative control table.

    *control_table* maps ``path.type`` → list of control steps.
    Each step is a tuple:
      - ``('empty',)``
      - ``('straight', from_attr, to_attr)``
      - ``('cc', circle_attr, config_attr, order)``
      - ``('hc', circle_attr, config_attr, order)``
      - ``('rs', circle_attr, config_attr, order)``

    Attribute names are resolved via ``getattr(path, name)``.
    """
    steps = control_table.get(path.type)
    if steps is None:
        return []
    controls = []
    _TURN_FNS = {'cc': cc_turn_controls, 'hc': hc_turn_controls, 'rs': rs_turn_controls}
    for step in steps:
        op = step[0]
        if op == 'empty':
            empty_controls(controls)
        elif op == 'straight':
            straight_controls(getattr(path, step[1]), getattr(path, step[2]), controls)
        else:
            _TURN_FNS[op](getattr(path, step[1]), getattr(path, step[2]), step[3], controls)
    return controls
