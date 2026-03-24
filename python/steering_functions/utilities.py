# Copyright (c) 2017 - for information on the respective copyright owner
# see the NOTICE file and/or the repository
#     https://github.com/hbanzhaf/steering_functions.git
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Derived from Continuous Curvature (CC) Steer.
# Copyright (c) 2016, Thierry Fraichard and Institut national de
# recherche en informatique et en automatique (Inria), licensed under
# the BSD license.

"""Math utility functions for the steering functions library.

Includes angle normalisation, Fresnel integrals (via Chebyshev approximation),
and geometric helpers for clothoid, circular-arc, and straight-line segments.
"""

import math
from typing import Tuple

PI = math.pi
HALF_PI = math.pi / 2.0
TWO_PI = 2.0 * math.pi
SQRT_PI = math.sqrt(math.pi)
SQRT_PI_INV = 1.0 / math.sqrt(math.pi)
SQRT_TWO_PI_INV = 1.0 / math.sqrt(2.0 * math.pi)

EPSILON = 1e-4


def get_epsilon() -> float:
    """Return numerical epsilon used for near-zero comparisons."""
    return EPSILON


def sgn(x: float) -> float:
    """Return the sign of *x* as ±1.0 (positive for zero)."""
    return -1.0 if x < 0 else 1.0


def point_distance(x1: float, y1: float, x2: float, y2: float) -> float:
    """Cartesian distance between two points."""
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def polar(x: float, y: float) -> Tuple[float, float]:
    """Return (r, theta) polar coordinates of the Cartesian point (x, y)."""
    r = math.sqrt(x * x + y * y)
    theta = math.atan2(y, x)
    return r, theta


def twopify(alpha: float) -> float:
    """Map *alpha* (radians) into [0, 2π)."""
    return alpha - TWO_PI * math.floor(alpha / TWO_PI)


def pify(alpha: float) -> float:
    """Map *alpha* (radians) into (−π, π]."""
    v = math.fmod(alpha, TWO_PI)
    if v < -PI:
        v += TWO_PI
    elif v > PI:
        v -= TWO_PI
    return v


# ---------------------------------------------------------------------------
# Fresnel integrals via Chebyshev polynomial approximation
# Source: Y. L. Luke, "Mathematical Functions and their Approximations,"
# Academic Press, 1975, pp. 139-142.
# ---------------------------------------------------------------------------

_CHEBEV_A = [
    0.76435138664186000189,
    -0.43135547547660179313,
    0.43288199979726653054,
    -0.26973310338387111029,
    0.08416045320876935378,
    -0.01546524484461381958,
    0.00187855423439822018,
    -0.00016264977618887547,
    0.00001057397656383260,
    -0.00000053609339889243,
    0.00000002181658454933,
    -0.00000000072901621186,
    0.00000000002037332548,
    -0.00000000000048344033,
    0.00000000000000986533,
    -0.00000000000000017502,
    0.00000000000000000272,
    -0.00000000000000000004,
]

_CHEBEV_B = [
    0.63041404314570539241,
    -0.42344511405705333544,
    0.37617172643343656625,
    -0.16249489154509567415,
    0.03822255778633008694,
    -0.00564563477132190899,
    0.00057454951976897367,
    -0.00004287071532102004,
    0.00000245120749923299,
    -0.00000011098841840868,
    0.00000000408249731696,
    -0.00000000012449830219,
    0.00000000000320048425,
    -0.00000000000007032416,
    0.00000000000000133638,
    -0.00000000000000002219,
    0.00000000000000000032,
]

_CHEBEV_E = [
    0.97462779093296822410,
    -0.02424701873969321371,
    0.00103400906842977317,
    -0.00008052450246908016,
    0.00000905962481966582,
    -0.00000131016996757743,
    0.00000022770820391497,
    -0.00000004558623552026,
    0.00000001021567537083,
    -0.00000000251114508133,
    0.00000000066704761275,
    -0.00000000018931512852,
    0.00000000005689898935,
    -0.00000000001798219359,
    0.00000000000594162963,
    -0.00000000000204285065,
    0.00000000000072797580,
    -0.00000000000026797428,
    0.00000000000010160694,
    -0.00000000000003958559,
    0.00000000000001581262,
    -0.00000000000000646411,
    0.00000000000000269981,
    -0.00000000000000115038,
    0.00000000000000049942,
    -0.00000000000000022064,
    0.00000000000000009910,
    -0.00000000000000004520,
    0.00000000000000002092,
    -0.00000000000000000982,
    0.00000000000000000467,
    -0.00000000000000000225,
    0.00000000000000000110,
    -0.00000000000000000054,
    0.00000000000000000027,
    -0.00000000000000000014,
    0.00000000000000000007,
    -0.00000000000000000004,
    0.00000000000000000002,
    -0.00000000000000000001,
    0.00000000000000000001,
]

_CHEBEV_F = [
    0.99461545179407928910,
    -0.00524276766084297210,
    0.00013325864229883909,
    -0.00000770856452642713,
    0.00000070848077032045,
    -0.00000008812517411602,
    0.00000001359784717148,
    -0.00000000246858295747,
    0.00000000050925789921,
    -0.00000000011653400634,
    0.00000000002906578309,
    -0.00000000000779847361,
    0.00000000000222802542,
    -0.00000000000067239338,
    0.00000000000021296411,
    -0.00000000000007041482,
    0.00000000000002419805,
    -0.00000000000000861080,
    0.00000000000000316287,
    -0.00000000000000119596,
    0.00000000000000046444,
    -0.00000000000000018485,
    0.00000000000000007527,
    -0.00000000000000003131,
    0.00000000000000001328,
    -0.00000000000000000574,
    0.00000000000000000252,
    -0.00000000000000000113,
    0.00000000000000000051,
    -0.00000000000000000024,
    0.00000000000000000011,
    -0.00000000000000000005,
    0.00000000000000000002,
    -0.00000000000000000001,
    0.00000000000000000001,
]


def _fresnel_0_8(x: float) -> Tuple[float, float]:
    """Fresnel integrals for x ∈ [0, 8] using Chebyshev polynomials."""
    quarter_x = 0.25 * x
    arg = 0.03125 * x * x - 1.0
    T0 = 1.0
    T1 = 0.125 * x
    T2 = arg
    T3 = quarter_x * T2 - T1
    A = _CHEBEV_A[0] * T0 + _CHEBEV_A[1] * T2
    B = _CHEBEV_B[0] * T1 + _CHEBEV_B[1] * T3
    T2n_m4 = T0
    T2n_m2 = T2
    T2n_m1 = T3
    for n in range(2, 17):
        T2n = 2.0 * arg * T2n_m2 - T2n_m4
        T2n_p1 = quarter_x * T2n - T2n_m1
        A += _CHEBEV_A[n] * T2n
        B += _CHEBEV_B[n] * T2n_p1
        T2n_m4 = T2n_m2
        T2n_m2 = T2n
        T2n_m1 = T2n_p1
    T34 = 2.0 * arg * T2n_m2 - T2n_m4
    A += _CHEBEV_A[17] * T34

    sqrt_x = math.sqrt(x)
    C_f = SQRT_TWO_PI_INV * sqrt_x * A
    S_f = SQRT_TWO_PI_INV * sqrt_x * B
    return S_f, C_f


def _fresnel_8_inf(x: float) -> Tuple[float, float]:
    """Fresnel integrals for x > 8 using asymptotic Chebyshev expansion."""
    arg = 128.0 / (x * x) - 1.0
    T0 = 1.0
    T2 = arg
    E = _CHEBEV_E[0] * T0 + _CHEBEV_E[1] * T2
    F = _CHEBEV_F[0] * T0 + _CHEBEV_F[1] * T2
    T2n_m4 = T0
    T2n_m2 = T2
    for n in range(2, 35):
        T2n = 2.0 * arg * T2n_m2 - T2n_m4
        E += _CHEBEV_E[n] * T2n
        F += _CHEBEV_F[n] * T2n
        T2n_m4 = T2n_m2
        T2n_m2 = T2n
    for n in range(35, 41):
        T2n = 2.0 * arg * T2n_m2 - T2n_m4
        E += _CHEBEV_E[n] * T2n
        T2n_m4 = T2n_m2
        T2n_m2 = T2n

    sin_x = math.sin(x)
    cos_x = math.cos(x)
    sqrt_x = math.sqrt(x)
    C_f = 0.5 - SQRT_TWO_PI_INV * (E * cos_x / (2 * x) - F * sin_x) / sqrt_x
    S_f = 0.5 - SQRT_TWO_PI_INV * (E * sin_x / (2 * x) + F * cos_x) / sqrt_x
    return S_f, C_f


def fresnel(s: float) -> Tuple[float, float]:
    """Compute Fresnel integrals S and C at *s*.

    Returns (S_f, C_f) where::

        S_f = ∫₀ˢ sin(π/2 · u²) du
        C_f = ∫₀ˢ cos(π/2 · u²) du

    Uses Chebyshev polynomial approximations for both small and large arguments.
    """
    x = HALF_PI * s * s
    if x <= 8.0:
        S_f, C_f = _fresnel_0_8(x)
    else:
        S_f, C_f = _fresnel_8_inf(x)
    if s < 0:
        S_f = -S_f
        C_f = -C_f
    return S_f, C_f


# ---------------------------------------------------------------------------
# Geometric primitives
# ---------------------------------------------------------------------------


def end_of_clothoid(
    x_i: float,
    y_i: float,
    theta_i: float,
    kappa_i: float,
    sigma: float,
    direction: float,
    length: float,
) -> Tuple[float, float, float, float]:
    """Propagate along a clothoid segment.

    Returns (x_f, y_f, theta_f, kappa_f).
    """
    sgn_sigma = sgn(sigma)
    abs_sigma = abs(sigma)
    sqrt_sigma_inv = 1.0 / math.sqrt(abs_sigma)
    k1 = theta_i - 0.5 * direction * kappa_i * kappa_i / sigma
    k2 = SQRT_PI_INV * sqrt_sigma_inv * (abs_sigma * length + sgn_sigma * kappa_i)
    k3 = SQRT_PI_INV * sqrt_sigma_inv * sgn_sigma * kappa_i
    cos_k1 = math.cos(k1)
    sin_k1 = math.sin(k1)
    fresnel_s_k2, fresnel_c_k2 = fresnel(k2)
    fresnel_s_k3, fresnel_c_k3 = fresnel(k3)
    x_f = x_i + SQRT_PI * sqrt_sigma_inv * (
        direction * cos_k1 * (fresnel_c_k2 - fresnel_c_k3)
        - sgn_sigma * sin_k1 * (fresnel_s_k2 - fresnel_s_k3)
    )
    y_f = y_i + SQRT_PI * sqrt_sigma_inv * (
        direction * sin_k1 * (fresnel_c_k2 - fresnel_c_k3)
        + sgn_sigma * cos_k1 * (fresnel_s_k2 - fresnel_s_k3)
    )
    theta_f = pify(
        theta_i
        + kappa_i * direction * length
        + 0.5 * sigma * direction * length * length
    )
    kappa_f = kappa_i + sigma * length
    return x_f, y_f, theta_f, kappa_f


def end_of_circular_arc(
    x_i: float,
    y_i: float,
    theta_i: float,
    kappa: float,
    direction: float,
    length: float,
) -> Tuple[float, float, float]:
    """Propagate along a circular arc segment.

    Returns (x_f, y_f, theta_f).
    """
    x_f = x_i + (1.0 / kappa) * (
        -math.sin(theta_i) + math.sin(theta_i + direction * length * kappa)
    )
    y_f = y_i + (1.0 / kappa) * (
        math.cos(theta_i) - math.cos(theta_i + direction * length * kappa)
    )
    theta_f = pify(theta_i + kappa * direction * length)
    return x_f, y_f, theta_f


def end_of_straight_line(
    x_i: float,
    y_i: float,
    theta: float,
    direction: float,
    length: float,
) -> Tuple[float, float]:
    """Propagate along a straight line segment.

    Returns (x_f, y_f).
    """
    x_f = x_i + direction * length * math.cos(theta)
    y_f = y_i + direction * length * math.sin(theta)
    return x_f, y_f


def global_frame_change(
    x: float, y: float, theta: float, local_x: float, local_y: float
) -> Tuple[float, float]:
    """Transform a point from a local frame to the global frame."""
    sin_th = math.sin(theta)
    cos_th = math.cos(theta)
    global_x = local_x * cos_th - local_y * sin_th + x
    global_y = local_x * sin_th + local_y * cos_th + y
    return global_x, global_y


def local_frame_change(
    x: float, y: float, theta: float, global_x: float, global_y: float
) -> Tuple[float, float]:
    """Transform a point from the global frame to a local frame."""
    sin_th = math.sin(theta)
    cos_th = math.cos(theta)
    local_x = (global_x - x) * cos_th + (global_y - y) * sin_th
    local_y = -(global_x - x) * sin_th + (global_y - y) * cos_th
    return local_x, local_y
