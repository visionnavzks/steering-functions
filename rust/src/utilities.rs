// Copyright (c) 2017 - for information on the respective copyright
// owner see the NOTICE file and/or the repository
//
//     https://github.com/hbanzhaf/steering_functions.git
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// This source code is derived from Continuous Curvature (CC) Steer.
// Copyright (c) 2016, Thierry Fraichard and Institut national de
// recherche en informatique et en automatique (Inria), licensed under
// the BSD license, cf. 3rd-party-licenses.txt file in the root
// directory of this source tree.

pub const PI: f64 = 3.1415926535897932384;
pub const HALF_PI: f64 = 1.5707963267948966192;
pub const TWO_PI: f64 = 6.2831853071795864770;
pub const SQRT_PI: f64 = 1.7724538509055160273;
pub const SQRT_PI_INV: f64 = 0.56418958354775628695;
pub const SQRT_TWO_PI_INV: f64 = 0.39894228040143267794;

/// Tolerance used for geometric comparisons.
pub const EPSILON: f64 = 1e-4;

// ---------------------------------------------------------------------------
// Chebyshev coefficient arrays for Fresnel integral approximation
// ---------------------------------------------------------------------------

#[rustfmt::skip]
const CHEBEV_A: [f64; 18] = [
    0.76435138664186000189, -0.43135547547660179313,  0.43288199979726653054,
   -0.26973310338387111029,  0.08416045320876935378, -0.01546524484461381958,
    0.00187855423439822018, -0.00016264977618887547,  0.00001057397656383260,
   -0.00000053609339889243,  0.00000002181658454933, -0.00000000072901621186,
    0.00000000002037332548, -0.00000000000048344033,  0.00000000000000986533,
   -0.00000000000000017502,  0.00000000000000000272, -0.00000000000000000004,
];

#[rustfmt::skip]
const CHEBEV_B: [f64; 17] = [
    0.63041404314570539241, -0.42344511405705333544,  0.37617172643343656625,
   -0.16249489154509567415,  0.03822255778633008694, -0.00564563477132190899,
    0.00057454951976897367, -0.00004287071532102004,  0.00000245120749923299,
   -0.00000011098841840868,  0.00000000408249731696, -0.00000000012449830219,
    0.00000000000320048425, -0.00000000000007032416,  0.00000000000000133638,
   -0.00000000000000002219,  0.00000000000000000032,
];

#[rustfmt::skip]
const CHEBEV_E: [f64; 41] = [
    0.97462779093296822410, -0.02424701873969321371,  0.00103400906842977317,
   -0.00008052450246908016,  0.00000905962481966582, -0.00000131016996757743,
    0.00000022770820391497, -0.00000004558623552026,  0.00000001021567537083,
   -0.00000000251114508133,  0.00000000066704761275, -0.00000000018931512852,
    0.00000000005689898935, -0.00000000001798219359,  0.00000000000594162963,
   -0.00000000000204285065,  0.00000000000072797580, -0.00000000000026797428,
    0.00000000000010160694, -0.00000000000003958559,  0.00000000000001581262,
   -0.00000000000000646411,  0.00000000000000269981, -0.00000000000000115038,
    0.00000000000000049942, -0.00000000000000022064,  0.00000000000000009910,
   -0.00000000000000004520,  0.00000000000000002092, -0.00000000000000000982,
    0.00000000000000000467, -0.00000000000000000225,  0.00000000000000000110,
   -0.00000000000000000054,  0.00000000000000000027, -0.00000000000000000014,
    0.00000000000000000007, -0.00000000000000000004,  0.00000000000000000002,
   -0.00000000000000000001,  0.00000000000000000001,
];

#[rustfmt::skip]
const CHEBEV_F: [f64; 35] = [
    0.99461545179407928910, -0.00524276766084297210,  0.00013325864229883909,
   -0.00000770856452642713,  0.00000070848077032045, -0.00000008812517411602,
    0.00000001359784717148, -0.00000000246858295747,  0.00000000050925789921,
   -0.00000000011653400634,  0.00000000002906578309, -0.00000000000779847361,
    0.00000000000222802542, -0.00000000000067239338,  0.00000000000021296411,
   -0.00000000000007041482,  0.00000000000002419805, -0.00000000000000861080,
    0.00000000000000316287, -0.00000000000000119596,  0.00000000000000046444,
   -0.00000000000000018485,  0.00000000000000007527, -0.00000000000000003131,
    0.00000000000000001328, -0.00000000000000000574,  0.00000000000000000252,
   -0.00000000000000000113,  0.00000000000000000051, -0.00000000000000000024,
    0.00000000000000000011, -0.00000000000000000005,  0.00000000000000000002,
   -0.00000000000000000001,  0.00000000000000000001,
];

// ---------------------------------------------------------------------------
// Basic math utilities
// ---------------------------------------------------------------------------

/// Return the epsilon tolerance.
#[inline]
pub fn get_epsilon() -> f64 {
    EPSILON
}

/// Return the sign of `x` (-1.0 if negative, 1.0 otherwise).
#[inline]
pub fn sgn(x: f64) -> f64 {
    if x < 0.0 {
        -1.0
    } else {
        1.0
    }
}

/// Cartesian distance between two points.
#[inline]
pub fn point_distance(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt()
}

/// Compute polar coordinates `(r, theta)` from Cartesian `(x, y)`.
#[inline]
pub fn polar(x: f64, y: f64) -> (f64, f64) {
    let r = (x * x + y * y).sqrt();
    let theta = y.atan2(x);
    (r, theta)
}

/// Convert an arbitrary angle (radians) into `[0, 2π)`.
#[inline]
pub fn twopify(alpha: f64) -> f64 {
    alpha - TWO_PI * (alpha / TWO_PI).floor()
}

/// Convert an arbitrary angle (radians) into `[-π, π)`.
#[inline]
pub fn pify(alpha: f64) -> f64 {
    let v = alpha % TWO_PI;
    if v < -PI {
        v + TWO_PI
    } else if v > PI {
        v - TWO_PI
    } else {
        v
    }
}

// ---------------------------------------------------------------------------
// Fresnel integrals
// ---------------------------------------------------------------------------

/// Fresnel integral approximation for `x ∈ [0, 8]` using Chebyshev polynomials.
/// Returns `(s_f, c_f)`.
fn fresnel_0_8(x: f64) -> (f64, f64) {
    let quarter_x = 0.25 * x;
    let arg = 0.03125 * x * x - 1.0;

    let t0: f64 = 1.0;
    let t1 = 0.125 * x;
    let t2 = arg;
    let t3 = quarter_x * t2 - t1;

    let mut a = CHEBEV_A[0] * t0 + CHEBEV_A[1] * t2;
    let mut b = CHEBEV_B[0] * t1 + CHEBEV_B[1] * t3;

    let mut t2n_m4 = t0;
    let mut t2n_m2 = t2;
    let mut t2n_m1 = t3;

    for n in 2..CHEBEV_B.len() {
        let t2n = 2.0 * arg * t2n_m2 - t2n_m4;
        let t2n_p1 = quarter_x * t2n - t2n_m1;
        a += CHEBEV_A[n] * t2n;
        b += CHEBEV_B[n] * t2n_p1;
        t2n_m4 = t2n_m2;
        t2n_m2 = t2n;
        t2n_m1 = t2n_p1;
    }

    let t_last = 2.0 * arg * t2n_m2 - t2n_m4;
    a += CHEBEV_A[CHEBEV_A.len() - 1] * t_last;

    let sqrt_x = x.sqrt();
    let c_f = SQRT_TWO_PI_INV * sqrt_x * a;
    let s_f = SQRT_TWO_PI_INV * sqrt_x * b;
    (s_f, c_f)
}

/// Fresnel integral approximation for `x ∈ (8, ∞)` using Chebyshev polynomials.
/// Returns `(s_f, c_f)`.
fn fresnel_8_inf(x: f64) -> (f64, f64) {
    let arg = 128.0 / (x * x) - 1.0;

    let t0: f64 = 1.0;
    let t2 = arg;

    let mut e = CHEBEV_E[0] * t0 + CHEBEV_E[1] * t2;
    let mut f = CHEBEV_F[0] * t0 + CHEBEV_F[1] * t2;

    let mut t2n_m4 = t0;
    let mut t2n_m2 = t2;

    for n in 2..CHEBEV_F.len() {
        let t2n = 2.0 * arg * t2n_m2 - t2n_m4;
        e += CHEBEV_E[n] * t2n;
        f += CHEBEV_F[n] * t2n;
        t2n_m4 = t2n_m2;
        t2n_m2 = t2n;
    }

    for n in CHEBEV_F.len()..CHEBEV_E.len() {
        let t2n = 2.0 * arg * t2n_m2 - t2n_m4;
        e += CHEBEV_E[n] * t2n;
        t2n_m4 = t2n_m2;
        t2n_m2 = t2n;
    }

    let sin_x = x.sin();
    let cos_x = x.cos();
    let sqrt_x = x.sqrt();

    let c_f = 0.5 - SQRT_TWO_PI_INV * (e * cos_x / (2.0 * x) - f * sin_x) / sqrt_x;
    let s_f = 0.5 - SQRT_TWO_PI_INV * (e * sin_x / (2.0 * x) + f * cos_x) / sqrt_x;
    (s_f, c_f)
}

/// Compute the Fresnel integrals for a given arc length `s`.
///
/// Returns `(s_f, c_f)` where:
/// - `s_f = ∫₀ˢ sin(π/2 · u²) du`
/// - `c_f = ∫₀ˢ cos(π/2 · u²) du`
pub fn fresnel(s: f64) -> (f64, f64) {
    let x = HALF_PI * s * s;
    let (mut s_f, mut c_f) = if x <= 8.0 {
        fresnel_0_8(x)
    } else {
        fresnel_8_inf(x)
    };
    if s < 0.0 {
        s_f = -s_f;
        c_f = -c_f;
    }
    (s_f, c_f)
}

// ---------------------------------------------------------------------------
// Geometry: end-of-segment computations
// ---------------------------------------------------------------------------

/// Compute the end configuration on a clothoid.
///
/// # Arguments
/// * `x_i`, `y_i`, `theta_i`, `kappa_i` – initial configuration
/// * `sigma` – sharpness of clothoid
/// * `direction` – driving direction (`-1.0` or `1.0`)
/// * `length` – length of clothoid (positive)
///
/// # Returns
/// `(x_f, y_f, theta_f, kappa_f)` – final configuration on clothoid.
pub fn end_of_clothoid(
    x_i: f64,
    y_i: f64,
    theta_i: f64,
    kappa_i: f64,
    sigma: f64,
    direction: f64,
    length: f64,
) -> (f64, f64, f64, f64) {
    let sgn_sigma = sgn(sigma);
    let abs_sigma = sigma.abs();
    let sqrt_sigma_inv = 1.0 / abs_sigma.sqrt();

    let k1 = theta_i - 0.5 * direction * kappa_i * kappa_i / sigma;
    let k2 = SQRT_PI_INV * sqrt_sigma_inv * (abs_sigma * length + sgn_sigma * kappa_i);
    let k3 = SQRT_PI_INV * sqrt_sigma_inv * sgn_sigma * kappa_i;

    let cos_k1 = k1.cos();
    let sin_k1 = k1.sin();

    let (fresnel_s_k2, fresnel_c_k2) = fresnel(k2);
    let (fresnel_s_k3, fresnel_c_k3) = fresnel(k3);

    let x_f = x_i
        + SQRT_PI
            * sqrt_sigma_inv
            * (direction * cos_k1 * (fresnel_c_k2 - fresnel_c_k3)
                - sgn_sigma * sin_k1 * (fresnel_s_k2 - fresnel_s_k3));

    let y_f = y_i
        + SQRT_PI
            * sqrt_sigma_inv
            * (direction * sin_k1 * (fresnel_c_k2 - fresnel_c_k3)
                + sgn_sigma * cos_k1 * (fresnel_s_k2 - fresnel_s_k3));

    let theta_f = pify(
        theta_i + kappa_i * direction * length + 0.5 * sigma * direction * length * length,
    );

    let kappa_f = kappa_i + sigma * length;

    (x_f, y_f, theta_f, kappa_f)
}

/// Compute the end configuration on a circular arc.
///
/// # Arguments
/// * `x_i`, `y_i`, `theta_i` – initial configuration
/// * `kappa` – curvature of circular arc
/// * `direction` – driving direction (`-1.0` or `1.0`)
/// * `length` – arc length (positive)
///
/// # Returns
/// `(x_f, y_f, theta_f)` – final configuration on circular arc.
pub fn end_of_circular_arc(
    x_i: f64,
    y_i: f64,
    theta_i: f64,
    kappa: f64,
    direction: f64,
    length: f64,
) -> (f64, f64, f64) {
    let x_f =
        x_i + (1.0 / kappa) * (-(theta_i.sin()) + (theta_i + direction * length * kappa).sin());
    let y_f =
        y_i + (1.0 / kappa) * (theta_i.cos() - (theta_i + direction * length * kappa).cos());
    let theta_f = pify(theta_i + kappa * direction * length);
    (x_f, y_f, theta_f)
}

/// Compute the end point on a straight line.
///
/// # Arguments
/// * `x_i`, `y_i` – initial position
/// * `theta` – angle of straight line
/// * `direction` – driving direction (`-1.0` or `1.0`)
/// * `length` – line length (positive)
///
/// # Returns
/// `(x_f, y_f)` – final position on straight line.
pub fn end_of_straight_line(
    x_i: f64,
    y_i: f64,
    theta: f64,
    direction: f64,
    length: f64,
) -> (f64, f64) {
    let x_f = x_i + direction * length * theta.cos();
    let y_f = y_i + direction * length * theta.sin();
    (x_f, y_f)
}

// ---------------------------------------------------------------------------
// Frame transforms
// ---------------------------------------------------------------------------

/// Transform `(local_x, local_y)` from a local coordinate system to the global frame.
///
/// # Returns
/// `(global_x, global_y)`
pub fn global_frame_change(
    x: f64,
    y: f64,
    theta: f64,
    local_x: f64,
    local_y: f64,
) -> (f64, f64) {
    let sin_th = theta.sin();
    let cos_th = theta.cos();
    let global_x = local_x * cos_th - local_y * sin_th + x;
    let global_y = local_x * sin_th + local_y * cos_th + y;
    (global_x, global_y)
}

/// Transform `(global_x, global_y)` from the global coordinate system to a local frame.
///
/// # Returns
/// `(local_x, local_y)`
pub fn local_frame_change(
    x: f64,
    y: f64,
    theta: f64,
    global_x: f64,
    global_y: f64,
) -> (f64, f64) {
    let sin_th = theta.sin();
    let cos_th = theta.cos();
    let local_x = (global_x - x) * cos_th + (global_y - y) * sin_th;
    let local_y = -(global_x - x) * sin_th + (global_y - y) * cos_th;
    (local_x, local_y)
}

// ---------------------------------------------------------------------------
// Array utilities
// ---------------------------------------------------------------------------

/// Find the index of the minimum value in a slice.
pub fn array_index_min(array: &[f64]) -> usize {
    assert!(!array.is_empty(), "array must not be empty");
    let mut min = array[0];
    let mut index_min = 0;
    for (i, &val) in array.iter().enumerate().skip(1) {
        if val < min {
            index_min = i;
            min = val;
        }
    }
    index_min
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS_FRESNEL: f64 = 1e-14;

    #[test]
    fn test_sgn() {
        assert_eq!(sgn(-5.0), -1.0);
        assert_eq!(sgn(0.0), 1.0);
        assert_eq!(sgn(3.0), 1.0);
    }

    #[test]
    fn test_point_distance() {
        let d = point_distance(0.0, 0.0, 3.0, 4.0);
        assert!((d - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_polar() {
        let (r, theta) = polar(1.0, 0.0);
        assert!((r - 1.0).abs() < 1e-12);
        assert!(theta.abs() < 1e-12);
    }

    #[test]
    fn test_twopify() {
        let a = twopify(-0.5);
        assert!(a >= 0.0 && a < TWO_PI);
        let b = twopify(7.0);
        assert!(b >= 0.0 && b < TWO_PI);
    }

    #[test]
    fn test_pify() {
        let v = pify(4.0);
        assert!(v >= -PI && v <= PI);
        let v2 = pify(-4.0);
        assert!(v2 >= -PI && v2 <= PI);
    }

    // Fresnel tests compared with scipy.special.fresnel values
    #[test]
    fn test_fresnel_negative_large() {
        let (s, c) = fresnel(-40.0);
        assert!((s - -0.4920422537902731).abs() < EPS_FRESNEL);
        assert!((c - -0.4999984168574456).abs() < EPS_FRESNEL);
    }

    #[test]
    fn test_fresnel_negative_medium() {
        let (s, c) = fresnel(-2.0);
        assert!((s - -0.34341567836369824).abs() < EPS_FRESNEL);
        assert!((c - -0.4882534060753408).abs() < EPS_FRESNEL);
    }

    #[test]
    fn test_fresnel_zero() {
        let (s, c) = fresnel(0.0);
        assert!(s.abs() < EPS_FRESNEL);
        assert!(c.abs() < EPS_FRESNEL);
    }

    #[test]
    fn test_fresnel_positive_small() {
        let (s, c) = fresnel(1.0);
        assert!((s - 0.4382591473903548).abs() < EPS_FRESNEL);
        assert!((c - 0.7798934003768228).abs() < EPS_FRESNEL);
    }

    #[test]
    fn test_fresnel_positive_large() {
        let (s, c) = fresnel(40.0);
        assert!((s - 0.4920422537902731).abs() < EPS_FRESNEL);
        assert!((c - 0.4999984168574456).abs() < EPS_FRESNEL);
    }

    #[test]
    fn test_end_of_straight_line() {
        let (x, y) = end_of_straight_line(0.0, 0.0, 0.0, 1.0, 5.0);
        assert!((x - 5.0).abs() < 1e-12);
        assert!(y.abs() < 1e-12);
    }

    #[test]
    fn test_global_local_frame_roundtrip() {
        let (gx, gy) = global_frame_change(1.0, 2.0, 0.5, 3.0, 4.0);
        let (lx, ly) = local_frame_change(1.0, 2.0, 0.5, gx, gy);
        assert!((lx - 3.0).abs() < 1e-12);
        assert!((ly - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_array_index_min() {
        assert_eq!(array_index_min(&[3.0, 1.0, 2.0]), 1);
        assert_eq!(array_index_min(&[5.0]), 0);
        assert_eq!(array_index_min(&[2.0, 2.0, 1.0, 3.0]), 2);
    }

    #[test]
    fn test_get_epsilon() {
        assert_eq!(get_epsilon(), EPSILON);
    }
}
