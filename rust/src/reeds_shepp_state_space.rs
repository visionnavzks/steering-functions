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

use crate::base_state_space::StateSpace;
use crate::state::{Control, State};
use crate::utilities::{pify, polar, PI};

const RS_EPS: f64 = 1e-6;
const RS_ZERO: f64 = 10.0 * f64::EPSILON;

// ---------------------------------------------------------------------------
// Reeds-Shepp path types and structures
// ---------------------------------------------------------------------------

/// Type of a single Reeds-Shepp path segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReedsSheppPathSegmentType {
    /// No operation (padding for shorter paths).
    Nop = 0,
    /// Turn left (positive curvature).
    Left = 1,
    /// Drive straight (zero curvature).
    Straight = 2,
    /// Turn right (negative curvature).
    Right = 3,
}

/// The 18 canonical Reeds-Shepp path families.
const REEDS_SHEPP_PATH_TYPE: [[ReedsSheppPathSegmentType; 5]; 18] = {
    use ReedsSheppPathSegmentType::*;
    [
        [Left, Right, Left, Nop, Nop],           // 0
        [Right, Left, Right, Nop, Nop],           // 1
        [Left, Right, Left, Right, Nop],          // 2
        [Right, Left, Right, Left, Nop],          // 3
        [Left, Right, Straight, Left, Nop],       // 4
        [Right, Left, Straight, Right, Nop],      // 5
        [Left, Straight, Right, Left, Nop],       // 6
        [Right, Straight, Left, Right, Nop],      // 7
        [Left, Right, Straight, Right, Nop],      // 8
        [Right, Left, Straight, Left, Nop],       // 9
        [Right, Straight, Right, Left, Nop],      // 10
        [Left, Straight, Left, Right, Nop],       // 11
        [Left, Straight, Right, Nop, Nop],        // 12
        [Right, Straight, Left, Nop, Nop],        // 13
        [Left, Straight, Left, Nop, Nop],         // 14
        [Right, Straight, Right, Nop, Nop],       // 15
        [Left, Right, Straight, Left, Right],     // 16
        [Right, Left, Straight, Right, Left],     // 17
    ]
};

/// A Reeds-Shepp path consisting of up to five segments.
///
/// Each segment is a left turn, right turn, straight line, or no-op.
/// Lengths are in *normalised* units (multiply by `1/kappa` to get actual
/// arc length). Negative lengths indicate reverse driving.
#[derive(Debug, Clone, Copy)]
pub struct ReedsSheppPath {
    /// Segment types (e.g. Left–Right–Left–Nop–Nop).
    pub segment_types: [ReedsSheppPathSegmentType; 5],
    /// Length of each segment in normalised units (may be negative for reverse).
    pub lengths: [f64; 5],
    /// Sum of absolute segment lengths.
    pub total_length: f64,
}

impl ReedsSheppPath {
    /// Create a new Reeds-Shepp path with up to five segment lengths.
    pub fn new(
        segment_types: [ReedsSheppPathSegmentType; 5],
        t: f64,
        u: f64,
        v: f64,
        w: f64,
        x: f64,
    ) -> Self {
        Self {
            segment_types,
            lengths: [t, u, v, w, x],
            total_length: t.abs() + u.abs() + v.abs() + w.abs() + x.abs(),
        }
    }

    /// Total normalised length of the path (sum of absolute segment lengths).
    pub fn length(&self) -> f64 {
        self.total_length
    }

    /// Get the normalised length of segment `index`.
    pub fn get_length(&self, index: usize) -> f64 {
        self.lengths[index]
    }
}

impl Default for ReedsSheppPath {
    /// Sentinel path with maximum length (represents "no valid path found").
    fn default() -> Self {
        Self {
            segment_types: REEDS_SHEPP_PATH_TYPE[0],
            lengths: [f64::MAX, 0.0, 0.0, 0.0, 0.0],
            total_length: f64::MAX,
        }
    }
}

// ---------------------------------------------------------------------------
// Helper functions for Reeds-Shepp path computation
// ---------------------------------------------------------------------------

#[inline]
fn tau_omega(u: f64, v: f64, xi: f64, eta: f64, phi: f64) -> (f64, f64) {
    let delta = pify(u - v);
    let a = u.sin() - delta.sin();
    let b = u.cos() - delta.cos() - 1.0;
    let t1 = (eta * a - xi * b).atan2(xi * a + eta * b);
    let t2 = 2.0 * (delta.cos() - v.cos() - u.cos()) + 3.0;
    let tau = if t2 < 0.0 {
        pify(t1 + PI)
    } else {
        pify(t1)
    };
    let omega = pify(tau - u + v - phi);
    (tau, omega)
}

// formula 8.1
#[inline]
fn lp_sp_lp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    let (u, t) = polar(x - phi.sin(), y - 1.0 + phi.cos());
    if t >= -RS_ZERO {
        let v = pify(phi - t);
        if v >= -RS_ZERO {
            return Some((t, u, v));
        }
    }
    None
}

// formula 8.2
#[inline]
fn lp_sp_rp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    let (u1, t1) = polar(x + phi.sin(), y - 1.0 - phi.cos());
    let u1_sq = u1 * u1;
    if u1_sq >= 4.0 {
        let u = (u1_sq - 4.0).sqrt();
        let theta = (2.0f64).atan2(u);
        let t = pify(t1 + theta);
        let v = pify(t - phi);
        if t >= -RS_ZERO && v >= -RS_ZERO {
            return Some((t, u, v));
        }
    }
    None
}

// formula 8.3/8.4
#[inline]
fn lp_rm_l(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    let xi = x - phi.sin();
    let eta = y - 1.0 + phi.cos();
    let (u1, theta) = polar(xi, eta);
    if u1 <= 4.0 {
        let u = -2.0 * (0.25 * u1).asin();
        let t = pify(theta + 0.5 * u + PI);
        let v = pify(phi - t + u);
        if t >= -RS_ZERO && u <= RS_ZERO {
            return Some((t, u, v));
        }
    }
    None
}

// formula 8.7
#[inline]
fn lp_rup_lum_rm(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    let xi = x + phi.sin();
    let eta = y - 1.0 - phi.cos();
    let rho = 0.25 * (2.0 + (xi * xi + eta * eta).sqrt());
    if rho <= 1.0 {
        let u = rho.acos();
        let (t, v) = tau_omega(u, -u, xi, eta, phi);
        if t >= -RS_ZERO && v <= RS_ZERO {
            return Some((t, u, v));
        }
    }
    None
}

// formula 8.8
#[inline]
fn lp_rum_lum_rp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    let xi = x + phi.sin();
    let eta = y - 1.0 - phi.cos();
    let rho = (20.0 - xi * xi - eta * eta) / 16.0;
    if (0.0..=1.0).contains(&rho) {
        let u = -rho.acos();
        if u >= -0.5 * PI {
            let (t, v) = tau_omega(u, u, xi, eta, phi);
            if t >= -RS_ZERO && v >= -RS_ZERO {
                return Some((t, u, v));
            }
        }
    }
    None
}

// formula 8.9
#[inline]
fn lp_rm_sm_lm(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    let xi = x - phi.sin();
    let eta = y - 1.0 + phi.cos();
    let (rho, theta) = polar(xi, eta);
    if rho >= 2.0 {
        let r = (rho * rho - 4.0).sqrt();
        let u = 2.0 - r;
        let t = pify(theta + r.atan2(-2.0));
        let v = pify(phi - 0.5 * PI - t);
        if t >= -RS_ZERO && u <= RS_ZERO && v <= RS_ZERO {
            return Some((t, u, v));
        }
    }
    None
}

// formula 8.10
#[inline]
fn lp_rm_sm_rm(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    let xi = x + phi.sin();
    let eta = y - 1.0 - phi.cos();
    let (rho, theta) = polar(-eta, xi);
    if rho >= 2.0 {
        let t = theta;
        let u = 2.0 - rho;
        let v = pify(t + 0.5 * PI - phi);
        if t >= -RS_ZERO && u <= RS_ZERO && v <= RS_ZERO {
            return Some((t, u, v));
        }
    }
    None
}

// formula 8.11
#[inline]
fn lp_rm_s_lm_rp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    let xi = x + phi.sin();
    let eta = y - 1.0 - phi.cos();
    let (rho, _theta) = polar(xi, eta);
    if rho >= 2.0 {
        let u = 4.0 - (rho * rho - 4.0).sqrt();
        if u <= RS_ZERO {
            let t = pify(((4.0 - u) * xi - 2.0 * eta).atan2(-2.0 * xi + (u - 4.0) * eta));
            let v = pify(t - phi);
            if t >= -RS_ZERO && v >= -RS_ZERO {
                return Some((t, u, v));
            }
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Path family collectors
// ---------------------------------------------------------------------------

/// Collect all CSC paths.
fn csc(x: f64, y: f64, phi: f64, paths: &mut Vec<ReedsSheppPath>) {
    // LpSpLp - type 14
    if let Some((t, u, v)) = lp_sp_lp(x, y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[14], t, u, v, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_sp_lp(-x, y, -phi) {
        // timeflip
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[14], -t, -u, -v, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_sp_lp(x, -y, -phi) {
        // reflect
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[15], t, u, v, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_sp_lp(-x, -y, phi) {
        // timeflip + reflect
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[15], -t, -u, -v, 0.0, 0.0));
    }

    // LpSpRp - type 12
    if let Some((t, u, v)) = lp_sp_rp(x, y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[12], t, u, v, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_sp_rp(-x, y, -phi) {
        // timeflip
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[12], -t, -u, -v, 0.0, 0.0));
    }
    // LpSpRp reflect - type 13
    if let Some((t, u, v)) = lp_sp_rp(x, -y, -phi) {
        // reflect
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[13], t, u, v, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_sp_rp(-x, -y, phi) {
        // timeflip + reflect
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[13], -t, -u, -v, 0.0, 0.0));
    }
}

/// Collect all CCC paths.
fn ccc(x: f64, y: f64, phi: f64, paths: &mut Vec<ReedsSheppPath>) {
    // type 0
    if let Some((t, u, v)) = lp_rm_l(x, y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[0], t, u, v, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_l(-x, y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[0], -t, -u, -v, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_l(x, -y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[1], t, u, v, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_l(-x, -y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[1], -t, -u, -v, 0.0, 0.0));
    }

    // backwards
    let xb = x * phi.cos() + y * phi.sin();
    let yb = x * phi.sin() - y * phi.cos();
    if let Some((t, u, v)) = lp_rm_l(xb, yb, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[0], v, u, t, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_l(-xb, yb, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[0], -v, -u, -t, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_l(xb, -yb, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[1], v, u, t, 0.0, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_l(-xb, -yb, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[1], -v, -u, -t, 0.0, 0.0));
    }
}

/// Collect all CCCC paths.
fn cccc(x: f64, y: f64, phi: f64, paths: &mut Vec<ReedsSheppPath>) {
    // LpRupLumRm - type 2
    if let Some((t, u, v)) = lp_rup_lum_rm(x, y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[2], t, u, -u, v, 0.0));
    }
    if let Some((t, u, v)) = lp_rup_lum_rm(-x, y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[2], -t, -u, u, -v, 0.0));
    }
    if let Some((t, u, v)) = lp_rup_lum_rm(x, -y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[3], t, u, -u, v, 0.0));
    }
    if let Some((t, u, v)) = lp_rup_lum_rm(-x, -y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[3], -t, -u, u, -v, 0.0));
    }

    // LpRumLumRp - type 2
    if let Some((t, u, v)) = lp_rum_lum_rp(x, y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[2], t, u, u, v, 0.0));
    }
    if let Some((t, u, v)) = lp_rum_lum_rp(-x, y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[2], -t, -u, -u, -v, 0.0));
    }
    if let Some((t, u, v)) = lp_rum_lum_rp(x, -y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[3], t, u, u, v, 0.0));
    }
    if let Some((t, u, v)) = lp_rum_lum_rp(-x, -y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[3], -t, -u, -u, -v, 0.0));
    }
}

/// Collect all CCSC paths.
fn ccsc(x: f64, y: f64, phi: f64, paths: &mut Vec<ReedsSheppPath>) {
    let half_pi = 0.5 * PI;

    // LpRmSmLm - types 4,5
    if let Some((t, u, v)) = lp_rm_sm_lm(x, y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[4], t, -half_pi, u, v, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_lm(-x, y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[4], -t, half_pi, -u, -v, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_lm(x, -y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[5], t, -half_pi, u, v, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_lm(-x, -y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[5], -t, half_pi, -u, -v, 0.0));
    }

    // LpRmSmRm - types 8,9
    if let Some((t, u, v)) = lp_rm_sm_rm(x, y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[8], t, -half_pi, u, v, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_rm(-x, y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[8], -t, half_pi, -u, -v, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_rm(x, -y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[9], t, -half_pi, u, v, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_rm(-x, -y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[9], -t, half_pi, -u, -v, 0.0));
    }

    // backwards
    let xb = x * phi.cos() + y * phi.sin();
    let yb = x * phi.sin() - y * phi.cos();

    // LpRmSmLm backwards - types 6,7
    if let Some((t, u, v)) = lp_rm_sm_lm(xb, yb, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[6], v, u, -half_pi, t, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_lm(-xb, yb, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[6], -v, -u, half_pi, -t, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_lm(xb, -yb, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[7], v, u, -half_pi, t, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_lm(-xb, -yb, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[7], -v, -u, half_pi, -t, 0.0));
    }

    // LpRmSmRm backwards - types 10,11
    if let Some((t, u, v)) = lp_rm_sm_rm(xb, yb, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[10], v, u, -half_pi, t, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_rm(-xb, yb, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[10], -v, -u, half_pi, -t, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_rm(xb, -yb, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[11], v, u, -half_pi, t, 0.0));
    }
    if let Some((t, u, v)) = lp_rm_sm_rm(-xb, -yb, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[11], -v, -u, half_pi, -t, 0.0));
    }
}

/// Collect all CCSCC paths.
fn ccscc(x: f64, y: f64, phi: f64, paths: &mut Vec<ReedsSheppPath>) {
    let half_pi = 0.5 * PI;

    // types 16,17
    if let Some((t, u, v)) = lp_rm_s_lm_rp(x, y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[16], t, -half_pi, u, -half_pi, v));
    }
    if let Some((t, u, v)) = lp_rm_s_lm_rp(-x, y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[16], -t, half_pi, -u, half_pi, -v));
    }
    if let Some((t, u, v)) = lp_rm_s_lm_rp(x, -y, -phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[17], t, -half_pi, u, -half_pi, v));
    }
    if let Some((t, u, v)) = lp_rm_s_lm_rp(-x, -y, phi) {
        paths.push(ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[17], -t, half_pi, -u, half_pi, -v));
    }
}

// ---------------------------------------------------------------------------
// Top-level Reeds-Shepp path computation (normalised coordinates)
// ---------------------------------------------------------------------------

/// Find the shortest among a set of collected candidate paths.
fn find_shortest(
    x: f64,
    y: f64,
    phi: f64,
    best: &mut ReedsSheppPath,
    collector: fn(f64, f64, f64, &mut Vec<ReedsSheppPath>),
) {
    let mut candidates = Vec::new();
    collector(x, y, phi, &mut candidates);
    let mut l_min = best.length();
    for p in &candidates {
        let l = p.length();
        if l < l_min {
            *best = *p;
            l_min = l;
        }
    }
}

/// Compute the shortest Reeds-Shepp path in normalised coordinates.
pub fn reeds_shepp(x: f64, y: f64, phi: f64) -> ReedsSheppPath {
    let mut path = ReedsSheppPath::default();
    find_shortest(x, y, phi, &mut path, csc);
    find_shortest(x, y, phi, &mut path, ccc);
    find_shortest(x, y, phi, &mut path, cccc);
    find_shortest(x, y, phi, &mut path, ccsc);
    find_shortest(x, y, phi, &mut path, ccscc);
    path
}

/// Collect all feasible Reeds-Shepp paths in normalised coordinates.
pub fn reeds_shepp_all(x: f64, y: f64, phi: f64) -> Vec<ReedsSheppPath> {
    let mut paths = Vec::new();
    csc(x, y, phi, &mut paths);
    ccc(x, y, phi, &mut paths);
    cccc(x, y, phi, &mut paths);
    ccsc(x, y, phi, &mut paths);
    ccscc(x, y, phi, &mut paths);
    paths
}

// ---------------------------------------------------------------------------
// ReedsSheppStateSpace
// ---------------------------------------------------------------------------

/// Reeds-Shepp state space for shortest-path planning with reverse driving.
///
/// Reeds-Shepp paths extend Dubins paths by allowing the vehicle to drive
/// backwards.  Paths consist of up to five segments: arcs of maximum curvature
/// and straight lines, in both forward and reverse directions.
pub struct ReedsSheppStateSpace {
    /// Maximum curvature.
    kappa: f64,
    /// Inverse of maximum curvature (minimum turning radius).
    kappa_inv: f64,
    /// Discretization step size for path integration.
    discretization: f64,
}

impl ReedsSheppStateSpace {
    /// Create a new Reeds-Shepp state space.
    ///
    /// # Arguments
    /// * `kappa` – maximum curvature (must be positive)
    /// * `discretization` – step size for path discretization (default 0.1)
    ///
    /// # Panics
    /// Panics if `kappa` is not positive.
    pub fn new(kappa: f64, discretization: f64) -> Self {
        assert!(kappa > 0.0, "kappa must be positive");
        Self {
            kappa,
            kappa_inv: 1.0 / kappa,
            discretization,
        }
    }

    /// Compute the normalised shortest Reeds-Shepp path between two states.
    pub fn reeds_shepp(&self, state1: &State, state2: &State) -> ReedsSheppPath {
        let dx = state2.x - state1.x;
        let dy = state2.y - state1.y;
        let dth = state2.theta - state1.theta;
        let c = state1.theta.cos();
        let s = state1.theta.sin();
        let x = c * dx + s * dy;
        let y = -s * dx + c * dy;
        reeds_shepp(x * self.kappa, y * self.kappa, dth)
    }

    /// Compute all feasible Reeds-Shepp paths between two states.
    pub fn get_all_rs_paths(&self, state1: &State, state2: &State) -> Vec<ReedsSheppPath> {
        let dx = state2.x - state1.x;
        let dy = state2.y - state1.y;
        let dth = state2.theta - state1.theta;
        let c = state1.theta.cos();
        let s = state1.theta.sin();
        let x = c * dx + s * dy;
        let y = -s * dx + c * dy;
        reeds_shepp_all(x * self.kappa, y * self.kappa, dth)
    }

    /// Compute the Reeds-Shepp distance (arc length) between two states.
    pub fn get_distance(&self, state1: &State, state2: &State) -> f64 {
        self.kappa_inv * self.reeds_shepp(state1, state2).length()
    }

    /// Extract control segments from a normalised Reeds-Shepp path.
    pub fn extract_controls_from_path(&self, path: &ReedsSheppPath) -> Vec<Control> {
        let mut controls = Vec::new();
        for i in 0..5 {
            match path.segment_types[i] {
                ReedsSheppPathSegmentType::Nop => return controls,
                ReedsSheppPathSegmentType::Left => {
                    let control = Control {
                        delta_s: self.kappa_inv * path.lengths[i],
                        kappa: self.kappa,
                        sigma: 0.0,
                    };
                    if control.delta_s.abs() > RS_EPS {
                        controls.push(control);
                    }
                }
                ReedsSheppPathSegmentType::Right => {
                    let control = Control {
                        delta_s: self.kappa_inv * path.lengths[i],
                        kappa: -self.kappa,
                        sigma: 0.0,
                    };
                    if control.delta_s.abs() > RS_EPS {
                        controls.push(control);
                    }
                }
                ReedsSheppPathSegmentType::Straight => {
                    let control = Control {
                        delta_s: self.kappa_inv * path.lengths[i],
                        kappa: 0.0,
                        sigma: 0.0,
                    };
                    if control.delta_s.abs() > RS_EPS {
                        controls.push(control);
                    }
                }
            }
        }
        controls
    }
}

impl StateSpace for ReedsSheppStateSpace {
    fn discretization(&self) -> f64 {
        self.discretization
    }

    fn get_controls(&self, state1: &State, state2: &State) -> Vec<Control> {
        let rs_path = self.reeds_shepp(state1, state2);
        self.extract_controls_from_path(&rs_path)
    }

    fn get_all_controls(&self, state1: &State, state2: &State) -> Vec<Vec<Control>> {
        if (state1.x - state2.x).abs() < RS_EPS
            && (state1.y - state2.y).abs() < RS_EPS
            && (state1.theta - state2.theta).abs() < RS_EPS
        {
            return vec![vec![Control {
                delta_s: 0.0,
                kappa: self.kappa,
                sigma: 0.0,
            }]];
        }

        let all_rs_paths = self.get_all_rs_paths(state1, state2);
        let mut all_controls = Vec::new();
        for rs_path in &all_rs_paths {
            let control = self.extract_controls_from_path(rs_path);
            if !control.is_empty() {
                all_controls.push(control);
            }
        }
        all_controls
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utilities::HALF_PI;

    const TEST_EPS: f64 = 1e-4;

    // -----------------------------------------------------------------------
    // ReedsSheppPath tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_reeds_shepp_path_length() {
        let p = ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[0], 1.0, -2.0, 3.0, 0.0, 0.0);
        assert!((p.length() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_reeds_shepp_path_default_is_max() {
        let p = ReedsSheppPath::default();
        assert!(p.length() > 1e30);
    }

    #[test]
    fn test_reeds_shepp_path_get_length() {
        let p = ReedsSheppPath::new(REEDS_SHEPP_PATH_TYPE[0], 1.5, -2.5, 3.5, 0.0, 0.0);
        assert!((p.get_length(0) - 1.5).abs() < 1e-12);
        assert!((p.get_length(1) - (-2.5)).abs() < 1e-12);
        assert!((p.get_length(2) - 3.5).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // Shortest-path selection tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_reeds_shepp_same_point_same_heading() {
        let path = reeds_shepp(0.0, 0.0, 0.0);
        assert!(path.length() < TEST_EPS);
    }

    #[test]
    fn test_reeds_shepp_straight_ahead() {
        let path = reeds_shepp(5.0, 0.0, 0.0);
        assert!((path.length() - 5.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_reeds_shepp_straight_reverse() {
        let path = reeds_shepp(-5.0, 0.0, 0.0);
        assert!((path.length() - 5.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_reeds_shepp_shorter_than_dubins() {
        // Going backwards should be shorter than going around
        let path = reeds_shepp(-3.0, 0.0, 0.0);
        assert!(path.length() < 3.0 + TEST_EPS);
    }

    #[test]
    fn test_reeds_shepp_all_returns_paths() {
        let paths = reeds_shepp_all(3.0, 2.0, 1.0);
        assert!(!paths.is_empty());
        // Shortest among all should match the single-best result
        let best = reeds_shepp(3.0, 2.0, 1.0);
        let min_all = paths.iter().map(|p| p.length()).fold(f64::MAX, f64::min);
        assert!((min_all - best.length()).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // ReedsSheppStateSpace tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_state_space_distance_same_point() {
        let ss = ReedsSheppStateSpace::new(1.0, 0.1);
        let s = State::new(0.0, 0.0, 0.0, 0.0);
        let d = ss.get_distance(&s, &s);
        assert!(d < TEST_EPS);
    }

    #[test]
    fn test_state_space_straight_distance() {
        let ss = ReedsSheppStateSpace::new(1.0, 0.1);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(5.0, 0.0, 0.0, 0.0);
        let d = ss.get_distance(&s1, &s2);
        assert!((d - 5.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_reverse_distance() {
        let ss = ReedsSheppStateSpace::new(1.0, 0.1);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(-5.0, 0.0, 0.0, 0.0);
        let d = ss.get_distance(&s1, &s2);
        // Reeds-Shepp can reverse, so distance should be 5.0
        assert!((d - 5.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_get_controls_straight() {
        let ss = ReedsSheppStateSpace::new(1.0, 0.1);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(5.0, 0.0, 0.0, 0.0);
        let controls = ss.get_controls(&s1, &s2);
        assert!(!controls.is_empty());
        let total: f64 = controls.iter().map(|c| c.delta_s.abs()).sum();
        assert!((total - 5.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_get_path_reaches_goal() {
        let ss = ReedsSheppStateSpace::new(1.0, 0.1);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(5.0, 0.0, 0.0, 0.0);
        let path = ss.get_path(&s1, &s2);
        assert!(!path.is_empty());
        let last = path.last().unwrap();
        assert!((last.x - s2.x).abs() < TEST_EPS);
        assert!((last.y - s2.y).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_get_path_with_turn() {
        let ss = ReedsSheppStateSpace::new(1.0, 0.1);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(2.0, 2.0, HALF_PI, 0.0);
        let path = ss.get_path(&s1, &s2);
        assert!(!path.is_empty());
        let last = path.last().unwrap();
        assert!((last.x - s2.x).abs() < TEST_EPS);
        assert!((last.y - s2.y).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_get_path_reverse() {
        let ss = ReedsSheppStateSpace::new(1.0, 0.1);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(-3.0, 0.0, 0.0, 0.0);
        let path = ss.get_path(&s1, &s2);
        assert!(!path.is_empty());
        let last = path.last().unwrap();
        assert!((last.x - s2.x).abs() < TEST_EPS);
        assert!((last.y - s2.y).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_all_controls_same_state() {
        let ss = ReedsSheppStateSpace::new(1.0, 0.1);
        let s = State::new(1.0, 2.0, 0.5, 0.0);
        let all = ss.get_all_controls(&s, &s);
        assert_eq!(all.len(), 1);
        assert_eq!(all[0].len(), 1);
        assert!((all[0][0].delta_s).abs() < 1e-12);
    }

    #[test]
    fn test_state_space_all_controls_multiple_paths() {
        let ss = ReedsSheppStateSpace::new(1.0, 0.1);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(3.0, 2.0, 1.0, 0.0);
        let all = ss.get_all_controls(&s1, &s2);
        assert!(!all.is_empty());
    }

    #[test]
    fn test_state_space_high_curvature() {
        let ss = ReedsSheppStateSpace::new(2.0, 0.05);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(1.0, 1.0, PI, 0.0);
        let path = ss.get_path(&s1, &s2);
        assert!(!path.is_empty());
        let last = path.last().unwrap();
        assert!((last.x - s2.x).abs() < TEST_EPS);
        assert!((last.y - s2.y).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_get_all_paths() {
        let ss = ReedsSheppStateSpace::new(1.0, 0.1);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(3.0, 2.0, 1.0, 0.0);
        let all_paths = ss.get_all_paths(&s1, &s2);
        assert!(!all_paths.is_empty());
        for p in &all_paths {
            assert!(!p.is_empty());
        }
    }
}
