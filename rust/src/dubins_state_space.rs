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

use crate::base_state_space::{integrate_ode, StateSpace};
use crate::state::{Control, State};
use crate::utilities::{twopify, TWO_PI};

const DUBINS_EPS: f64 = 1e-6;
const DUBINS_ZERO: f64 = -1e-9;

// ---------------------------------------------------------------------------
// Dubins path types and structures
// ---------------------------------------------------------------------------

/// Type of a single Dubins path segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DubinsPathSegmentType {
    /// Turn left (positive curvature).
    Left = 0,
    /// Drive straight (zero curvature).
    Straight = 1,
    /// Turn right (negative curvature).
    Right = 2,
}

/// The six canonical Dubins path families: LSL, RSR, RSL, LSR, RLR, LRL.
const DUBINS_PATH_TYPE: [[DubinsPathSegmentType; 3]; 6] = {
    use DubinsPathSegmentType::*;
    [
        [Left, Straight, Left],   // LSL
        [Right, Straight, Right], // RSR
        [Right, Straight, Left],  // RSL
        [Left, Straight, Right],  // LSR
        [Right, Left, Right],     // RLR
        [Left, Right, Left],      // LRL
    ]
};

/// A Dubins path consisting of three segments.
///
/// Each segment is either a left turn, a straight line, or a right turn.
/// Lengths are in *normalised* units (multiply by `1/kappa` to get actual
/// arc length).
#[derive(Debug, Clone, Copy)]
pub struct DubinsPath {
    /// Segment types (e.g. Left–Straight–Left).
    pub segment_types: [DubinsPathSegmentType; 3],
    /// Length of each segment in normalised units.
    pub lengths: [f64; 3],
}

impl DubinsPath {
    /// Create a new Dubins path with the given segment types and lengths.
    pub fn new(segment_types: [DubinsPathSegmentType; 3], t: f64, p: f64, q: f64) -> Self {
        Self {
            segment_types,
            lengths: [t, p, q],
        }
    }

    /// Total normalised length of the path.
    pub fn length(&self) -> f64 {
        self.lengths[0] + self.lengths[1] + self.lengths[2]
    }
}

impl Default for DubinsPath {
    /// Sentinel path with maximum length (represents "no valid path found").
    fn default() -> Self {
        Self {
            segment_types: DUBINS_PATH_TYPE[0],
            lengths: [f64::MAX, 0.0, 0.0],
        }
    }
}

// ---------------------------------------------------------------------------
// The six Dubins path family computations
// ---------------------------------------------------------------------------

fn dubins_lsl(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sa - sb));
    if tmp >= DUBINS_ZERO {
        let theta = (cb - ca).atan2(d + sa - sb);
        let t = twopify(-alpha + theta);
        let p = tmp.max(0.0).sqrt();
        let q = twopify(beta - theta);
        DubinsPath::new(DUBINS_PATH_TYPE[0], t, p, q)
    } else {
        DubinsPath::default()
    }
}

fn dubins_rsr(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sb - sa));
    if tmp >= DUBINS_ZERO {
        let theta = (ca - cb).atan2(d - sa + sb);
        let t = twopify(alpha - theta);
        let p = tmp.max(0.0).sqrt();
        let q = twopify(-beta + theta);
        DubinsPath::new(DUBINS_PATH_TYPE[1], t, p, q)
    } else {
        DubinsPath::default()
    }
}

fn dubins_rsl(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let tmp = d * d - 2.0 + 2.0 * (ca * cb + sa * sb - d * (sa + sb));
    if tmp >= DUBINS_ZERO {
        let p = tmp.max(0.0).sqrt();
        let theta = (-ca - cb).atan2(d - sa - sb) - p.atan2(-2.0);
        let t = twopify(alpha - theta);
        let q = twopify(beta - theta);
        DubinsPath::new(DUBINS_PATH_TYPE[2], t, p, q)
    } else {
        DubinsPath::default()
    }
}

fn dubins_lsr(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let tmp = -2.0 + d * d + 2.0 * (ca * cb + sa * sb + d * (sa + sb));
    if tmp >= DUBINS_ZERO {
        let p = tmp.max(0.0).sqrt();
        let theta = (ca + cb).atan2(d + sa + sb) - (-2.0f64).atan2(p);
        let t = twopify(-alpha + theta);
        let q = twopify(-beta + theta);
        DubinsPath::new(DUBINS_PATH_TYPE[3], t, p, q)
    } else {
        DubinsPath::default()
    }
}

fn dubins_rlr(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let tmp = (6.0 - d * d + 2.0 * (ca * cb + sa * sb + d * (sa - sb))) / 8.0;
    if tmp.abs() <= 1.0 {
        let p = TWO_PI - tmp.acos();
        let theta = (ca - cb).atan2(d - sa + sb);
        let t = twopify(alpha - theta + p / 2.0);
        let q = twopify(alpha - beta - t + p);
        DubinsPath::new(DUBINS_PATH_TYPE[4], t, p, q)
    } else {
        DubinsPath::default()
    }
}

fn dubins_lrl(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let tmp = (6.0 - d * d + 2.0 * (ca * cb + sa * sb - d * (sa - sb))) / 8.0;
    if tmp.abs() <= 1.0 {
        let p = TWO_PI - tmp.acos();
        let theta = (-ca + cb).atan2(d + sa - sb);
        let t = twopify(-alpha + theta + p / 2.0);
        let q = twopify(beta - alpha - t + p);
        DubinsPath::new(DUBINS_PATH_TYPE[5], t, p, q)
    } else {
        DubinsPath::default()
    }
}

/// Find the shortest Dubins path in normalised coordinates.
///
/// # Arguments
/// * `d` – normalised distance between start and goal
/// * `alpha` – normalised start heading
/// * `beta` – normalised goal heading
pub fn dubins_shortest_path(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    if d < DUBINS_EPS && (alpha - beta).abs() < DUBINS_EPS {
        return DubinsPath::new(DUBINS_PATH_TYPE[0], 0.0, d, 0.0);
    }

    let mut best = dubins_lsl(d, alpha, beta);
    let mut min_length = best.length();

    let candidates = [
        dubins_rsr(d, alpha, beta),
        dubins_rsl(d, alpha, beta),
        dubins_lsr(d, alpha, beta),
        dubins_rlr(d, alpha, beta),
        dubins_lrl(d, alpha, beta),
    ];

    for candidate in &candidates {
        let len = candidate.length();
        if len < min_length {
            min_length = len;
            best = *candidate;
        }
    }

    best
}

// ---------------------------------------------------------------------------
// DubinsStateSpace
// ---------------------------------------------------------------------------

/// Dubins state space for curvature-bounded shortest-path planning.
///
/// Dubins paths consist of at most three segments: arcs of maximum curvature
/// and straight lines.  This state space computes the shortest such path
/// between two configurations (x, y, θ).
pub struct DubinsStateSpace {
    /// Maximum curvature.
    kappa: f64,
    /// Inverse of maximum curvature (minimum turning radius).
    kappa_inv: f64,
    /// Discretization step size for path integration.
    discretization: f64,
    /// Whether the vehicle drives forwards (`true`) or backwards (`false`).
    forwards: bool,
}

impl DubinsStateSpace {
    /// Create a new Dubins state space.
    ///
    /// # Arguments
    /// * `kappa` – maximum curvature (must be positive)
    /// * `discretization` – step size for path discretization (default 0.1)
    /// * `forwards` – driving direction
    ///
    /// # Panics
    /// Panics if `kappa` is not positive.
    pub fn new(kappa: f64, discretization: f64, forwards: bool) -> Self {
        assert!(kappa > 0.0, "kappa must be positive");
        Self {
            kappa,
            kappa_inv: 1.0 / kappa,
            discretization,
            forwards,
        }
    }

    /// Compute the normalised Dubins path between two states.
    pub fn dubins(&self, state1: &State, state2: &State) -> DubinsPath {
        let dx = state2.x - state1.x;
        let dy = state2.y - state1.y;
        let th = dy.atan2(dx);
        let d = (dx * dx + dy * dy).sqrt() * self.kappa;
        let alpha = twopify(state1.theta - th);
        let beta = twopify(state2.theta - th);
        dubins_shortest_path(d, alpha, beta)
    }

    /// Compute the Dubins distance (arc length) between two states.
    pub fn get_distance(&self, state1: &State, state2: &State) -> f64 {
        if self.forwards {
            self.kappa_inv * self.dubins(state1, state2).length()
        } else {
            self.kappa_inv * self.dubins(state2, state1).length()
        }
    }

    /// Core control computation shared between forwards and reverse.
    fn get_controls_core(
        &self,
        state1: &State,
        state2: &State,
        reverse_direction: bool,
    ) -> Vec<Control> {
        let path = if !reverse_direction {
            self.dubins(state1, state2)
        } else {
            self.dubins(state2, state1)
        };

        let mut controls = Vec::new();
        for i in 0..3 {
            let control = match path.segment_types[i] {
                DubinsPathSegmentType::Left => Control {
                    delta_s: self.kappa_inv * path.lengths[i],
                    kappa: self.kappa,
                    sigma: 0.0,
                },
                DubinsPathSegmentType::Straight => Control {
                    delta_s: self.kappa_inv * path.lengths[i],
                    kappa: 0.0,
                    sigma: 0.0,
                },
                DubinsPathSegmentType::Right => Control {
                    delta_s: self.kappa_inv * path.lengths[i],
                    kappa: -self.kappa,
                    sigma: 0.0,
                },
            };
            if control.delta_s.abs() > DUBINS_EPS {
                controls.push(control);
            }
        }

        if reverse_direction {
            controls.reverse();
            for control in &mut controls {
                control.delta_s = -control.delta_s;
            }
        }

        controls
    }

    /// Compute controls for reverse driving.
    pub fn get_controls_reverse(&self, state1: &State, state2: &State) -> Vec<Control> {
        self.get_controls_core(state1, state2, self.forwards)
    }

    /// Interpolate along a path at multiple sorted parameter values.
    ///
    /// This is more efficient than calling [`StateSpace::interpolate`] repeatedly
    /// because it traverses the control sequence only once.
    ///
    /// **The values in `ts` must be sorted in non-decreasing order.**
    pub fn interpolate_multi(
        &self,
        state: &State,
        controls: &[Control],
        ts: &[f64],
    ) -> Vec<State> {
        if controls.is_empty() || ts.is_empty() {
            return Vec::new();
        }

        let s_path: f64 = controls.iter().map(|c| c.delta_s.abs()).sum();
        let s_inters: Vec<f64> = ts.iter().map(|&t| t.clamp(0.0, 1.0) * s_path).collect();

        let mut state_inters = Vec::with_capacity(ts.len());
        let mut state_curr = *state;
        let mut s = 0.0;
        let mut idx = 0;

        for control in controls {
            if idx >= s_inters.len() {
                break;
            }
            let abs_delta_s = control.delta_s.abs();
            while idx < s_inters.len() {
                let s_inter = s_inters[idx];
                if s_inter - s > abs_delta_s {
                    break;
                }
                let state_inter = integrate_ode(&state_curr, control, s_inter - s);
                state_inters.push(state_inter);
                idx += 1;
            }
            state_curr = integrate_ode(&state_curr, control, abs_delta_s);
            s += abs_delta_s;
        }

        state_inters
    }
}

impl StateSpace for DubinsStateSpace {
    fn discretization(&self) -> f64 {
        self.discretization
    }

    fn get_controls(&self, state1: &State, state2: &State) -> Vec<Control> {
        self.get_controls_core(state1, state2, !self.forwards)
    }

    fn get_all_controls(&self, state1: &State, state2: &State) -> Vec<Vec<Control>> {
        vec![self.get_controls(state1, state2)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utilities::{HALF_PI, PI};

    const TEST_EPS: f64 = 1e-4;

    // -----------------------------------------------------------------------
    // DubinsPath tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_dubins_path_length() {
        let p = DubinsPath::new(DUBINS_PATH_TYPE[0], 1.0, 2.0, 3.0);
        assert!((p.length() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_dubins_path_default_is_max() {
        let p = DubinsPath::default();
        assert!(p.length() > 1e30);
    }

    // -----------------------------------------------------------------------
    // Shortest-path selection tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_dubins_same_point_same_heading() {
        let path = dubins_shortest_path(0.0, 0.0, 0.0);
        assert!(path.length() < DUBINS_EPS);
    }

    #[test]
    fn test_dubins_straight_ahead() {
        // Going straight: d > 0, alpha = 0, beta = 0 → LSL with t=0, p=d, q=0
        let d = 5.0;
        let path = dubins_shortest_path(d, 0.0, 0.0);
        assert!((path.length() - d).abs() < TEST_EPS);
    }

    #[test]
    fn test_dubins_lsl_is_found() {
        // A case that should produce an LSL path
        let path = dubins_shortest_path(4.0, 0.1, 0.1);
        assert!(path.length() < f64::MAX);
    }

    // -----------------------------------------------------------------------
    // DubinsStateSpace tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_state_space_distance_same_point() {
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s = State::new(0.0, 0.0, 0.0, 0.0);
        let d = ss.get_distance(&s, &s);
        assert!(d < TEST_EPS);
    }

    #[test]
    fn test_state_space_straight_distance() {
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(5.0, 0.0, 0.0, 0.0);
        let d = ss.get_distance(&s1, &s2);
        assert!((d - 5.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_distance_symmetry() {
        // For forwards, dubins(s1, s2) is used; distance may differ from dubins(s2, s1)
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(3.0, 2.0, HALF_PI, 0.0);
        let d_fwd = ss.get_distance(&s1, &s2);
        assert!(d_fwd > 0.0);
    }

    #[test]
    fn test_state_space_get_controls_straight() {
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(5.0, 0.0, 0.0, 0.0);
        let controls = ss.get_controls(&s1, &s2);
        assert!(!controls.is_empty());
        let total: f64 = controls.iter().map(|c| c.delta_s.abs()).sum();
        assert!((total - 5.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_get_path_reaches_goal() {
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(5.0, 0.0, 0.0, 0.0);
        let path = ss.get_path(&s1, &s2);
        assert!(!path.is_empty());
        let last = path.last().unwrap();
        assert!((last.x - s2.x).abs() < TEST_EPS);
        assert!((last.y - s2.y).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_get_path_turn() {
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(2.0, 2.0, HALF_PI, 0.0);
        let path = ss.get_path(&s1, &s2);
        assert!(!path.is_empty());
        let last = path.last().unwrap();
        assert!((last.x - s2.x).abs() < TEST_EPS);
        assert!((last.y - s2.y).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_interpolate_boundaries() {
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(5.0, 0.0, 0.0, 0.0);
        let controls = ss.get_controls(&s1, &s2);

        let at_start = ss.interpolate(&s1, &controls, 0.0);
        assert!(at_start.x.abs() < TEST_EPS);

        let at_end = ss.interpolate(&s1, &controls, 1.0);
        assert!((at_end.x - 5.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_interpolate_multi() {
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(4.0, 0.0, 0.0, 0.0);
        let controls = ss.get_controls(&s1, &s2);

        let ts = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let states = ss.interpolate_multi(&s1, &controls, &ts);
        assert_eq!(states.len(), 5);
        assert!(states[0].x.abs() < TEST_EPS);
        assert!((states[2].x - 2.0).abs() < TEST_EPS);
        assert!((states[4].x - 4.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_state_space_reverse() {
        let ss = DubinsStateSpace::new(1.0, 0.1, false);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(5.0, 0.0, 0.0, 0.0);
        let d = ss.get_distance(&s1, &s2);
        assert!(d > 0.0);
        let controls = ss.get_controls(&s1, &s2);
        assert!(!controls.is_empty());
        // Reverse controls should have negative delta_s
        for c in &controls {
            assert!(c.delta_s < 0.0 || c.delta_s.abs() < DUBINS_EPS);
        }
    }

    #[test]
    fn test_state_space_high_curvature() {
        let ss = DubinsStateSpace::new(2.0, 0.05, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(1.0, 1.0, PI, 0.0);
        let path = ss.get_path(&s1, &s2);
        assert!(!path.is_empty());
        let last = path.last().unwrap();
        assert!((last.x - s2.x).abs() < TEST_EPS);
        assert!((last.y - s2.y).abs() < TEST_EPS);
    }

    #[test]
    fn test_dubins_all_six_families() {
        // Exercise all six path families by choosing configurations that
        // favour each one as shortest.
        let configs: [(f64, f64, f64); 6] = [
            (4.0, 0.3, 0.3),   // LSL
            (4.0, -0.3, -0.3), // RSR
            (4.0, -1.0, 1.0),  // RSL
            (4.0, 1.0, -1.0),  // LSR
            (0.5, 2.5, -2.5),  // RLR
            (0.5, -2.5, 2.5),  // LRL
        ];
        for (d, alpha, beta) in &configs {
            let alpha = twopify(*alpha);
            let beta = twopify(*beta);
            let path = dubins_shortest_path(*d, alpha, beta);
            assert!(
                path.length() < f64::MAX,
                "No valid path for d={}, alpha={}, beta={}",
                d,
                alpha,
                beta
            );
        }
    }

    #[test]
    fn test_get_controls_reverse_vs_forward() {
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(3.0, 2.0, HALF_PI, 0.0);

        let fwd = ss.get_controls(&s1, &s2);
        let rev = ss.get_controls_reverse(&s1, &s2);

        // Both should be non-empty but produce different sign on delta_s
        assert!(!fwd.is_empty());
        assert!(!rev.is_empty());

        // Forward controls have positive delta_s
        for c in &fwd {
            assert!(c.delta_s > 0.0);
        }
        // Reverse controls have negative delta_s
        for c in &rev {
            assert!(c.delta_s < 0.0);
        }
    }

    #[test]
    fn test_semicircle_path() {
        // Semicircle: (0,0,0) → (0,2,π) with kappa=1 is a left turn of π
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(0.0, 2.0, PI, 0.0);
        let d = ss.get_distance(&s1, &s2);
        assert!((d - PI).abs() < TEST_EPS);
    }

    #[test]
    fn test_interpolate_multi_empty() {
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s = State::new(0.0, 0.0, 0.0, 0.0);
        assert!(ss.interpolate_multi(&s, &[], &[0.5]).is_empty());
        assert!(ss
            .interpolate_multi(
                &s,
                &[Control {
                    delta_s: 1.0,
                    kappa: 0.0,
                    sigma: 0.0
                }],
                &[]
            )
            .is_empty());
    }

    #[test]
    fn test_get_all_controls() {
        let ss = DubinsStateSpace::new(1.0, 0.1, true);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(3.0, 0.0, 0.0, 0.0);
        let all = ss.get_all_controls(&s1, &s2);
        assert_eq!(all.len(), 1);
        assert!(!all[0].is_empty());
    }
}
