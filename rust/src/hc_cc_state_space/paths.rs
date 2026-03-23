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

use std::fmt;

use super::configuration::Configuration;
use super::hc_cc_circle::HcCcCircle;
use crate::state::{Control, State};
use crate::utilities::{get_epsilon, point_distance, sgn, twopify};

// ---------------------------------------------------------------------------
// Path type enumerations
// ---------------------------------------------------------------------------

/// CC-Dubins path types: E (Empty), S (Straight), T (Turn).
pub const NB_CC_DUBINS_PATHS: usize = 7;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CcDubinsPathType {
    E,
    S,
    T,
    TT,
    /// Dubins families
    TST,
    TTT,
    TTTT,
}

impl fmt::Display for CcDubinsPathType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::E => write!(f, "E"),
            Self::S => write!(f, "S"),
            Self::T => write!(f, "T"),
            Self::TT => write!(f, "TT"),
            Self::TST => write!(f, "TST"),
            Self::TTT => write!(f, "TTT"),
            Self::TTTT => write!(f, "TTTT"),
        }
    }
}

/// HC/CC-Reeds-Shepp path types: E (Empty), S (Straight), T (Turn), c (Cusp).
pub const NB_HC_CC_RS_PATHS: usize = 18;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HcCcRsPathType {
    E,
    S,
    T,
    TT,
    TcT,
    /// Reeds-Shepp families
    TcTcT,
    TcTT,
    TTcT,
    TST,
    TSTcT,
    TcTST,
    TcTSTcT,
    TTcTT,
    TcTTcT,
    TTT,
    TcST,
    TScT,
    TcScT,
}

impl fmt::Display for HcCcRsPathType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::E => write!(f, "E"),
            Self::S => write!(f, "S"),
            Self::T => write!(f, "T"),
            Self::TT => write!(f, "TT"),
            Self::TcT => write!(f, "TcT"),
            Self::TcTcT => write!(f, "TcTcT"),
            Self::TcTT => write!(f, "TcTT"),
            Self::TTcT => write!(f, "TTcT"),
            Self::TST => write!(f, "TST"),
            Self::TSTcT => write!(f, "TSTcT"),
            Self::TcTST => write!(f, "TcTST"),
            Self::TcTSTcT => write!(f, "TcTSTcT"),
            Self::TTcTT => write!(f, "TTcTT"),
            Self::TcTTcT => write!(f, "TcTTcT"),
            Self::TTT => write!(f, "TTT"),
            Self::TcST => write!(f, "TcST"),
            Self::TScT => write!(f, "TScT"),
            Self::TcScT => write!(f, "TcScT"),
        }
    }
}

// ---------------------------------------------------------------------------
// Path structs
// ---------------------------------------------------------------------------

/// A CC-Dubins path.
#[derive(Debug, Clone)]
pub struct CcDubinsPath {
    /// Start configuration.
    pub start: Configuration,
    /// End configuration.
    pub end: Configuration,
    /// Path type.
    pub path_type: CcDubinsPathType,
    /// Maximum curvature (unsigned).
    pub kappa: f64,
    /// Maximum sharpness (unsigned).
    pub sigma: f64,
    /// Intermediate configurations.
    pub qi1: Option<Configuration>,
    pub qi2: Option<Configuration>,
    pub qi3: Option<Configuration>,
    pub qi4: Option<Configuration>,
    /// Start, end, and intermediate circles.
    pub cstart: Option<HcCcCircle>,
    pub cend: Option<HcCcCircle>,
    pub ci1: Option<HcCcCircle>,
    pub ci2: Option<HcCcCircle>,
    /// Path length.
    pub length: f64,
}

impl CcDubinsPath {
    /// Create a new CC-Dubins path.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        start: Configuration,
        end: Configuration,
        path_type: CcDubinsPathType,
        kappa: f64,
        sigma: f64,
        qi1: Option<Configuration>,
        qi2: Option<Configuration>,
        qi3: Option<Configuration>,
        qi4: Option<Configuration>,
        cstart: Option<HcCcCircle>,
        cend: Option<HcCcCircle>,
        ci1: Option<HcCcCircle>,
        ci2: Option<HcCcCircle>,
        length: f64,
    ) -> Self {
        Self {
            start,
            end,
            path_type,
            kappa,
            sigma,
            qi1,
            qi2,
            qi3,
            qi4,
            cstart,
            cend,
            ci1,
            ci2,
            length,
        }
    }
}

impl fmt::Display for CcDubinsPath {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "CC_Dubins_Path: type {}, length {}, configurations {}",
            self.path_type, self.length, self.start
        )?;
        for qi in [&self.qi1, &self.qi2, &self.qi3, &self.qi4] {
            if let Some(q) = qi {
                write!(f, " -> {}", q)?;
            }
        }
        write!(f, " -> {}", self.end)
    }
}

/// An HC/CC-Reeds-Shepp path.
#[derive(Debug, Clone)]
pub struct HcCcRsPath {
    /// Start configuration.
    pub start: Configuration,
    /// End configuration.
    pub end: Configuration,
    /// Path type.
    pub path_type: HcCcRsPathType,
    /// Maximum curvature (unsigned).
    pub kappa: f64,
    /// Maximum sharpness (unsigned).
    pub sigma: f64,
    /// Intermediate configurations.
    pub qi1: Option<Configuration>,
    pub qi2: Option<Configuration>,
    pub qi3: Option<Configuration>,
    pub qi4: Option<Configuration>,
    /// Start, end, and intermediate circles.
    pub cstart: Option<HcCcCircle>,
    pub cend: Option<HcCcCircle>,
    pub ci1: Option<HcCcCircle>,
    pub ci2: Option<HcCcCircle>,
    /// Path length.
    pub length: f64,
}

impl HcCcRsPath {
    /// Create a new HC/CC-Reeds-Shepp path.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        start: Configuration,
        end: Configuration,
        path_type: HcCcRsPathType,
        kappa: f64,
        sigma: f64,
        qi1: Option<Configuration>,
        qi2: Option<Configuration>,
        qi3: Option<Configuration>,
        qi4: Option<Configuration>,
        cstart: Option<HcCcCircle>,
        cend: Option<HcCcCircle>,
        ci1: Option<HcCcCircle>,
        ci2: Option<HcCcCircle>,
        length: f64,
    ) -> Self {
        Self {
            start,
            end,
            path_type,
            kappa,
            sigma,
            qi1,
            qi2,
            qi3,
            qi4,
            cstart,
            cend,
            ci1,
            ci2,
            length,
        }
    }
}

impl fmt::Display for HcCcRsPath {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "HC_CC_RS_Path: type {}, length {}, configurations {}",
            self.path_type, self.length, self.start
        )?;
        for qi in [&self.qi1, &self.qi2, &self.qi3, &self.qi4] {
            if let Some(q) = qi {
                write!(f, " -> {}", q)?;
            }
        }
        write!(f, " -> {}", self.end)
    }
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Check whether two states are equal within tolerance.
pub fn state_equal(state1: &State, state2: &State) -> bool {
    if (state2.kappa - state1.kappa).abs() > get_epsilon() {
        return false;
    }
    if (twopify(state2.theta) - twopify(state1.theta)).abs() > get_epsilon() {
        return false;
    }
    if point_distance(state1.x, state1.y, state2.x, state2.y) > get_epsilon() {
        return false;
    }
    true
}

/// Reverse a control (negate delta_s, update kappa, negate sigma).
pub fn reverse_control(control: &mut Control) {
    control.delta_s = -control.delta_s;
    control.kappa += control.delta_s.abs() * control.sigma;
    control.sigma = -control.sigma;
}

/// Subtract `control2` from `control1` (same sigma direction required).
pub fn subtract_control(control1: &Control, control2: &Control) -> Control {
    debug_assert!(
        (sgn(control1.delta_s) * control1.sigma - sgn(control2.delta_s) * control2.sigma).abs()
            < get_epsilon()
    );
    Control {
        delta_s: control1.delta_s - control2.delta_s,
        kappa: control1.kappa,
        sigma: control1.sigma,
    }
}

/// Append a zero-length control segment.
pub fn empty_controls(controls: &mut Vec<Control>) {
    controls.push(Control {
        delta_s: 0.0,
        kappa: 0.0,
        sigma: 0.0,
    });
}

/// Append a straight-line control from `q1` to `q2`.
pub fn straight_controls(q1: &Configuration, q2: &Configuration, controls: &mut Vec<Control>) {
    let length = point_distance(q1.x, q1.y, q2.x, q2.y);
    let dot_product = q1.theta.cos() * (q2.x - q1.x) + q1.theta.sin() * (q2.y - q1.y);
    let d = sgn(dot_product);
    controls.push(Control {
        delta_s: d * length,
        kappa: 0.0,
        sigma: 0.0,
    });
}

/// Compute the driving direction sign from forward/order flags.
fn direction(forward: bool, order: bool) -> f64 {
    if (forward && order) || (!forward && !order) {
        1.0
    } else {
        -1.0
    }
}

/// Append controls for a Reeds-Shepp turn on circle `c` to reach `q`.
pub fn rs_turn_controls(
    c: &HcCcCircle,
    q: &Configuration,
    order: bool,
    controls: &mut Vec<Control>,
) {
    debug_assert!((c.param.kappa.abs() - q.kappa.abs()).abs() < get_epsilon());
    debug_assert!((c.param.sigma.abs() - f64::MAX).abs() < get_epsilon());
    let delta = c.deflection(q);
    let length_arc = c.param.kappa_inv.abs() * c.rs_circular_deflection(delta);
    let d = direction(c.forward, order);
    controls.push(Control {
        delta_s: d * length_arc,
        kappa: c.param.kappa,
        sigma: 0.0,
    });
}

/// Append controls for an HC turn on circle `c` to reach `q`.
pub fn hc_turn_controls(
    c: &HcCcCircle,
    q: &Configuration,
    order: bool,
    controls: &mut Vec<Control>,
) {
    debug_assert!((c.param.kappa.abs() - q.kappa.abs()).abs() < get_epsilon());
    let delta = c.deflection(q);
    let length_min = (c.param.kappa / c.param.sigma).abs();
    let length_arc = c.param.kappa_inv.abs() * c.hc_circular_deflection(delta);
    let d = direction(c.forward, order);

    if order {
        controls.push(Control {
            delta_s: d * length_min,
            kappa: 0.0,
            sigma: c.param.sigma,
        });
    }

    controls.push(Control {
        delta_s: d * length_arc,
        kappa: c.param.kappa,
        sigma: 0.0,
    });

    if !order {
        controls.push(Control {
            delta_s: d * length_min,
            kappa: c.param.kappa,
            sigma: -c.param.sigma,
        });
    }
}

/// Append controls for a CC elementary path (two clothoids, no arc).
///
/// Returns `true` if an elementary path exists, `false` otherwise.
pub fn cc_elementary_controls(
    c: &HcCcCircle,
    q: &Configuration,
    delta: f64,
    order: bool,
    controls: &mut Vec<Control>,
) -> bool {
    if let Some(sigma0) = c.cc_elementary_sharpness(q, delta) {
        let length = (delta / sigma0.abs()).sqrt();
        let d = direction(c.forward, order);

        controls.push(Control {
            delta_s: d * length,
            kappa: 0.0,
            sigma: sigma0,
        });
        controls.push(Control {
            delta_s: d * length,
            kappa: sigma0 * length,
            sigma: -sigma0,
        });
        true
    } else {
        false
    }
}

/// Append controls for a default CC turn (wind-up clothoid, arc, wind-down clothoid).
pub fn cc_default_controls(
    c: &HcCcCircle,
    q: &Configuration,
    delta: f64,
    order: bool,
    controls: &mut Vec<Control>,
) {
    let _ = q; // q is unused; delta already encodes the deflection
    let length_min = (c.param.kappa / c.param.sigma).abs();
    let length_arc = c.param.kappa_inv.abs() * c.cc_circular_deflection(delta);
    let d = direction(c.forward, order);

    controls.push(Control {
        delta_s: d * length_min,
        kappa: 0.0,
        sigma: c.param.sigma,
    });
    controls.push(Control {
        delta_s: d * length_arc,
        kappa: c.param.kappa,
        sigma: 0.0,
    });
    controls.push(Control {
        delta_s: d * length_min,
        kappa: c.param.kappa,
        sigma: -c.param.sigma,
    });
}

/// Append controls for a CC turn on circle `c` to reach `q`.
///
/// Selects between elementary, default, or straight controls depending on the
/// deflection angle.
pub fn cc_turn_controls(
    c: &HcCcCircle,
    q: &Configuration,
    order: bool,
    controls: &mut Vec<Control>,
) {
    debug_assert!(q.kappa.abs() < get_epsilon());
    let delta = c.deflection(q);

    // delta ≈ 0: straight line
    if delta < get_epsilon() {
        if order {
            straight_controls(&c.start, q, controls);
        } else {
            straight_controls(q, &c.start, controls);
        }
        return;
    }

    // 0 < delta < 2 * delta_min: try elementary, fall back to default
    if delta < 2.0 * c.param.delta_min {
        let mut controls_elementary = Vec::new();
        if cc_elementary_controls(c, q, delta, order, &mut controls_elementary) {
            let mut controls_default = Vec::new();
            cc_default_controls(c, q, delta, order, &mut controls_default);

            let length_elementary: f64 =
                controls_elementary.iter().map(|c| c.delta_s.abs()).sum();
            let length_default: f64 = controls_default.iter().map(|c| c.delta_s.abs()).sum();

            if length_elementary < length_default {
                controls.extend(controls_elementary);
            } else {
                controls.extend(controls_default);
            }
            return;
        }
        cc_default_controls(c, q, delta, order, controls);
        return;
    }

    // delta >= 2 * delta_min
    cc_default_controls(c, q, delta, order, controls);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hc_cc_state_space::hc_cc_circle::HcCcCircleParam;
    use crate::utilities::EPSILON;

    fn make_param() -> HcCcCircleParam {
        let kappa = 1.0;
        let sigma = 1.0;
        let mu = 0.5 * kappa * kappa / sigma;
        HcCcCircleParam {
            kappa,
            kappa_inv: 1.0 / kappa,
            sigma,
            radius: 1.2,
            mu,
            sin_mu: mu.sin(),
            cos_mu: mu.cos(),
            delta_min: mu,
        }
    }

    #[test]
    fn test_cc_dubins_path_type_display() {
        assert_eq!(format!("{}", CcDubinsPathType::E), "E");
        assert_eq!(format!("{}", CcDubinsPathType::TST), "TST");
        assert_eq!(format!("{}", CcDubinsPathType::TTTT), "TTTT");
    }

    #[test]
    fn test_hc_cc_rs_path_type_display() {
        assert_eq!(format!("{}", HcCcRsPathType::E), "E");
        assert_eq!(format!("{}", HcCcRsPathType::TcTSTcT), "TcTSTcT");
        assert_eq!(format!("{}", HcCcRsPathType::TcScT), "TcScT");
    }

    #[test]
    fn test_cc_dubins_path_display() {
        let start = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let end = Configuration::new(1.0, 0.0, 0.0, 0.0);
        let path = CcDubinsPath::new(
            start,
            end,
            CcDubinsPathType::S,
            1.0,
            1.0,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            1.0,
        );
        let s = format!("{}", path);
        assert!(s.contains("CC_Dubins_Path"));
        assert!(s.contains("S"));
    }

    #[test]
    fn test_hc_cc_rs_path_display() {
        let start = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let end = Configuration::new(1.0, 0.0, 0.0, 0.0);
        let path = HcCcRsPath::new(
            start,
            end,
            HcCcRsPathType::TST,
            1.0,
            1.0,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            2.0,
        );
        let s = format!("{}", path);
        assert!(s.contains("HC_CC_RS_Path"));
        assert!(s.contains("TST"));
    }

    #[test]
    fn test_state_equal_identical() {
        let s1 = State::new(1.0, 2.0, 0.5, 0.1);
        let s2 = State::new(1.0, 2.0, 0.5, 0.1);
        assert!(state_equal(&s1, &s2));
    }

    #[test]
    fn test_state_equal_different_kappa() {
        let s1 = State::new(1.0, 2.0, 0.5, 0.1);
        let s2 = State::new(1.0, 2.0, 0.5, 0.5);
        assert!(!state_equal(&s1, &s2));
    }

    #[test]
    fn test_state_equal_different_position() {
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(1.0, 0.0, 0.0, 0.0);
        assert!(!state_equal(&s1, &s2));
    }

    #[test]
    fn test_reverse_control() {
        let mut c = Control {
            delta_s: 1.0,
            kappa: 0.0,
            sigma: 2.0,
        };
        reverse_control(&mut c);
        assert!((c.delta_s - (-1.0)).abs() < 1e-12);
        assert!((c.kappa - 2.0).abs() < 1e-12);
        assert!((c.sigma - (-2.0)).abs() < 1e-12);
    }

    #[test]
    fn test_subtract_control() {
        let c1 = Control {
            delta_s: 3.0,
            kappa: 0.5,
            sigma: 1.0,
        };
        let c2 = Control {
            delta_s: 1.0,
            kappa: 0.5,
            sigma: 1.0,
        };
        let result = subtract_control(&c1, &c2);
        assert!((result.delta_s - 2.0).abs() < 1e-12);
        assert!((result.kappa - 0.5).abs() < 1e-12);
        assert!((result.sigma - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_empty_controls() {
        let mut controls = Vec::new();
        empty_controls(&mut controls);
        assert_eq!(controls.len(), 1);
        assert!((controls[0].delta_s).abs() < 1e-12);
    }

    #[test]
    fn test_straight_controls_forward() {
        let q1 = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let q2 = Configuration::new(2.0, 0.0, 0.0, 0.0);
        let mut controls = Vec::new();
        straight_controls(&q1, &q2, &mut controls);
        assert_eq!(controls.len(), 1);
        assert!((controls[0].delta_s - 2.0).abs() < EPSILON);
    }

    #[test]
    fn test_straight_controls_backward() {
        let q1 = Configuration::new(2.0, 0.0, 0.0, 0.0);
        let q2 = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let mut controls = Vec::new();
        straight_controls(&q1, &q2, &mut controls);
        assert_eq!(controls.len(), 1);
        assert!((controls[0].delta_s - (-2.0)).abs() < EPSILON);
    }

    #[test]
    fn test_direction() {
        assert!((direction(true, true) - 1.0).abs() < 1e-12);
        assert!((direction(true, false) - (-1.0)).abs() < 1e-12);
        assert!((direction(false, true) - (-1.0)).abs() < 1e-12);
        assert!((direction(false, false) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_cc_turn_controls_zero_deflection() {
        let param = make_param();
        let start = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let c = HcCcCircle::from_configuration(start, true, true, true, &param);
        // q has same theta → deflection ≈ 0
        let q = Configuration::new(1.0, 0.0, 0.0, 0.0);
        let mut controls = Vec::new();
        cc_turn_controls(&c, &q, true, &mut controls);
        assert!(!controls.is_empty());
    }

    #[test]
    fn test_nb_path_counts() {
        assert_eq!(NB_CC_DUBINS_PATHS, 7);
        assert_eq!(NB_HC_CC_RS_PATHS, 18);
    }
}
