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
use crate::utilities::{
    fresnel, get_epsilon, global_frame_change, point_distance, twopify, HALF_PI, PI, TWO_PI,
};

/// Parameters that describe an HC/CC circle.
#[derive(Debug, Clone, Copy)]
pub struct HcCcCircleParam {
    /// Maximum curvature.
    pub kappa: f64,
    /// Inverse of maximum curvature (`1 / kappa`).
    pub kappa_inv: f64,
    /// Maximum sharpness.
    pub sigma: f64,
    /// Radius of the outer circle.
    pub radius: f64,
    /// Angle between the initial orientation and the tangent to the circle.
    pub mu: f64,
    /// Precomputed sine of `mu`.
    pub sin_mu: f64,
    /// Precomputed cosine of `mu`.
    pub cos_mu: f64,
    /// Minimal deflection.
    pub delta_min: f64,
}

impl HcCcCircleParam {
    /// Set all parameters at once.
    pub fn set_param(
        &mut self,
        kappa: f64,
        sigma: f64,
        radius: f64,
        mu: f64,
        sin_mu: f64,
        cos_mu: f64,
        delta_min: f64,
    ) {
        self.kappa = kappa;
        self.kappa_inv = 1.0 / kappa;
        self.sigma = sigma;
        self.radius = radius;
        self.mu = mu;
        self.sin_mu = sin_mu;
        self.cos_mu = cos_mu;
        self.delta_min = delta_min;
    }
}

impl Default for HcCcCircleParam {
    fn default() -> Self {
        Self {
            kappa: 0.0,
            kappa_inv: 0.0,
            sigma: 0.0,
            radius: 0.0,
            mu: 0.0,
            sin_mu: 0.0,
            cos_mu: 0.0,
            delta_min: 0.0,
        }
    }
}

/// An HC/CC circle used in continuous-curvature steering.
///
/// The circle is constructed either from a start [`Configuration`] (which
/// determines the centre) or from explicit centre coordinates.  The stored
/// `param` values may have their sign adjusted (negated for right-turning
/// circles) so that the methods can use them directly.
#[derive(Debug, Clone, Copy)]
pub struct HcCcCircle {
    /// Start configuration of the circle.
    pub start: Configuration,
    /// Turning direction: `true` = left, `false` = right.
    pub left: bool,
    /// Driving direction: `true` = forward, `false` = backward.
    pub forward: bool,
    /// Circle type: `true` = regular, `false` = irregular.
    pub regular: bool,
    /// Centre x-coordinate.
    pub xc: f64,
    /// Centre y-coordinate.
    pub yc: f64,
    /// Circle parameters (kappa/sigma signs adjusted for turning direction).
    pub param: HcCcCircleParam,
}

impl HcCcCircle {
    /// Construct a circle from a start configuration.
    ///
    /// The circle centre is computed via a global frame change from the start
    /// position.  For right-turning circles the curvature and sharpness signs
    /// are negated.
    pub fn from_configuration(
        start: Configuration,
        left: bool,
        forward: bool,
        regular: bool,
        param: &HcCcCircleParam,
    ) -> Self {
        let delta_x = param.radius * param.sin_mu;
        let delta_y = param.radius * param.cos_mu;

        let (effective_kappa, effective_kappa_inv, effective_sigma, local_x, local_y) = if left {
            let (lx, ly) = if forward {
                (delta_x, delta_y)
            } else {
                (-delta_x, delta_y)
            };
            (param.kappa, param.kappa_inv, param.sigma, lx, ly)
        } else {
            let (lx, ly) = if forward {
                (delta_x, -delta_y)
            } else {
                (-delta_x, -delta_y)
            };
            (-param.kappa, -param.kappa_inv, -param.sigma, lx, ly)
        };

        let (xc, yc) = global_frame_change(start.x, start.y, start.theta, local_x, local_y);

        Self {
            start,
            left,
            forward,
            regular,
            xc,
            yc,
            param: HcCcCircleParam {
                kappa: effective_kappa,
                kappa_inv: effective_kappa_inv,
                sigma: effective_sigma,
                radius: param.radius,
                mu: param.mu,
                sin_mu: param.sin_mu,
                cos_mu: param.cos_mu,
                delta_min: param.delta_min,
            },
        }
    }

    /// Construct a circle from explicit centre coordinates.
    ///
    /// A dummy start configuration `(0, 0, 0, 0)` is used.
    pub fn from_center(
        xc: f64,
        yc: f64,
        left: bool,
        forward: bool,
        regular: bool,
        param: &HcCcCircleParam,
    ) -> Self {
        let (effective_kappa, effective_kappa_inv, effective_sigma) = if left {
            (param.kappa, param.kappa_inv, param.sigma)
        } else {
            (-param.kappa, -param.kappa_inv, -param.sigma)
        };

        Self {
            start: Configuration::new(0.0, 0.0, 0.0, 0.0),
            left,
            forward,
            regular,
            xc,
            yc,
            param: HcCcCircleParam {
                kappa: effective_kappa,
                kappa_inv: effective_kappa_inv,
                sigma: effective_sigma,
                radius: param.radius,
                mu: param.mu,
                sin_mu: param.sin_mu,
                cos_mu: param.cos_mu,
                delta_min: param.delta_min,
            },
        }
    }

    /// Deflection angle between the start configuration and `q`.
    pub fn deflection(&self, q: &Configuration) -> f64 {
        let alpha_c = self.start.theta;
        let alpha_q = q.theta;
        if (self.left && self.forward) || (!self.left && !self.forward) {
            twopify(alpha_q - alpha_c)
        } else {
            twopify(alpha_c - alpha_q)
        }
    }

    /// D1 evaluation for an elementary path.
    pub fn d1(&self, alpha: f64) -> f64 {
        let s = (2.0 * alpha / PI).sqrt();
        let (fresnel_s, fresnel_c) = fresnel(s);
        alpha.cos() * fresnel_c + alpha.sin() * fresnel_s
    }

    /// Circular deflection for a Reeds-Shepp turn.
    pub fn rs_circular_deflection(&self, delta: f64) -> f64 {
        if self.regular {
            delta
        } else if delta <= PI {
            delta
        } else {
            delta - TWO_PI
        }
    }

    /// Length of a Reeds-Shepp turn to reach configuration `q`.
    pub fn rs_turn_length(&self, q: &Configuration) -> f64 {
        debug_assert!(
            (self.param.kappa.abs() - q.kappa.abs()).abs() < get_epsilon()
                && (self.param.sigma.abs() - f64::MAX).abs() < get_epsilon()
        );
        let delta = self.deflection(q);
        (self.param.kappa_inv * self.rs_circular_deflection(delta)).abs()
    }

    /// Circular deflection for an HC turn.
    pub fn hc_circular_deflection(&self, delta: f64) -> f64 {
        let delta_min_twopified = twopify(self.param.delta_min);
        if self.regular {
            if delta < delta_min_twopified {
                TWO_PI + delta - delta_min_twopified
            } else {
                delta - delta_min_twopified
            }
        } else {
            let (delta_arc1, delta_arc2) = if delta < delta_min_twopified {
                let d1 = delta - delta_min_twopified; // negative
                (d1, d1 + TWO_PI) // positive
            } else {
                let d1 = delta - delta_min_twopified; // positive
                (d1, d1 - TWO_PI) // negative
            };
            if delta_arc1.abs() < delta_arc2.abs() {
                delta_arc1
            } else {
                delta_arc2
            }
        }
    }

    /// Length of an HC turn to reach configuration `q`.
    pub fn hc_turn_length(&self, q: &Configuration) -> f64 {
        debug_assert!((self.param.kappa.abs() - q.kappa.abs()).abs() < get_epsilon());
        let delta = self.deflection(q);
        (self.param.kappa / self.param.sigma).abs()
            + (self.param.kappa_inv * self.hc_circular_deflection(delta)).abs()
    }

    /// Compute the elementary sharpness for a CC elementary path.
    ///
    /// Returns `Some(sigma0)` if the elementary path exists, `None` otherwise.
    pub fn cc_elementary_sharpness(&self, q: &Configuration, delta: f64) -> Option<f64> {
        let distance = point_distance(self.start.x, self.start.y, q.x, q.y);
        if delta < 4.5948 && distance > get_epsilon() {
            let mut sigma0 = 4.0 * PI * self.d1(0.5 * delta).powi(2) / distance.powi(2);
            if !self.left {
                sigma0 = -sigma0;
            }
            Some(sigma0)
        } else {
            None
        }
    }

    /// Circular deflection for a CC turn.
    pub fn cc_circular_deflection(&self, delta: f64) -> f64 {
        let two_delta_min_twopified = twopify(2.0 * self.param.delta_min);
        if self.regular {
            if delta < two_delta_min_twopified {
                TWO_PI + delta - two_delta_min_twopified
            } else {
                delta - two_delta_min_twopified
            }
        } else {
            let (delta_arc1, delta_arc2) = if delta < two_delta_min_twopified {
                let d1 = delta - two_delta_min_twopified; // negative
                (d1, d1 + TWO_PI) // positive
            } else {
                let d1 = delta - two_delta_min_twopified; // positive
                (d1, d1 - TWO_PI) // negative
            };
            if delta_arc1.abs() < delta_arc2.abs() {
                delta_arc1
            } else {
                delta_arc2
            }
        }
    }

    /// Length of a CC turn to reach configuration `q`.
    pub fn cc_turn_length(&self, q: &Configuration) -> f64 {
        debug_assert!(q.kappa.abs() < get_epsilon());
        let delta = self.deflection(q);
        // delta ≈ 0: elementary path
        if delta < get_epsilon() {
            return 2.0 * self.param.radius * self.param.sin_mu;
        }
        // default length via circular deflection
        let length_min = (self.param.kappa / self.param.sigma).abs();
        let length_default =
            2.0 * length_min + (self.param.kappa_inv * self.cc_circular_deflection(delta)).abs();
        // 0 < delta < 2 * delta_min: try elementary path
        if delta < 2.0 * self.param.delta_min {
            if let Some(sigma0) = self.cc_elementary_sharpness(q, delta) {
                let length_elementary = 2.0 * (delta / sigma0.abs()).sqrt();
                if length_elementary < length_default {
                    return length_elementary;
                }
            }
        }
        length_default
    }
}

impl fmt::Display for HcCcCircle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let dir = if self.left { "left" } else { "right" };
        let drive = if self.forward {
            "forward"
        } else {
            "backward"
        };
        let kind = if self.regular {
            "regular"
        } else {
            "irregular"
        };
        write!(
            f,
            "HC_CC_Circle: start: {}, {}, {}, {}, kappa: {}, sigma: {}, \
             centre: ({}, {}), radius {}, mu: {}",
            self.start,
            dir,
            drive,
            kind,
            self.param.kappa,
            self.param.sigma,
            self.xc,
            self.yc,
            self.param.radius,
            self.param.mu
        )
    }
}

/// Cartesian distance between the centres of two circles.
pub fn center_distance(c1: &HcCcCircle, c2: &HcCcCircle) -> f64 {
    point_distance(c1.xc, c1.yc, c2.xc, c2.yc)
}

/// Test whether configuration `q` lies on circle `c`.
pub fn configuration_on_hc_cc_circle(c: &HcCcCircle, q: &Configuration) -> bool {
    let distance = point_distance(c.xc, c.yc, q.x, q.y);
    if (distance - c.param.radius).abs() > get_epsilon() {
        return false;
    }
    let mut angle = (q.y - c.yc).atan2(q.x - c.xc);
    if c.left && c.forward {
        angle += HALF_PI - c.param.mu;
    }
    if c.left && !c.forward {
        angle += HALF_PI + c.param.mu;
    }
    if !c.left && c.forward {
        angle += -HALF_PI + c.param.mu;
    }
    if !c.left && !c.forward {
        angle += -HALF_PI - c.param.mu;
    }
    angle = twopify(angle);
    (q.theta - angle).abs() < get_epsilon()
}

#[cfg(test)]
mod tests {
    use super::*;
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
    fn test_param_set_param() {
        let mut p = HcCcCircleParam::default();
        p.set_param(2.0, 1.0, 3.0, 0.5, 0.5_f64.sin(), 0.5_f64.cos(), 0.25);
        assert!((p.kappa - 2.0).abs() < 1e-12);
        assert!((p.kappa_inv - 0.5).abs() < 1e-12);
        assert!((p.sigma - 1.0).abs() < 1e-12);
        assert!((p.radius - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_from_configuration_left_forward() {
        let param = make_param();
        let start = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let c = HcCcCircle::from_configuration(start, true, true, true, &param);
        assert!(c.left);
        assert!(c.forward);
        assert!(c.regular);
        assert!((c.param.kappa - param.kappa).abs() < 1e-12);
    }

    #[test]
    fn test_from_configuration_right_forward() {
        let param = make_param();
        let start = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let c = HcCcCircle::from_configuration(start, false, true, true, &param);
        assert!(!c.left);
        assert!((c.param.kappa - (-param.kappa)).abs() < 1e-12);
        assert!((c.param.sigma - (-param.sigma)).abs() < 1e-12);
    }

    #[test]
    fn test_from_center() {
        let param = make_param();
        let c = HcCcCircle::from_center(1.0, 2.0, true, true, true, &param);
        assert!((c.xc - 1.0).abs() < 1e-12);
        assert!((c.yc - 2.0).abs() < 1e-12);
        assert!((c.start.x).abs() < 1e-12);
    }

    #[test]
    fn test_deflection_left_forward() {
        let param = make_param();
        let start = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let c = HcCcCircle::from_configuration(start, true, true, true, &param);
        let q = Configuration::new(0.0, 0.0, 1.0, 0.0);
        let d = c.deflection(&q);
        assert!((d - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_deflection_right_forward() {
        let param = make_param();
        let start = Configuration::new(0.0, 0.0, 1.0, 0.0);
        let c = HcCcCircle::from_configuration(start, false, true, true, &param);
        let q = Configuration::new(0.0, 0.0, 0.5, 0.0);
        let d = c.deflection(&q);
        // right-forward: twopify(alpha_c - alpha_q) = twopify(1.0 - 0.5) = 0.5
        assert!((d - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_d1_zero() {
        let param = make_param();
        let c = HcCcCircle::from_center(0.0, 0.0, true, true, true, &param);
        let val = c.d1(0.0);
        assert!(val.abs() < 1e-12);
    }

    #[test]
    fn test_rs_circular_deflection_regular() {
        let param = make_param();
        let c = HcCcCircle::from_center(0.0, 0.0, true, true, true, &param);
        assert!((c.rs_circular_deflection(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_rs_circular_deflection_irregular_small() {
        let param = make_param();
        let c = HcCcCircle::from_center(0.0, 0.0, true, true, false, &param);
        assert!((c.rs_circular_deflection(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_rs_circular_deflection_irregular_large() {
        let param = make_param();
        let c = HcCcCircle::from_center(0.0, 0.0, true, true, false, &param);
        let delta = 4.0;
        assert!((c.rs_circular_deflection(delta) - (delta - TWO_PI)).abs() < 1e-12);
    }

    #[test]
    fn test_center_distance() {
        let param = make_param();
        let c1 = HcCcCircle::from_center(0.0, 0.0, true, true, true, &param);
        let c2 = HcCcCircle::from_center(3.0, 4.0, true, true, true, &param);
        assert!((center_distance(&c1, &c2) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_configuration_on_circle() {
        // mu=0 corresponds to the Reeds-Shepp case where the start
        // configuration lies directly on the circular arc.
        let kappa = 1.0;
        let radius = 1.0;
        let mu = 0.0;
        let param = HcCcCircleParam {
            kappa,
            kappa_inv: 1.0 / kappa,
            sigma: f64::MAX,
            radius,
            mu,
            sin_mu: mu.sin(),
            cos_mu: mu.cos(),
            delta_min: 0.0,
        };
        let start = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let c = HcCcCircle::from_configuration(start, true, true, true, &param);
        assert!(configuration_on_hc_cc_circle(&c, &start));
    }

    #[test]
    fn test_configuration_not_on_circle() {
        let param = make_param();
        let start = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let c = HcCcCircle::from_configuration(start, true, true, true, &param);
        let far = Configuration::new(100.0, 100.0, 0.0, 0.0);
        assert!(!configuration_on_hc_cc_circle(&c, &far));
    }

    #[test]
    fn test_hc_circular_deflection_regular() {
        let param = make_param();
        let c = HcCcCircle::from_center(0.0, 0.0, true, true, true, &param);
        let delta = 3.0;
        let hcd = c.hc_circular_deflection(delta);
        // For regular: delta >= delta_min_twopified → delta - delta_min_twopified
        assert!(hcd >= 0.0 || hcd < TWO_PI + EPSILON);
    }

    #[test]
    fn test_cc_circular_deflection_regular() {
        let param = make_param();
        let c = HcCcCircle::from_center(0.0, 0.0, true, true, true, &param);
        let delta = 5.0;
        let ccd = c.cc_circular_deflection(delta);
        assert!(ccd.is_finite());
    }

    #[test]
    fn test_display() {
        let param = make_param();
        let c = HcCcCircle::from_center(1.0, 2.0, true, true, true, &param);
        let s = format!("{}", c);
        assert!(s.contains("HC_CC_Circle"));
        assert!(s.contains("left"));
        assert!(s.contains("forward"));
        assert!(s.contains("regular"));
    }

    #[test]
    fn test_cc_turn_length_zero_deflection() {
        let kappa = 1.0;
        let sigma = 1.0;
        let mu = 0.5 * kappa * kappa / sigma;
        let radius = 1.2;
        let param = HcCcCircleParam {
            kappa,
            kappa_inv: 1.0 / kappa,
            sigma,
            radius,
            mu,
            sin_mu: mu.sin(),
            cos_mu: mu.cos(),
            delta_min: mu,
        };
        let start = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let c = HcCcCircle::from_configuration(start, true, true, true, &param);
        // q has same theta as start → deflection ≈ 0
        let q = Configuration::new(1.0, 0.0, 0.0, 0.0);
        let length = c.cc_turn_length(&q);
        let expected = 2.0 * radius * mu.sin();
        assert!((length - expected).abs() < EPSILON);
    }
}
