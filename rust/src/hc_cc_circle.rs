use std::f64::consts::PI;
use crate::utilities::{TWO_PI, HALF_PI, get_epsilon, point_distance, twopify, fresnel, global_frame_change};
use crate::configuration::Configuration;

/// Parameters shared by all HC/CC circles.
#[derive(Clone, Copy, Debug, Default)]
pub struct HcCcCircleParam {
    pub kappa: f64,
    pub kappa_inv: f64,
    pub sigma: f64,
    pub radius: f64,
    pub mu: f64,
    pub sin_mu: f64,
    pub cos_mu: f64,
    pub delta_min: f64,
}

impl HcCcCircleParam {
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

/// HC/CC circle defined by a start configuration and turning parameters.
#[derive(Clone, Debug)]
pub struct HcCcCircle {
    // inherited from HcCcCircleParam
    pub kappa: f64,
    pub kappa_inv: f64,
    pub sigma: f64,
    pub radius: f64,
    pub mu: f64,
    pub sin_mu: f64,
    pub cos_mu: f64,
    pub delta_min: f64,
    // own fields
    pub start: Configuration,
    pub left: bool,
    pub forward: bool,
    pub regular: bool,
    pub xc: f64,
    pub yc: f64,
}

impl HcCcCircle {
    /// Construct from a start Configuration.
    pub fn from_configuration(
        start: Configuration,
        left: bool,
        forward: bool,
        regular: bool,
        param: &HcCcCircleParam,
    ) -> Self {
        let delta_x = param.radius * param.sin_mu;
        let delta_y = param.radius * param.cos_mu;

        let (kappa, kappa_inv, sigma, xc, yc) = if left {
            let (xc, yc) = if forward {
                global_frame_change(start.x, start.y, start.theta, delta_x, delta_y)
            } else {
                global_frame_change(start.x, start.y, start.theta, -delta_x, delta_y)
            };
            (param.kappa, param.kappa_inv, param.sigma, xc, yc)
        } else {
            let (xc, yc) = if forward {
                global_frame_change(start.x, start.y, start.theta, delta_x, -delta_y)
            } else {
                global_frame_change(start.x, start.y, start.theta, -delta_x, -delta_y)
            };
            (-param.kappa, -param.kappa_inv, -param.sigma, xc, yc)
        };

        Self {
            kappa,
            kappa_inv,
            sigma,
            radius: param.radius,
            mu: param.mu,
            sin_mu: param.sin_mu,
            cos_mu: param.cos_mu,
            delta_min: param.delta_min,
            start,
            left,
            forward,
            regular,
            xc,
            yc,
        }
    }

    /// Construct from circle center coordinates.
    pub fn from_center(
        xc: f64,
        yc: f64,
        left: bool,
        forward: bool,
        regular: bool,
        param: &HcCcCircleParam,
    ) -> Self {
        let (kappa, kappa_inv, sigma) = if left {
            (param.kappa, param.kappa_inv, param.sigma)
        } else {
            (-param.kappa, -param.kappa_inv, -param.sigma)
        };
        Self {
            kappa,
            kappa_inv,
            sigma,
            radius: param.radius,
            mu: param.mu,
            sin_mu: param.sin_mu,
            cos_mu: param.cos_mu,
            delta_min: param.delta_min,
            start: Configuration::new(0.0, 0.0, 0.0, 0.0),
            left,
            forward,
            regular,
            xc,
            yc,
        }
    }

    /// Angle between start configuration of circle and configuration q.
    pub fn deflection(&self, q: &Configuration) -> f64 {
        let alpha_c = self.start.theta;
        let alpha_q = q.theta;
        if (self.left && self.forward) || (!self.left && !self.forward) {
            twopify(alpha_q - alpha_c)
        } else {
            twopify(alpha_c - alpha_q)
        }
    }

    fn d1(&self, alpha: f64) -> f64 {
        let s = (2.0 * alpha / PI).sqrt();
        let (fresnel_s, fresnel_c) = fresnel(s);
        alpha.cos() * fresnel_c + alpha.sin() * fresnel_s
    }

    pub fn rs_circular_deflection(&self, delta: f64) -> f64 {
        if self.regular {
            delta
        } else if delta <= PI {
            delta
        } else {
            delta - TWO_PI
        }
    }

    pub fn rs_turn_length(&self, q: &Configuration) -> f64 {
        debug_assert!((self.kappa.abs() - q.kappa.abs()).abs() < get_epsilon());
        let delta = self.deflection(q);
        (self.kappa_inv * self.rs_circular_deflection(delta)).abs()
    }

    pub fn hc_circular_deflection(&self, delta: f64) -> f64 {
        let delta_min_twopified = twopify(self.delta_min);
        if self.regular {
            if delta < delta_min_twopified {
                TWO_PI + delta - delta_min_twopified
            } else {
                delta - delta_min_twopified
            }
        } else {
            let delta_arc1 = delta - delta_min_twopified;
            let delta_arc2 = if delta < delta_min_twopified {
                delta_arc1 + TWO_PI
            } else {
                delta_arc1 - TWO_PI
            };
            if delta_arc1.abs() < delta_arc2.abs() { delta_arc1 } else { delta_arc2 }
        }
    }

    pub fn hc_turn_length(&self, q: &Configuration) -> f64 {
        debug_assert!((self.kappa.abs() - q.kappa.abs()).abs() < get_epsilon());
        let delta = self.deflection(q);
        (self.kappa / self.sigma).abs() + (self.kappa_inv * self.hc_circular_deflection(delta)).abs()
    }

    /// Returns (success, sigma0).
    pub fn cc_elementary_sharpness(&self, q: &Configuration, delta: f64) -> (bool, f64) {
        let distance = point_distance(self.start.x, self.start.y, q.x, q.y);
        if delta < 4.5948 && distance > get_epsilon() {
            let mut sigma0 = 4.0 * PI * self.d1(0.5 * delta).powi(2) / distance.powi(2);
            if !self.left {
                sigma0 = -sigma0;
            }
            (true, sigma0)
        } else {
            (false, 0.0)
        }
    }

    pub fn cc_circular_deflection(&self, delta: f64) -> f64 {
        let two_delta_min_twopified = twopify(2.0 * self.delta_min);
        if self.regular {
            if delta < two_delta_min_twopified {
                TWO_PI + delta - two_delta_min_twopified
            } else {
                delta - two_delta_min_twopified
            }
        } else {
            let delta_arc1 = delta - two_delta_min_twopified;
            let delta_arc2 = if delta < two_delta_min_twopified {
                delta_arc1 + TWO_PI
            } else {
                delta_arc1 - TWO_PI
            };
            if delta_arc1.abs() < delta_arc2.abs() { delta_arc1 } else { delta_arc2 }
        }
    }

    pub fn cc_turn_length(&self, q: &Configuration) -> f64 {
        debug_assert!(q.kappa.abs() < get_epsilon());
        let delta = self.deflection(q);
        if delta < get_epsilon() {
            return 2.0 * self.radius * self.sin_mu;
        }
        let length_min = (self.kappa / self.sigma).abs();
        let length_default =
            2.0 * length_min + (self.kappa_inv * self.cc_circular_deflection(delta)).abs();
        if delta < 2.0 * self.delta_min {
            let (success, sigma0) = self.cc_elementary_sharpness(q, delta);
            if success {
                let length_elementary = 2.0 * (delta / sigma0.abs()).sqrt();
                return length_elementary.min(length_default);
            }
        }
        length_default
    }
}

/// Cartesian distance between the centres of two circles.
pub fn center_distance(c1: &HcCcCircle, c2: &HcCcCircle) -> f64 {
    ((c2.xc - c1.xc).powi(2) + (c2.yc - c1.yc).powi(2)).sqrt()
}

/// Check whether configuration q lies on circle c.
pub fn configuration_on_hc_cc_circle(c: &HcCcCircle, q: &Configuration) -> bool {
    let distance = point_distance(c.xc, c.yc, q.x, q.y);
    if (distance - c.radius).abs() > get_epsilon() {
        return false;
    }
    let mut angle = (q.y - c.yc).atan2(q.x - c.xc);
    if c.left && c.forward {
        angle = angle + HALF_PI - c.mu;
    } else if c.left && !c.forward {
        angle = angle + HALF_PI + c.mu;
    } else if !c.left && c.forward {
        angle = angle - HALF_PI + c.mu;
    } else {
        angle = angle - HALF_PI - c.mu;
    }
    let angle = twopify(angle);
    (q.theta - angle).abs() < get_epsilon()
}
