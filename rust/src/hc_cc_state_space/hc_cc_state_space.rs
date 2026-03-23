// Copyright (c) 2017 Robert Bosch GmbH.
// All rights reserved.
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

use super::hc_cc_circle::HcCcCircleParam;
use crate::base_state_space::StateSpace;
use crate::state::{Control, State};
use crate::utilities::{end_of_clothoid, get_epsilon, point_distance};

/// Base data for HC/CC state spaces.
///
/// This struct holds the common parameters (curvature, sharpness,
/// discretisation step, and the derived circle parameters) shared by all
/// HC/CC steering variants.  Concrete state spaces (CC-Dubins, HC-RS,
/// CC-RS, …) embed an `HcCcStateSpace` and implement the [`StateSpace`]
/// trait by delegating to its helper methods.
#[derive(Debug, Clone)]
pub struct HcCcStateSpace {
    /// Maximum curvature.
    pub kappa: f64,
    /// Maximum sharpness.
    pub sigma: f64,
    /// Discretisation step for path integration.
    pub discretization: f64,
    /// Derived HC/CC circle parameters.
    pub hc_cc_circle_param: HcCcCircleParam,
}

impl HcCcStateSpace {
    /// Create a new `HcCcStateSpace`.
    ///
    /// # Panics
    ///
    /// Panics if `kappa`, `sigma`, or `discretization` is not positive.
    pub fn new(kappa: f64, sigma: f64, discretization: f64) -> Self {
        assert!(
            kappa > 0.0 && sigma > 0.0 && discretization > 0.0,
            "kappa, sigma and discretization must be positive"
        );

        // Intermediate configuration after the first clothoid.
        let length_min = kappa / sigma;
        let (x_i, y_i, theta_i) = if length_min > get_epsilon() {
            let (x, y, th, _kappa) = end_of_clothoid(0.0, 0.0, 0.0, 0.0, sigma, 1.0, length_min);
            (x, y, th)
        } else {
            (0.0, 0.0, 0.0)
        };

        // Centre of the circle tangent to the clothoid end-point.
        let xc = x_i - theta_i.sin() / kappa;
        let yc = y_i + theta_i.cos() / kappa;
        let radius = point_distance(xc, yc, 0.0, 0.0);

        let mu = (xc / yc).abs().atan();
        let sin_mu = mu.sin();
        let cos_mu = mu.cos();
        let delta_min = 0.5 * kappa * kappa / sigma;

        let mut hc_cc_circle_param = HcCcCircleParam::default();
        hc_cc_circle_param.set_param(kappa, sigma, radius, mu, sin_mu, cos_mu, delta_min);

        Self {
            kappa,
            sigma,
            discretization,
            hc_cc_circle_param,
        }
    }

    /// Filter negligible control segments and integrate.
    ///
    /// This is the HC/CC override of `get_path`: it obtains unfiltered
    /// controls via `raw_controls`, removes segments with
    /// `|delta_s| ≤ 1e-6`, and integrates the remainder.
    pub fn get_path_filtered<F>(
        &self,
        state1: &State,
        state2: &State,
        raw_controls_fn: F,
    ) -> (Vec<State>, Vec<Control>)
    where
        F: FnOnce(&State, &State) -> Vec<Control>,
    {
        let unfiltered = raw_controls_fn(state1, state2);
        let controls: Vec<Control> = unfiltered
            .into_iter()
            .filter(|c| c.delta_s.abs() > 1e-6)
            .collect();
        let path = self.integrate_controls(state1, &controls);
        (path, controls)
    }

    /// Integrate a control sequence using the stored discretisation.
    ///
    /// This is a convenience wrapper so that embedders don't have to
    /// construct a full [`StateSpace`] implementor just to call `integrate`.
    fn integrate_controls(&self, state: &State, controls: &[Control]) -> Vec<State> {
        /// A thin adapter that exposes `self.discretization` via the trait.
        struct Adapter(f64);
        impl StateSpace for Adapter {
            fn discretization(&self) -> f64 {
                self.0
            }
            fn get_controls(&self, _: &State, _: &State) -> Vec<Control> {
                Vec::new()
            }
            fn get_all_controls(&self, _: &State, _: &State) -> Vec<Vec<Control>> {
                Vec::new()
            }
        }
        let adapter = Adapter(self.discretization);
        adapter.integrate(state, controls, self.discretization)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_basic() {
        let ss = HcCcStateSpace::new(1.0, 1.0, 0.1);
        assert!((ss.kappa - 1.0).abs() < 1e-12);
        assert!((ss.sigma - 1.0).abs() < 1e-12);
        assert!((ss.discretization - 0.1).abs() < 1e-12);
        assert!(ss.hc_cc_circle_param.radius > 0.0);
        assert!(ss.hc_cc_circle_param.delta_min > 0.0);
    }

    #[test]
    #[should_panic]
    fn test_new_zero_kappa() {
        HcCcStateSpace::new(0.0, 1.0, 0.1);
    }

    #[test]
    #[should_panic]
    fn test_new_negative_sigma() {
        HcCcStateSpace::new(1.0, -1.0, 0.1);
    }

    #[test]
    fn test_param_consistency() {
        let ss = HcCcStateSpace::new(2.0, 1.0, 0.1);
        let p = &ss.hc_cc_circle_param;
        assert!((p.kappa - 2.0).abs() < 1e-12);
        assert!((p.sigma - 1.0).abs() < 1e-12);
        assert!((p.kappa_inv - 0.5).abs() < 1e-12);
        assert!((p.delta_min - 2.0).abs() < 1e-12); // 0.5 * 4.0 / 1.0
    }

    #[test]
    fn test_get_path_filtered_removes_tiny() {
        let ss = HcCcStateSpace::new(1.0, 1.0, 0.1);
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(1.0, 0.0, 0.0, 0.0);
        let (path, controls) = ss.get_path_filtered(&s1, &s2, |_, _| {
            vec![
                Control {
                    delta_s: 1e-8,
                    kappa: 0.0,
                    sigma: 0.0,
                },
                Control {
                    delta_s: 1.0,
                    kappa: 0.0,
                    sigma: 0.0,
                },
            ]
        });
        // The tiny segment should be filtered out
        assert_eq!(controls.len(), 1);
        assert!(!path.is_empty());
    }
}
