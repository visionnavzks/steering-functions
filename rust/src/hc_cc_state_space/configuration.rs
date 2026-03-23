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

use crate::utilities::{get_epsilon, point_distance, twopify};

/// A configuration in the plane with position, orientation, and curvature.
///
/// Orientation `theta` is always stored in `[0, 2π)`.
#[derive(Debug, Clone, Copy)]
pub struct Configuration {
    /// Position x.
    pub x: f64,
    /// Position y.
    pub y: f64,
    /// Orientation in radians, normalised to `[0, 2π)`.
    pub theta: f64,
    /// Curvature.
    pub kappa: f64,
}

impl Configuration {
    /// Create a new configuration. `theta` is normalised to `[0, 2π)`.
    pub fn new(x: f64, y: f64, theta: f64, kappa: f64) -> Self {
        Self {
            x,
            y,
            theta: twopify(theta),
            kappa,
        }
    }
}

impl Default for Configuration {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            theta: 0.0,
            kappa: 0.0,
        }
    }
}

impl fmt::Display for Configuration {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {}, {}, {})", self.x, self.y, self.theta, self.kappa)
    }
}

/// Cartesian distance between two configurations.
pub fn configuration_distance(q1: &Configuration, q2: &Configuration) -> f64 {
    point_distance(q1.x, q1.y, q2.x, q2.y)
}

/// Test whether two configurations are aligned (same orientation and collinear).
pub fn configuration_aligned(q1: &Configuration, q2: &Configuration) -> bool {
    if (q2.theta - q1.theta).abs() > get_epsilon() {
        return false;
    }
    let angle = twopify((q2.y - q1.y).atan2(q2.x - q1.x));
    (angle - q1.theta).abs() <= get_epsilon()
}

/// Test whether two configurations are equal (same position and orientation
/// within [`EPSILON`](crate::utilities::EPSILON)).
pub fn configuration_equal(q1: &Configuration, q2: &Configuration) -> bool {
    if (q2.theta - q1.theta).abs() > get_epsilon() {
        return false;
    }
    if configuration_distance(q1, q2) > get_epsilon() {
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utilities::{EPSILON, PI, TWO_PI};

    #[test]
    fn test_configuration_default() {
        let c = Configuration::default();
        assert_eq!(c.x, 0.0);
        assert_eq!(c.y, 0.0);
        assert_eq!(c.theta, 0.0);
        assert_eq!(c.kappa, 0.0);
    }

    #[test]
    fn test_configuration_new_normalises_theta() {
        let c = Configuration::new(1.0, 2.0, -1.0, 0.5);
        assert!(c.theta >= 0.0 && c.theta < TWO_PI);
    }

    #[test]
    fn test_configuration_new_large_theta() {
        let c = Configuration::new(0.0, 0.0, 10.0, 0.0);
        assert!(c.theta >= 0.0 && c.theta < TWO_PI);
    }

    #[test]
    fn test_configuration_display() {
        let c = Configuration::new(1.0, 2.0, 0.0, 0.5);
        let s = format!("{}", c);
        assert!(s.starts_with('('));
        assert!(s.contains("1"));
    }

    #[test]
    fn test_configuration_distance() {
        let q1 = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let q2 = Configuration::new(3.0, 4.0, 0.0, 0.0);
        assert!((configuration_distance(&q1, &q2) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_configuration_distance_same_point() {
        let q = Configuration::new(1.0, 2.0, PI, 0.0);
        assert!(configuration_distance(&q, &q) < 1e-12);
    }

    #[test]
    fn test_configuration_equal_identical() {
        let q1 = Configuration::new(1.0, 2.0, 0.5, 0.1);
        let q2 = Configuration::new(1.0, 2.0, 0.5, 0.1);
        assert!(configuration_equal(&q1, &q2));
    }

    #[test]
    fn test_configuration_equal_different_theta() {
        let q1 = Configuration::new(1.0, 2.0, 0.5, 0.1);
        let q2 = Configuration::new(1.0, 2.0, 1.0, 0.1);
        assert!(!configuration_equal(&q1, &q2));
    }

    #[test]
    fn test_configuration_equal_different_position() {
        let q1 = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let q2 = Configuration::new(1.0, 0.0, 0.0, 0.0);
        assert!(!configuration_equal(&q1, &q2));
    }

    #[test]
    fn test_configuration_aligned_collinear() {
        let q1 = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let q2 = Configuration::new(1.0, 0.0, 0.0, 0.0);
        assert!(configuration_aligned(&q1, &q2));
    }

    #[test]
    fn test_configuration_aligned_different_theta() {
        let q1 = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let q2 = Configuration::new(1.0, 0.0, 1.0, 0.0);
        assert!(!configuration_aligned(&q1, &q2));
    }

    #[test]
    fn test_configuration_aligned_not_collinear() {
        let q1 = Configuration::new(0.0, 0.0, 0.0, 0.0);
        let q2 = Configuration::new(0.0, 1.0, 0.0, 0.0);
        assert!(!configuration_aligned(&q1, &q2));
    }

    #[test]
    fn test_configuration_equal_within_epsilon() {
        let q1 = Configuration::new(1.0, 2.0, 0.5, 0.0);
        let q2 = Configuration::new(1.0 + EPSILON * 0.5, 2.0, 0.5, 0.0);
        assert!(configuration_equal(&q1, &q2));
    }
}
