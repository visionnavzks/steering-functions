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

use std::fmt;

/// Tolerance for floating-point equality comparisons.
const STATE_EPSILON: f64 = 1e-6;

/// Description of a kinematic car's state.
#[derive(Debug, Clone, Copy)]
pub struct State {
    /// Position in x of the robot.
    pub x: f64,
    /// Position in y of the robot.
    pub y: f64,
    /// Orientation of the robot.
    pub theta: f64,
    /// Curvature at position (x, y).
    pub kappa: f64,
    /// Sharpness (derivative of curvature w.r.t. arc length).
    pub sigma: f64,
    /// Driving direction {-1, 0, 1}.
    pub d: f64,
    /// Arc length.
    pub s: f64,
    /// Velocity.
    pub vel: f64,
    /// Acceleration.
    pub acc: f64,
    /// Time.
    pub time: f64,
    /// Fork y position.
    pub fork_y: f64,
}

impl State {
    /// Create a new state with position, orientation, and optional curvature.
    pub fn new(x: f64, y: f64, theta: f64, kappa: f64) -> Self {
        Self {
            x,
            y,
            theta,
            kappa,
            sigma: 0.0,
            d: 0.0,
            s: 0.0,
            vel: 0.0,
            acc: 0.0,
            time: 0.0,
            fork_y: 0.0,
        }
    }
}

impl Default for State {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            theta: 0.0,
            kappa: 0.0,
            sigma: 0.0,
            d: 0.0,
            s: 0.0,
            vel: 0.0,
            acc: 0.0,
            time: 0.0,
            fork_y: 0.0,
        }
    }
}

impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        (self.x - other.x).abs() < STATE_EPSILON
            && (self.y - other.y).abs() < STATE_EPSILON
            && (self.theta - other.theta).abs() < STATE_EPSILON
            && (self.kappa - other.kappa).abs() < STATE_EPSILON
            && (self.sigma - other.sigma).abs() < STATE_EPSILON
            && (self.d - other.d).abs() < STATE_EPSILON
            && (self.s - other.s).abs() < STATE_EPSILON
            && (self.vel - other.vel).abs() < STATE_EPSILON
            && (self.acc - other.acc).abs() < STATE_EPSILON
            && (self.time - other.time).abs() < STATE_EPSILON
            && (self.fork_y - other.fork_y).abs() < STATE_EPSILON
    }
}

impl fmt::Display for State {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "State(x: {}, y: {}, theta: {}, kappa: {}, sigma: {}, d: {}, s: {}, vel: {}, acc: {}, time: {}, fork_y: {})",
            self.x, self.y, self.theta, self.kappa, self.sigma, self.d, self.s, self.vel, self.acc,
            self.time, self.fork_y
        )
    }
}

/// Description of a path segment with its corresponding control inputs.
#[derive(Debug, Clone, Copy)]
pub struct Control {
    /// Signed arc length of a segment.
    pub delta_s: f64,
    /// Curvature at the beginning of a segment.
    pub kappa: f64,
    /// Sharpness (derivative of curvature w.r.t. arc length) of a segment.
    pub sigma: f64,
}

impl fmt::Display for Control {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Control Segment [delta_s: {}, curvature: {}, sharpness: {}]",
            self.delta_s, self.kappa, self.sigma
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_state_default() {
        let s = State::default();
        assert_eq!(s.x, 0.0);
        assert_eq!(s.y, 0.0);
        assert_eq!(s.theta, 0.0);
        assert_eq!(s.kappa, 0.0);
        assert_eq!(s.sigma, 0.0);
        assert_eq!(s.d, 0.0);
        assert_eq!(s.s, 0.0);
        assert_eq!(s.vel, 0.0);
        assert_eq!(s.acc, 0.0);
        assert_eq!(s.time, 0.0);
        assert_eq!(s.fork_y, 0.0);
    }

    #[test]
    fn test_state_new() {
        let s = State::new(1.0, 2.0, 3.0, 4.0);
        assert_eq!(s.x, 1.0);
        assert_eq!(s.y, 2.0);
        assert_eq!(s.theta, 3.0);
        assert_eq!(s.kappa, 4.0);
        assert_eq!(s.sigma, 0.0);
    }

    #[test]
    fn test_state_equality_exact() {
        let a = State::new(1.0, 2.0, 3.0, 4.0);
        let b = State::new(1.0, 2.0, 3.0, 4.0);
        assert_eq!(a, b);
    }

    #[test]
    fn test_state_equality_within_epsilon() {
        let a = State::new(1.0, 2.0, 3.0, 4.0);
        let b = State::new(1.0 + 5e-7, 2.0 - 5e-7, 3.0, 4.0);
        assert_eq!(a, b);
    }

    #[test]
    fn test_state_inequality() {
        let a = State::new(1.0, 2.0, 3.0, 4.0);
        let b = State::new(1.0 + 1e-5, 2.0, 3.0, 4.0);
        assert_ne!(a, b);
    }

    #[test]
    fn test_state_display() {
        let s = State::new(1.0, 2.0, 3.0, 4.0);
        let display = format!("{}", s);
        assert!(display.starts_with("State("));
        assert!(display.contains("x: 1"));
    }

    #[test]
    fn test_control_display() {
        let c = Control {
            delta_s: 1.5,
            kappa: 0.1,
            sigma: 0.01,
        };
        let display = format!("{}", c);
        assert!(display.starts_with("Control Segment"));
        assert!(display.contains("delta_s: 1.5"));
    }
}
