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

use crate::state::{Control, State};
use crate::utilities::{
    end_of_circular_arc, end_of_clothoid, end_of_straight_line, get_epsilon, sgn,
};

/// Integrate the kinematic ODE for a single control step.
///
/// Given a state and a control input, computes the next state after traveling
/// `integration_step` distance. The control's `sigma`, `kappa`, and `delta_s`
/// determine whether a clothoid, circular arc, or straight line segment is used.
pub fn integrate_ode(state: &State, control: &Control, integration_step: f64) -> State {
    let kappa = control.kappa;
    let sigma = control.sigma;
    let d = sgn(control.delta_s);

    let mut state_next = State::default();

    if sigma.abs() > get_epsilon() {
        // Clothoid segment
        let (x_f, y_f, theta_f, kappa_f) = end_of_clothoid(
            state.x,
            state.y,
            state.theta,
            state.kappa,
            sigma,
            d,
            integration_step,
        );
        state_next.x = x_f;
        state_next.y = y_f;
        state_next.theta = theta_f;
        state_next.kappa = kappa_f;
        state_next.sigma = sigma;
    } else if kappa.abs() > get_epsilon() {
        // Circular arc segment
        let (x_f, y_f, theta_f) =
            end_of_circular_arc(state.x, state.y, state.theta, kappa, d, integration_step);
        state_next.x = x_f;
        state_next.y = y_f;
        state_next.theta = theta_f;
        state_next.kappa = kappa;
    } else {
        // Straight line segment
        let (x_f, y_f) =
            end_of_straight_line(state.x, state.y, state.theta, d, integration_step);
        state_next.x = x_f;
        state_next.y = y_f;
        state_next.theta = state.theta;
    }

    state_next.d = d;
    state_next.s = state.s + integration_step;
    state_next
}

/// Trait for state-space implementations that compute controls and paths.
///
/// Implementors must provide [`get_controls`], [`get_all_controls`], and
/// [`discretization`].  Default implementations of path-related helpers are
/// built on top of these three methods together with [`integrate_ode`].
pub trait StateSpace {
    /// Discretization step size used for path integration.
    fn discretization(&self) -> f64;

    /// Compute the shortest control sequence to steer from `state1` to `state2`.
    fn get_controls(&self, state1: &State, state2: &State) -> Vec<Control>;

    /// Compute *all* feasible control sequences from `state1` to `state2`.
    fn get_all_controls(&self, state1: &State, state2: &State) -> Vec<Vec<Control>>;

    /// Compute the shortest path from `state1` to `state2`.
    fn get_path(&self, state1: &State, state2: &State) -> Vec<State> {
        let controls = self.get_controls(state1, state2);
        self.integrate(state1, &controls, self.discretization())
    }

    /// Compute the shortest path and return it together with the controls used.
    fn get_path_with_controls(&self, state1: &State, state2: &State) -> (Vec<State>, Vec<Control>) {
        let controls = self.get_controls(state1, state2);
        let path = self.integrate(state1, &controls, self.discretization());
        (path, controls)
    }

    /// Interpolate along a control sequence at parameter `t` ∈ \[0, 1\].
    ///
    /// Returns the state at arc-length fraction `t` of the total path defined
    /// by `controls`, starting from `state`.
    fn interpolate(&self, state: &State, controls: &[Control], t: f64) -> State {
        let state_inter = *state;
        if controls.is_empty() {
            return state_inter;
        }

        let s_path: f64 = controls.iter().map(|c| c.delta_s.abs()).sum();
        let t = t.clamp(0.0, 1.0);
        let s_inter = t * s_path;

        let mut s_accumulated = 0.0;
        let mut state_curr = *state;

        for control in controls {
            let abs_delta_s = control.delta_s.abs();
            if s_inter - s_accumulated > abs_delta_s {
                state_curr = integrate_ode(&state_curr, control, abs_delta_s);
                s_accumulated += abs_delta_s;
            } else {
                return integrate_ode(&state_curr, control, s_inter - s_accumulated);
            }
        }

        state_inter
    }

    /// Integrate a control sequence into a discretized path.
    ///
    /// The path starts at `state` and follows each control in order, placing
    /// intermediate states at intervals of approximately `discretization`.
    fn integrate(&self, state: &State, controls: &[Control], discretization: f64) -> Vec<State> {
        if controls.is_empty() {
            return Vec::new();
        }

        // Pre-allocate with estimated capacity
        let mut n_states: usize = 0;
        for control in controls {
            let abs_delta_s = control.delta_s.abs();
            n_states += (abs_delta_s / discretization).ceil() as usize;
        }

        let mut path = Vec::with_capacity(n_states + 13);

        let mut state_curr = State {
            x: state.x,
            y: state.y,
            theta: state.theta,
            kappa: controls[0].kappa,
            sigma: controls[0].sigma,
            d: sgn(controls[0].delta_s),
            ..State::default()
        };
        path.push(state_curr);

        for control in controls {
            let abs_delta_s = control.delta_s.abs();
            let n = 2.0f64.max((abs_delta_s / discretization).ceil()) as usize;
            let step_size = abs_delta_s / n as f64;
            assert!(step_size > 0.0);

            let mut s_seg = 0.0;
            for _ in 0..n {
                s_seg += step_size;
                let integration_step = if s_seg > abs_delta_s {
                    let overshoot = s_seg - abs_delta_s;
                    s_seg = abs_delta_s;
                    step_size - overshoot
                } else {
                    step_size
                };

                let state_next = integrate_ode(&state_curr, control, integration_step);
                path.push(state_next);
                state_curr = state_next;
            }
        }

        path
    }

    /// Integrate a control sequence using the default discretization.
    fn integrate_default(&self, state: &State, controls: &[Control]) -> Vec<State> {
        self.integrate(state, controls, self.discretization())
    }

    /// Compute all feasible paths from `state1` to `state2`.
    fn get_all_paths(&self, state1: &State, state2: &State) -> Vec<Vec<State>> {
        let all_controls = self.get_all_controls(state1, state2);
        let disc = self.discretization();
        all_controls
            .iter()
            .map(|controls| self.integrate(state1, controls, disc))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utilities::HALF_PI;

    const TEST_EPS: f64 = 1e-6;

    // A minimal state space for testing the default trait methods.
    struct TestStateSpace {
        discretization: f64,
    }

    impl StateSpace for TestStateSpace {
        fn discretization(&self) -> f64 {
            self.discretization
        }

        fn get_controls(&self, _state1: &State, _state2: &State) -> Vec<Control> {
            // Drive straight for 1 unit
            vec![Control {
                delta_s: 1.0,
                kappa: 0.0,
                sigma: 0.0,
            }]
        }

        fn get_all_controls(&self, state1: &State, state2: &State) -> Vec<Vec<Control>> {
            vec![self.get_controls(state1, state2)]
        }
    }

    #[test]
    fn test_integrate_ode_straight_line() {
        let state = State::new(0.0, 0.0, 0.0, 0.0);
        let control = Control {
            delta_s: 1.0,
            kappa: 0.0,
            sigma: 0.0,
        };
        let next = integrate_ode(&state, &control, 1.0);
        assert!((next.x - 1.0).abs() < TEST_EPS);
        assert!(next.y.abs() < TEST_EPS);
        assert!(next.theta.abs() < TEST_EPS);
        assert!((next.d - 1.0).abs() < TEST_EPS);
        assert!((next.s - 1.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_integrate_ode_circular_arc() {
        let state = State::new(0.0, 0.0, 0.0, 0.0);
        let control = Control {
            delta_s: HALF_PI,
            kappa: 1.0,
            sigma: 0.0,
        };
        // Quarter circle of radius 1: should end at (1, 1) heading π/2
        let next = integrate_ode(&state, &control, HALF_PI);
        assert!((next.x - 1.0).abs() < 1e-4);
        assert!((next.y - 1.0).abs() < 1e-4);
        assert!((next.theta - HALF_PI).abs() < 1e-4);
    }

    #[test]
    fn test_integrate_ode_reverse_straight() {
        let state = State::new(1.0, 0.0, 0.0, 0.0);
        let control = Control {
            delta_s: -1.0,
            kappa: 0.0,
            sigma: 0.0,
        };
        let next = integrate_ode(&state, &control, 1.0);
        assert!((next.x).abs() < TEST_EPS);
        assert!(next.y.abs() < TEST_EPS);
        assert!((next.d - (-1.0)).abs() < TEST_EPS);
    }

    #[test]
    fn test_interpolate_at_boundaries() {
        let ss = TestStateSpace {
            discretization: 0.1,
        };
        let state = State::new(0.0, 0.0, 0.0, 0.0);
        let controls = vec![Control {
            delta_s: 2.0,
            kappa: 0.0,
            sigma: 0.0,
        }];

        // t = 0 should return the start
        let s0 = ss.interpolate(&state, &controls, 0.0);
        assert!(s0.x.abs() < TEST_EPS);

        // t = 1 should return the end
        let s1 = ss.interpolate(&state, &controls, 1.0);
        assert!((s1.x - 2.0).abs() < TEST_EPS);

        // t = 0.5 should be halfway
        let s_mid = ss.interpolate(&state, &controls, 0.5);
        assert!((s_mid.x - 1.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_interpolate_empty_controls() {
        let ss = TestStateSpace {
            discretization: 0.1,
        };
        let state = State::new(5.0, 3.0, 1.0, 0.0);
        let s = ss.interpolate(&state, &[], 0.5);
        assert!((s.x - 5.0).abs() < TEST_EPS);
        assert!((s.y - 3.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_integrate_produces_path() {
        let ss = TestStateSpace {
            discretization: 0.1,
        };
        let state = State::new(0.0, 0.0, 0.0, 0.0);
        let controls = vec![Control {
            delta_s: 1.0,
            kappa: 0.0,
            sigma: 0.0,
        }];

        let path = ss.integrate(&state, &controls, 0.5);
        // With step size 0.5 and length 1.0, n = max(2, ceil(1/0.5)) = 2
        // So we get initial state + 2 integrated states = 3
        assert_eq!(path.len(), 3);
        assert!(path[0].x.abs() < TEST_EPS);
        assert!((path.last().unwrap().x - 1.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_integrate_empty_controls() {
        let ss = TestStateSpace {
            discretization: 0.1,
        };
        let state = State::new(0.0, 0.0, 0.0, 0.0);
        let path = ss.integrate(&state, &[], 0.1);
        assert!(path.is_empty());
    }

    #[test]
    fn test_get_path() {
        let ss = TestStateSpace {
            discretization: 0.1,
        };
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(1.0, 0.0, 0.0, 0.0);
        let path = ss.get_path(&s1, &s2);
        assert!(!path.is_empty());
        assert!((path.last().unwrap().x - 1.0).abs() < TEST_EPS);
    }

    #[test]
    fn test_get_all_paths() {
        let ss = TestStateSpace {
            discretization: 0.1,
        };
        let s1 = State::new(0.0, 0.0, 0.0, 0.0);
        let s2 = State::new(1.0, 0.0, 0.0, 0.0);
        let all_paths = ss.get_all_paths(&s1, &s2);
        assert_eq!(all_paths.len(), 1);
        assert!(!all_paths[0].is_empty());
    }

    #[test]
    fn test_interpolate_multi_segment() {
        let ss = TestStateSpace {
            discretization: 0.1,
        };
        let state = State::new(0.0, 0.0, 0.0, 0.0);
        let controls = vec![
            Control {
                delta_s: 1.0,
                kappa: 0.0,
                sigma: 0.0,
            },
            Control {
                delta_s: 1.0,
                kappa: 0.0,
                sigma: 0.0,
            },
        ];

        // t=0.5 is halfway through total path length of 2.0, so x≈1.0
        let s_mid = ss.interpolate(&state, &controls, 0.5);
        assert!((s_mid.x - 1.0).abs() < TEST_EPS);

        // t=1.0 should reach the end at x≈2.0
        let s_end = ss.interpolate(&state, &controls, 1.0);
        assert!((s_end.x - 2.0).abs() < TEST_EPS);
    }
}
