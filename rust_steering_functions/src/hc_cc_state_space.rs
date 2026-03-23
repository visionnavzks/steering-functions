use crate::utilities::{
    end_of_circular_arc, end_of_clothoid, end_of_straight_line, point_distance, sgn,
};
use crate::{Control, State, StateSpace};
use std::f64;

#[derive(Debug, Clone, Copy, Default)]
pub struct Configuration {
    pub x: f64,
    pub y: f64,
    pub theta: f64,
    pub kappa: f64,
}

impl Configuration {
    pub fn new(x: f64, y: f64, theta: f64, kappa: f64) -> Self {
        Self { x, y, theta, kappa }
    }
}

#[derive(Debug, Clone, Copy, Default)]
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
    pub fn new(
        kappa: f64,
        sigma: f64,
        radius: f64,
        mu: f64,
        sin_mu: f64,
        cos_mu: f64,
        delta_min: f64,
    ) -> Self {
        Self {
            kappa,
            kappa_inv: 1.0 / kappa,
            sigma,
            radius,
            mu,
            sin_mu,
            cos_mu,
            delta_min,
        }
    }
}

pub struct HcCcStateSpace {
    pub kappa: f64,
    pub sigma: f64,
    pub discretization: f64,
    pub hc_cc_circle_param: HcCcCircleParam,
}

impl HcCcStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64) -> Self {
        assert!(kappa > 0.0 && sigma > 0.0 && discretization > 0.0);

        let length_min = kappa / sigma;
        let mut x_i = 0.0;
        let mut y_i = 0.0;
        let mut theta_i = 0.0;
        let mut kappa_i = 0.0;

        if length_min > crate::EPSILON {
            let (xf, yf, tf, kf) = end_of_clothoid(0.0, 0.0, 0.0, 0.0, sigma, 1.0, length_min);
            x_i = xf;
            y_i = yf;
            theta_i = tf;
            kappa_i = kf;
        }

        let xc = x_i - theta_i.sin() / kappa;
        let yc = y_i + theta_i.cos() / kappa;
        let radius = point_distance(xc, yc, 0.0, 0.0);

        let mu = (xc / yc).abs().atan();
        let sin_mu = mu.sin();
        let cos_mu = mu.cos();

        let delta_min = 0.5 * kappa.powi(2) / sigma;

        Self {
            kappa,
            sigma,
            discretization,
            hc_cc_circle_param: HcCcCircleParam::new(
                kappa, sigma, radius, mu, sin_mu, cos_mu, delta_min,
            ),
        }
    }

    #[inline]
    pub fn integrate_ode(&self, state: &State, control: &Control, integration_step: f64) -> State {
        let mut state_next = *state;
        let kappa = control.kappa;
        let sigma = control.sigma;
        let d = sgn(control.delta_s);

        if sigma.abs() > crate::EPSILON {
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
        } else {
            if kappa.abs() > crate::EPSILON {
                let (x_f, y_f, theta_f) =
                    end_of_circular_arc(state.x, state.y, state.theta, kappa, d, integration_step);
                state_next.x = x_f;
                state_next.y = y_f;
                state_next.theta = theta_f;
                state_next.kappa = kappa;
            } else {
                let (x_f, y_f) =
                    end_of_straight_line(state.x, state.y, state.theta, d, integration_step);
                state_next.x = x_f;
                state_next.y = y_f;
                state_next.theta = state.theta;
            }
        }
        state_next.d = d;
        state_next.s = state.s + integration_step;
        state_next
    }
}

impl StateSpace for HcCcStateSpace {
    fn get_controls(&self, _state1: &State, _state2: &State) -> Vec<Control> {
        // This is a placeholder as the CC/HC implementation involves many turn types similar to RS.
        // It's quite long to fully translate CC/HC paths generation from hc_cc_state_space without
        // access to `cc_steer.cpp` which often contains the actual get_controls logic in CC.
        // Since we are translating what's available, we will return empty or throw unimplemented if it requires it.
        // Let's return empty for now since the original C++ header declares `get_controls` as purely virtual (`= 0`)
        // in `HC_CC_State_Space`, meaning HC/CC specific classes derived from this implement it.
        unimplemented!("get_controls is pure virtual in HC_CC_State_Space and depends on CC/HC specific steer implementations")
    }

    fn get_all_controls(&self, _state1: &State, _state2: &State) -> Vec<Vec<Control>> {
        unimplemented!("get_all_controls is pure virtual in HC_CC_State_Space and depends on CC/HC specific steer implementations")
    }

    fn get_path(&self, state1: &State, state2: &State) -> Vec<State> {
        let un_filtered_controls = self.get_controls(state1, state2);
        let mut controls = Vec::new();
        for control in un_filtered_controls {
            if control.delta_s.abs() > 1e-6 {
                controls.push(control);
            }
        }
        self.integrate(state1, &controls)
    }

    fn integrate(&self, state: &State, controls: &[Control]) -> Vec<State> {
        let mut path = Vec::new();
        if controls.is_empty() {
            return path;
        }

        let mut n_states = 0;
        for control in controls {
            n_states += (control.delta_s.abs() / self.discretization).ceil() as usize;
        }
        path.reserve(n_states + 13);

        let mut state_curr = *state;
        state_curr.kappa = controls[0].kappa;
        state_curr.sigma = controls[0].sigma;
        state_curr.d = sgn(controls[0].delta_s);
        path.push(state_curr);

        for control in controls {
            let abs_delta_s = control.delta_s.abs();
            let n = f64::max(2.0, (abs_delta_s / self.discretization).ceil()) as usize;
            let discretization = abs_delta_s / (n as f64);

            for _ in 0..n {
                let state_next = self.integrate_ode(&state_curr, control, discretization);
                path.push(state_next);
                state_curr = state_next;
            }
        }
        path
    }

    fn interpolate(&self, state: &State, controls: &[Control], mut t: f64) -> State {
        let mut state_curr = *state;

        let mut s_path = 0.0;
        for control in controls {
            s_path += control.delta_s.abs();
        }
        t = t.clamp(0.0, 1.0);
        let s_inter = t * s_path;

        let mut s = 0.0;
        for control in controls {
            let abs_delta_s = control.delta_s.abs();

            if s_inter - s > abs_delta_s {
                state_curr = self.integrate_ode(&state_curr, control, abs_delta_s);
                s += abs_delta_s;
            } else {
                return self.integrate_ode(&state_curr, control, s_inter - s);
            }
        }
        state_curr
    }
}
