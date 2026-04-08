use crate::state::{State, Control};
use crate::utilities::{
    sgn, get_epsilon,
    end_of_clothoid, end_of_circular_arc, end_of_straight_line,
};

/// Trait implemented by all steering state spaces.
pub trait StateSpace {
    /// Return the control sequence for the shortest path from `s1` to `s2`.
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control>;

    /// Return all possible control sequences from `s1` to `s2`.
    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>>;

    /// Discretisation step length.
    fn discretization(&self) -> f64;

    /// Integrate a single ODE step.
    fn integrate_ode(state: &State, control: &Control, integration_step: f64) -> State {
        let mut next = State::default();
        let kappa = control.kappa;
        let sigma = control.sigma;
        let d = sgn(control.delta_s);

        if sigma.abs() > get_epsilon() {
            let (xf, yf, tf, kf) = end_of_clothoid(
                state.x, state.y, state.theta, state.kappa,
                sigma, d, integration_step,
            );
            next.x = xf;
            next.y = yf;
            next.theta = tf;
            next.kappa = kf;
            next.sigma = sigma;
        } else if kappa.abs() > get_epsilon() {
            let (xf, yf, tf) = end_of_circular_arc(
                state.x, state.y, state.theta,
                kappa, d, integration_step,
            );
            next.x = xf;
            next.y = yf;
            next.theta = tf;
            next.kappa = kappa;
        } else {
            let (xf, yf) = end_of_straight_line(
                state.x, state.y, state.theta,
                d, integration_step,
            );
            next.x = xf;
            next.y = yf;
            next.theta = state.theta;
        }

        next.d = d;
        next.s = state.s + integration_step;
        next
    }

    /// Compute the discretised path from `state1` to `state2`.
    fn get_path(&self, state1: &State, state2: &State) -> Vec<State> {
        let controls = self.get_controls(state1, state2);
        self.integrate(state1, &controls)
    }

    /// Integrate controls producing a discretised path.
    fn integrate(&self, state: &State, controls: &[Control]) -> Vec<State> {
        self.integrate_with_disc(state, controls, self.discretization())
    }

    fn integrate_with_disc(&self, state: &State, controls: &[Control], disc: f64) -> Vec<State> {
        let mut path = Vec::new();
        if controls.is_empty() {
            return path;
        }

        let mut curr = State {
            x: state.x,
            y: state.y,
            theta: state.theta,
            kappa: controls[0].kappa,
            sigma: controls[0].sigma,
            d: sgn(controls[0].delta_s),
            s: state.s,
            ..State::default()
        };
        path.push(curr);

        for control in controls {
            let abs_ds = control.delta_s.abs();
            let n = ((abs_ds / disc).ceil() as usize).max(2);
            let step = abs_ds / n as f64;
            let mut s_seg = 0.0_f64;

            for _ in 0..n {
                s_seg += step;
                let integration_step = if s_seg > abs_ds {
                    step - (s_seg - abs_ds)
                } else {
                    step
                };
                let next = Self::integrate_ode(&curr, control, integration_step);
                path.push(next);
                curr = next;
            }
        }
        path
    }

    /// Interpolate at normalised parameter `t` ∈ [0, 1].
    fn interpolate(&self, state: &State, controls: &[Control], t: f64) -> State {
        let t = t.clamp(0.0, 1.0);
        if controls.is_empty() {
            return *state;
        }

        let s_path: f64 = controls.iter().map(|c| c.delta_s.abs()).sum();
        let s_inter = t * s_path;
        let mut s_accum = 0.0_f64;
        let mut curr = *state;
        let mut result = *state;

        for control in controls {
            let abs_ds = control.delta_s.abs();
            if s_inter - s_accum > abs_ds {
                curr = Self::integrate_ode(&curr, control, abs_ds);
                s_accum += abs_ds;
            } else {
                result = Self::integrate_ode(&curr, control, s_inter - s_accum);
                return result;
            }
        }
        result
    }

    /// All discretised paths between two states.
    fn get_all_paths(&self, s1: &State, s2: &State) -> Vec<Vec<State>> {
        self.get_all_controls(s1, s2)
            .iter()
            .map(|cs| self.integrate(s1, cs))
            .collect()
    }
}
