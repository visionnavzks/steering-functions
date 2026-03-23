use crate::steering_functions::{BaseStateSpace, Control, State};

#[derive(Debug, Clone, Copy)]
pub struct DubinsStateSpace {
    pub kappa: f64,
    pub discretization: f64,
}

impl DubinsStateSpace {
    pub fn new(kappa: f64, discretization: f64) -> Self {
        Self {
            kappa,
            discretization,
        }
    }
}

impl BaseStateSpace for DubinsStateSpace {
    fn get_path(&self, state1: &State, state2: &State, controls: &mut Vec<Control>) -> Vec<State> {
        controls.clear();
        controls.push(Control {
            delta_s: ((state2.x - state1.x).powi(2) + (state2.y - state1.y).powi(2)).sqrt(),
            kappa: 0.0,
            sigma: 0.0,
        });
        vec![*state1, *state2]
    }

    fn interpolate(&self, state: &State, controls: &[Control], t: f64) -> State {
        if controls.is_empty() {
            return *state;
        }
        let clamped_t = t.clamp(0.0, 1.0);
        State {
            x: state.x + controls[0].delta_s * clamped_t * state.theta.cos(),
            y: state.y + controls[0].delta_s * clamped_t * state.theta.sin(),
            theta: state.theta,
            kappa: state.kappa,
            d: state.d,
        }
    }

    fn integrate(&self, state: &State, controls: &[Control], _discretization: f64) -> Vec<State> {
        if controls.is_empty() {
            return vec![*state];
        }
        vec![*state, self.interpolate(state, controls, 1.0)]
    }

    fn get_all_controls(&self, state1: &State, state2: &State) -> Vec<Vec<Control>> {
        let mut controls = Vec::new();
        let _ = self.get_path(state1, state2, &mut controls);
        vec![controls]
    }
}
