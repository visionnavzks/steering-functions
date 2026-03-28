use crate::state::{State, Control};
use crate::base_state_space::StateSpace;
use crate::hc_cc_state_space::HcCcStateSpaceParams;

/// CC-Dubins state space (zero start/end curvature variant).
pub struct CcDubinsStateSpace {
    params_: HcCcStateSpaceParams,
    discretization_: f64,
    forward_: bool,
}

impl CcDubinsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64, forward: bool) -> Self {
        Self {
            params_: HcCcStateSpaceParams::new(kappa, sigma),
            discretization_: discretization,
            forward_: forward,
        }
    }
}

impl StateSpace for CcDubinsStateSpace {
    fn get_controls(&self, _s1: &State, _s2: &State) -> Vec<Control> {
        Vec::new()
    }
    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        vec![self.get_controls(s1, s2)]
    }
    fn discretization(&self) -> f64 { self.discretization_ }
}
