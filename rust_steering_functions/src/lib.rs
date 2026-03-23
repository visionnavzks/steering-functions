pub mod dubins_state_space;
pub mod hc_cc_state_space;
pub mod reeds_shepp_state_space;
pub mod state;
pub mod utilities;

pub use dubins_state_space::DubinsStateSpace;
pub use hc_cc_state_space::HcCcStateSpace;
pub use reeds_shepp_state_space::ReedsSheppStateSpace;
pub use state::{Control, State};
pub use utilities::*;

pub trait StateSpace {
    fn get_path(&self, state1: &State, state2: &State) -> Vec<State>;
    fn get_controls(&self, state1: &State, state2: &State) -> Vec<Control>;
    fn get_all_controls(&self, state1: &State, state2: &State) -> Vec<Vec<Control>>;
    fn integrate(&self, state: &State, controls: &[Control]) -> Vec<State>;
    fn interpolate(&self, state: &State, controls: &[Control], t: f64) -> State;
}

/// Helper struct for steering paths
pub struct SteeringFunctions;

impl SteeringFunctions {
    pub fn get_dubins_path(
        state1: &State,
        state2: &State,
        kappa: f64,
        discretization: f64,
        forwards: bool,
    ) -> Vec<State> {
        let space = DubinsStateSpace::new(kappa, discretization, forwards);
        space.get_path(state1, state2)
    }

    pub fn get_reeds_shepp_path(
        state1: &State,
        state2: &State,
        kappa: f64,
        discretization: f64,
    ) -> Vec<State> {
        let space = ReedsSheppStateSpace::new(kappa, discretization);
        space.get_path(state1, state2)
    }

    pub fn get_hc_cc_path(
        state1: &State,
        state2: &State,
        kappa: f64,
        sigma: f64,
        discretization: f64,
    ) -> Vec<State> {
        let space = HcCcStateSpace::new(kappa, sigma, discretization);
        space.get_path(state1, state2)
    }
}
