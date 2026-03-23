use crate::steering_functions::{Control, State, StateWithCovariance};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Configuration {
    pub x: f64,
    pub y: f64,
    pub theta: f64,
    pub kappa: f64,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Path {
    pub controls: Vec<Control>,
}

pub trait HcCcStateSpace {
    fn get_controls(&self, state1: &State, state2: &State) -> Vec<Control>;

    fn integrate(&self, state: &State, controls: &[Control], discretization: f64) -> Vec<State>;

    fn integrate_with_covariance(
        &self,
        state: &State,
        controls: &[Control],
        discretization: f64,
    ) -> Vec<StateWithCovariance>;
}

macro_rules! define_hc_cc_space {
    ($name:ident) => {
        #[derive(Debug, Clone, Copy)]
        pub struct $name;

        impl HcCcStateSpace for $name {
            fn get_controls(&self, _state1: &State, _state2: &State) -> Vec<Control> {
                Vec::new()
            }

            fn integrate(
                &self,
                state: &State,
                _controls: &[Control],
                _discretization: f64,
            ) -> Vec<State> {
                vec![*state]
            }

            fn integrate_with_covariance(
                &self,
                state: &State,
                _controls: &[Control],
                _discretization: f64,
            ) -> Vec<StateWithCovariance> {
                vec![StateWithCovariance {
                    state: *state,
                    ..StateWithCovariance::default()
                }]
            }
        }
    };
}

define_hc_cc_space!(Cc00DubinsStateSpace);
define_hc_cc_space!(Cc0pmDubinsStateSpace);
define_hc_cc_space!(Ccpm0DubinsStateSpace);
define_hc_cc_space!(CcpmpmDubinsStateSpace);
define_hc_cc_space!(CcDubinsStateSpace);

define_hc_cc_space!(HcReedsSheppStateSpace);
define_hc_cc_space!(Hc00ReedsSheppStateSpace);
define_hc_cc_space!(Hc0pmReedsSheppStateSpace);
define_hc_cc_space!(Hcpm0ReedsSheppStateSpace);
define_hc_cc_space!(HcpmpmReedsSheppStateSpace);
define_hc_cc_space!(Cc00ReedsSheppStateSpace);
