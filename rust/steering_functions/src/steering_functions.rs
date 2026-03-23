#[derive(Debug, Clone, Copy, PartialEq)]
pub struct State {
    pub x: f64,
    pub y: f64,
    pub theta: f64,
    pub kappa: f64,
    pub d: i32,
}

impl Default for State {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            theta: 0.0,
            kappa: 0.0,
            d: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StateWithCovariance {
    pub state: State,
    pub sigma: [f64; 16],
    pub lambda: [f64; 16],
    pub covariance: [f64; 16],
}

impl Default for StateWithCovariance {
    fn default() -> Self {
        Self {
            state: State::default(),
            sigma: [0.0; 16],
            lambda: [0.0; 16],
            covariance: [0.0; 16],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Control {
    pub delta_s: f64,
    pub kappa: f64,
    pub sigma: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MotionNoise {
    pub alpha1: f64,
    pub alpha2: f64,
    pub alpha3: f64,
    pub alpha4: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MeasurementNoise {
    pub std_x: f64,
    pub std_y: f64,
    pub std_theta: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Controller {
    pub k1: f64,
    pub k2: f64,
    pub k3: f64,
}

pub trait BaseStateSpace {
    fn get_path(&self, state1: &State, state2: &State, controls: &mut Vec<Control>) -> Vec<State>;

    fn interpolate(&self, state: &State, controls: &[Control], t: f64) -> State;

    fn integrate(&self, state: &State, controls: &[Control], discretization: f64) -> Vec<State>;

    fn get_all_controls(&self, state1: &State, state2: &State) -> Vec<Vec<Control>>;
}
