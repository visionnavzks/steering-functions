/// Kinematic car state.
#[derive(Clone, Copy, Debug, Default)]
pub struct State {
    pub x: f64,
    pub y: f64,
    pub theta: f64,
    pub kappa: f64,
    pub sigma: f64,
    pub d: f64,
    pub s: f64,
    pub vel: f64,
    pub acc: f64,
    pub time: f64,
    pub fork_y: f64,
}

impl State {
    pub fn nearly_equal(&self, other: &State) -> bool {
        const EPS: f64 = 1e-6;
        (self.x - other.x).abs() < EPS
            && (self.y - other.y).abs() < EPS
            && (self.theta - other.theta).abs() < EPS
            && (self.kappa - other.kappa).abs() < EPS
            && (self.sigma - other.sigma).abs() < EPS
            && (self.d - other.d).abs() < EPS
            && (self.s - other.s).abs() < EPS
            && (self.vel - other.vel).abs() < EPS
            && (self.acc - other.acc).abs() < EPS
            && (self.time - other.time).abs() < EPS
            && (self.fork_y - other.fork_y).abs() < EPS
    }
}

impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        self.nearly_equal(other)
    }
}

/// Path segment control input.
#[derive(Clone, Copy, Debug, Default)]
pub struct Control {
    pub delta_s: f64,
    pub kappa: f64,
    pub sigma: f64,
}
