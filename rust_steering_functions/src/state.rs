#[derive(Debug, Clone, Copy, Default)]
pub struct State {
    pub x: f64,
    pub y: f64,
    pub theta: f64,
    pub kappa: f64,
    pub sigma: f64,
    pub d: f64, // Direction
    pub s: f64, // Path length
    // the following fields are omitted or optional, but since user provided them I include them
    pub vel: f64,
    pub acc: f64,
    pub time: f64,
    pub fork_y: f64,
}

impl State {
    pub fn new(x: f64, y: f64, theta: f64, kappa: f64) -> Self {
        Self {
            x,
            y,
            theta,
            kappa,
            ..Default::default()
        }
    }

    pub fn eq_with_tol(&self, other: &Self, tol: f64) -> bool {
        (self.x - other.x).abs() < tol
            && (self.y - other.y).abs() < tol
            && (self.theta - other.theta).abs() < tol
            && (self.kappa - other.kappa).abs() < tol
            && (self.sigma - other.sigma).abs() < tol
            && (self.d - other.d).abs() < tol
            && (self.s - other.s).abs() < tol
    }
}

impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        self.eq_with_tol(other, 1e-6)
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Control {
    pub delta_s: f64,
    pub kappa: f64,
    pub sigma: f64,
}

impl Control {
    pub fn new(delta_s: f64, kappa: f64, sigma: f64) -> Self {
        Self {
            delta_s,
            kappa,
            sigma,
        }
    }
}
