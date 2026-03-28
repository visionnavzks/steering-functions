use crate::state::{State, Control};
use crate::base_state_space::StateSpace;
use crate::utilities::twopify;

const DUBINS_LEFT: u8 = 0;
const DUBINS_STRAIGHT: u8 = 1;
const DUBINS_RIGHT: u8 = 2;

const DUBINS_PATH_TYPE: [[u8; 3]; 6] = [
    [DUBINS_LEFT, DUBINS_STRAIGHT, DUBINS_LEFT],
    [DUBINS_RIGHT, DUBINS_STRAIGHT, DUBINS_RIGHT],
    [DUBINS_RIGHT, DUBINS_STRAIGHT, DUBINS_LEFT],
    [DUBINS_LEFT, DUBINS_STRAIGHT, DUBINS_RIGHT],
    [DUBINS_RIGHT, DUBINS_LEFT, DUBINS_RIGHT],
    [DUBINS_LEFT, DUBINS_RIGHT, DUBINS_LEFT],
];

const DUBINS_EPS: f64 = 1e-6;
const DUBINS_ZERO: f64 = -1e-9;

#[derive(Clone, Debug)]
struct DubinsPath {
    type_: [u8; 3],
    length_: [f64; 3],
}

impl DubinsPath {
    fn new(type_: [u8; 3], length_: [f64; 3]) -> Self {
        Self { type_, length_ }
    }
    fn invalid() -> Self {
        Self {
            type_: DUBINS_PATH_TYPE[0],
            length_: [0.0, f64::INFINITY, 0.0],
        }
    }
    fn length(&self) -> f64 {
        self.length_[0] + self.length_[1] + self.length_[2]
    }
}

fn dubins_lsl(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let (ca, sa) = (alpha.cos(), alpha.sin());
    let (cb, sb) = (beta.cos(), beta.sin());
    let tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sa - sb));
    if tmp >= DUBINS_ZERO {
        let theta = (cb - ca).atan2(d + sa - sb);
        let t = twopify(-alpha + theta);
        let p = tmp.max(0.0).sqrt();
        let q = twopify(beta - theta);
        DubinsPath::new(DUBINS_PATH_TYPE[0], [t, p, q])
    } else {
        DubinsPath::invalid()
    }
}

fn dubins_rsr(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let (ca, sa) = (alpha.cos(), alpha.sin());
    let (cb, sb) = (beta.cos(), beta.sin());
    let tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sb - sa));
    if tmp >= DUBINS_ZERO {
        let theta = (ca - cb).atan2(d - sa + sb);
        let t = twopify(alpha - theta);
        let p = tmp.max(0.0).sqrt();
        let q = twopify(-beta + theta);
        DubinsPath::new(DUBINS_PATH_TYPE[1], [t, p, q])
    } else {
        DubinsPath::invalid()
    }
}

fn dubins_rsl(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let (ca, sa) = (alpha.cos(), alpha.sin());
    let (cb, sb) = (beta.cos(), beta.sin());
    let tmp = d * d - 2.0 + 2.0 * (ca * cb + sa * sb - d * (sa + sb));
    if tmp >= DUBINS_ZERO {
        let p = tmp.max(0.0).sqrt();
        let theta = (ca + cb).atan2(d - sa - sb) - (2.0_f64).atan2(p);
        let t = twopify(alpha - theta);
        let q = twopify(beta - theta);
        DubinsPath::new(DUBINS_PATH_TYPE[2], [t, p, q])
    } else {
        DubinsPath::invalid()
    }
}

fn dubins_lsr(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let (ca, sa) = (alpha.cos(), alpha.sin());
    let (cb, sb) = (beta.cos(), beta.sin());
    let tmp = -2.0 + d * d + 2.0 * (ca * cb + sa * sb + d * (sa + sb));
    if tmp >= DUBINS_ZERO {
        let p = tmp.max(0.0).sqrt();
        let theta = (-ca - cb).atan2(d + sa + sb) - (-2.0_f64).atan2(p);
        let t = twopify(-alpha + theta);
        let q = twopify(-beta + theta);
        DubinsPath::new(DUBINS_PATH_TYPE[3], [t, p, q])
    } else {
        DubinsPath::invalid()
    }
}

fn dubins_rlr(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let (ca, sa) = (alpha.cos(), alpha.sin());
    let (cb, sb) = (beta.cos(), beta.sin());
    let tmp = 0.125 * (6.0 - d * d + 2.0 * (ca * cb + sa * sb + d * (sa - sb)));
    if tmp.abs() <= 1.0 {
        let p = twopify(std::f64::consts::TAU - tmp.acos());
        let theta = (ca - cb).atan2(d - sa + sb);
        let t = twopify(alpha - theta + 0.5 * p);
        let q = twopify(alpha - beta - t + p);
        DubinsPath::new(DUBINS_PATH_TYPE[4], [t, p, q])
    } else {
        DubinsPath::invalid()
    }
}

fn dubins_lrl(d: f64, alpha: f64, beta: f64) -> DubinsPath {
    let (ca, sa) = (alpha.cos(), alpha.sin());
    let (cb, sb) = (beta.cos(), beta.sin());
    let tmp = 0.125 * (6.0 - d * d + 2.0 * (ca * cb + sa * sb - d * (sa - sb)));
    if tmp.abs() <= 1.0 {
        let p = twopify(std::f64::consts::TAU - tmp.acos());
        let theta = (-ca + cb).atan2(d + sa - sb);
        let t = twopify(-alpha + theta + 0.5 * p);
        let q = twopify(beta - alpha - t + p);
        DubinsPath::new(DUBINS_PATH_TYPE[5], [t, p, q])
    } else {
        DubinsPath::invalid()
    }
}

fn dubins_word(path_type: usize, d: f64, alpha: f64, beta: f64) -> DubinsPath {
    match path_type {
        0 => dubins_lsl(d, alpha, beta),
        1 => dubins_rsr(d, alpha, beta),
        2 => dubins_rsl(d, alpha, beta),
        3 => dubins_lsr(d, alpha, beta),
        4 => dubins_rlr(d, alpha, beta),
        5 => dubins_lrl(d, alpha, beta),
        _ => DubinsPath::invalid(),
    }
}

fn shortest_dubins_path(
    q0: &State, q1: &State, rho: f64, forward: bool,
) -> (DubinsPath, f64) {
    let dx = q1.x - q0.x;
    let dy = q1.y - q0.y;
    let d = dx.hypot(dy) / rho;
    let (theta0, theta1) = if forward {
        (q0.theta, q1.theta)
    } else {
        (q0.theta + std::f64::consts::PI, q1.theta + std::f64::consts::PI)
    };
    let alpha = twopify(theta0 - dy.atan2(dx));
    let beta  = twopify(theta1 - dy.atan2(dx));

    let mut best = DubinsPath::invalid();
    let mut best_len = f64::INFINITY;
    for i in 0..6 {
        let p = dubins_word(i, d, alpha, beta);
        let l = p.length();
        if l < best_len {
            best_len = l;
            best = p;
        }
    }
    (best, best_len * rho)
}

fn controls_from_dubins(path: &DubinsPath, rho: f64, forward: bool, dx: f64, dy: f64) -> Vec<Control> {
    let d = if forward { 1.0 } else { -1.0 };
    let mut controls = Vec::new();
    for i in 0..3 {
        let seg_type = path.type_[i];
        let length = path.length_[i] * rho;
        if length.abs() < DUBINS_EPS { continue; }
        let (kappa, sigma) = match seg_type {
            t if t == DUBINS_LEFT     => (1.0 / rho,  0.0),
            t if t == DUBINS_STRAIGHT => (0.0,         0.0),
            _                         => (-1.0 / rho, 0.0),
        };
        controls.push(Control { delta_s: d * length, kappa, sigma });
    }
    controls
}

/// Dubins state space.
pub struct DubinsStateSpace {
    kappa_: f64,
    discretization_: f64,
    forward_: bool,
}

impl DubinsStateSpace {
    pub fn new(kappa: f64, discretization: f64, forward: bool) -> Self {
        assert!(kappa > 0.0 && discretization > 0.0);
        Self { kappa_: kappa, discretization_: discretization, forward_: forward }
    }
}

impl StateSpace for DubinsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let rho = 1.0 / self.kappa_;
        let dx = s2.x - s1.x;
        let dy = s2.y - s1.y;
        let (path, _) = shortest_dubins_path(s1, s2, rho, self.forward_);
        controls_from_dubins(&path, rho, self.forward_, dx, dy)
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        let rho = 1.0 / self.kappa_;
        let dx = s2.x - s1.x;
        let dy = s2.y - s1.y;
        let mut all = Vec::new();
        let d = dx.hypot(dy) / rho;
        let alpha = twopify(s1.theta - dy.atan2(dx));
        let beta  = twopify(s2.theta - dy.atan2(dx));
        for i in 0..6 {
            let p = dubins_word(i, d, alpha, beta);
            if p.length().is_finite() {
                all.push(controls_from_dubins(&p, rho, self.forward_, dx, dy));
            }
        }
        all
    }

    fn discretization(&self) -> f64 { self.discretization_ }
}
