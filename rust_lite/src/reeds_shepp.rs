use std::f64::consts::PI;
use crate::state::{State, Control};
use crate::base_state_space::StateSpace;
use crate::utilities::pify;

const RS_NOP: u8 = 0;
const RS_LEFT: u8 = 1;
const RS_STRAIGHT: u8 = 2;
const RS_RIGHT: u8 = 3;

const RS_PATH_TYPE: [[u8; 5]; 18] = [
    [RS_LEFT, RS_RIGHT, RS_LEFT, RS_NOP, RS_NOP],
    [RS_RIGHT, RS_LEFT, RS_RIGHT, RS_NOP, RS_NOP],
    [RS_LEFT, RS_RIGHT, RS_LEFT, RS_RIGHT, RS_NOP],
    [RS_RIGHT, RS_LEFT, RS_RIGHT, RS_LEFT, RS_NOP],
    [RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_NOP],
    [RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_NOP],
    [RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_LEFT, RS_NOP],
    [RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_RIGHT, RS_NOP],
    [RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP],
    [RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_NOP],
    [RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_LEFT, RS_NOP],
    [RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_RIGHT, RS_NOP],
    [RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_NOP, RS_NOP],
    [RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_NOP, RS_NOP],
    [RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_NOP, RS_NOP],
    [RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP, RS_NOP],
    [RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_RIGHT],
    [RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_LEFT],
];

const RS_EPS: f64 = 1e-6;
const RS_ZERO: f64 = 10.0 * f64::EPSILON;

#[derive(Clone, Debug)]
struct RsPath {
    type_: [u8; 5],
    length_: [f64; 5],
    total_: f64,
}

impl RsPath {
    fn new(type_: [u8; 5], t: f64, u: f64, v: f64, w: f64, x: f64) -> Self {
        let total_ = t.abs() + u.abs() + v.abs() + w.abs() + x.abs();
        Self { type_, length_: [t, u, v, w, x], total_ }
    }
    fn invalid() -> Self {
        Self {
            type_: RS_PATH_TYPE[0],
            length_: [f64::INFINITY, 0.0, 0.0, 0.0, 0.0],
            total_: f64::INFINITY,
        }
    }
    fn length(&self) -> f64 { self.total_ }
}

fn tau_omega(u: f64, v: f64, xi: f64, eta: f64, phi: f64) -> (f64, f64) {
    let delta = pify(u - v);
    let a = u.sin() - delta.sin();
    let b = u.cos() - delta.cos() - 1.0;
    let t1 = (eta * a - xi * b).atan2(xi * a + eta * b);
    let t2 = 2.0 * (delta.cos() - v.cos() - u.cos()) + 3.0;
    let tau = if t2 < 0.0 { pify(t1 + PI) } else { pify(t1) };
    let omega = pify(tau - u + v - phi);
    (tau, omega)
}

fn lp_sp_lp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    use crate::utilities::polar;
    let (u, t) = polar(x - phi.sin(), y - 1.0 + phi.cos());
    if t >= -RS_ZERO {
        let v = pify(phi - t);
        if v >= -RS_ZERO { return Some((t, u, v)); }
    }
    None
}

fn lp_sp_rp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    use crate::utilities::polar;
    let (u1, t1) = polar(x + phi.sin(), y - 1.0 - phi.cos());
    let u1sq = u1 * u1;
    if u1sq >= 4.0 {
        let u = (u1sq - 4.0).sqrt();
        let theta = (2.0_f64).atan2(u);
        let t = pify(t1 + theta);
        let v = pify(t - phi);
        if t >= -RS_ZERO && v >= -RS_ZERO { return Some((t, u, v)); }
    }
    None
}

fn lp_rm_l(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    use crate::utilities::polar;
    let xi = x - phi.sin();
    let eta = y - 1.0 + phi.cos();
    let (u1, theta) = polar(xi, eta);
    if u1 <= 4.0 {
        let u = -2.0 * (0.25 * u1).asin();
        let t = pify(theta + 0.5 * u + PI);
        let v = pify(phi - t + u);
        if t >= -RS_ZERO && u <= RS_ZERO { return Some((t, u, v)); }
    }
    None
}

fn lp_rup_lum_rm(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    let xi = x + phi.sin();
    let eta = y - 1.0 - phi.cos();
    let rho = 0.25 * (2.0 + (xi * xi + eta * eta).sqrt());
    if rho <= 1.0 {
        let u = rho.acos();
        let (t, v) = tau_omega(u, -u, xi, eta, phi);
        if t >= -RS_ZERO && v <= RS_ZERO { return Some((t, u, v)); }
    }
    None
}

fn lp_rum_lum_rp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    let xi = x + phi.sin();
    let eta = y - 1.0 - phi.cos();
    let rho = (20.0 - xi * xi - eta * eta) / 16.0;
    if rho >= 0.0 && rho <= 1.0 {
        let u = -rho.acos();
        if u >= -0.5 * PI {
            let (t, v) = tau_omega(u, u, xi, eta, phi);
            if t >= -RS_ZERO && v >= -RS_ZERO { return Some((t, u, v)); }
        }
    }
    None
}

fn lp_rm_sm_lm(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    use crate::utilities::polar;
    let xi = x - phi.sin();
    let eta = y - 1.0 + phi.cos();
    let (rho, theta) = polar(xi, eta);
    if rho >= 2.0 {
        let r = (rho * rho - 4.0).sqrt();
        let u = 2.0 - r;
        let t = pify(theta + r.atan2(-2.0));
        let v = pify(phi - 0.5 * PI - t);
        if t >= -RS_ZERO && u <= RS_ZERO && v <= RS_ZERO { return Some((t, u, v)); }
    }
    None
}

fn lp_rm_sm_rm(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    use crate::utilities::polar;
    let xi = x + phi.sin();
    let eta = y - 1.0 - phi.cos();
    let (rho, theta) = polar(-eta, xi);
    if rho >= 2.0 {
        let t = theta;
        let u = 2.0 - rho;
        let v = pify(t + 0.5 * PI - phi);
        if t >= -RS_ZERO && u <= RS_ZERO && v <= RS_ZERO { return Some((t, u, v)); }
    }
    None
}

fn lp_rm_s_lm_rp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
    use crate::utilities::polar;
    let xi = x + phi.sin();
    let eta = y - 1.0 - phi.cos();
    let (rho, _theta) = polar(xi, eta);
    if rho >= 2.0 {
        let u = 4.0 - (rho * rho - 4.0).sqrt();
        if u <= RS_ZERO {
            let t = pify(((4.0 - u) * xi - 2.0 * eta).atan2(-2.0 * xi + (u - 4.0) * eta));
            let v = pify(t - phi);
            if t >= -RS_ZERO && v >= -RS_ZERO { return Some((t, u, v)); }
        }
    }
    None
}

fn collect_all_paths(x: f64, y: f64, phi: f64) -> Vec<RsPath> {
    let mut paths: Vec<RsPath> = Vec::new();
    let xb = x * phi.cos() + y * phi.sin();
    let yb = x * phi.sin() - y * phi.cos();

    // CSC
    for &(nx, ny, nphi, pidx, sign) in &[
        (x, y, phi, 14usize, 1.0_f64),
        (-x, y, -phi, 14, -1.0),
        (x, -y, -phi, 15, 1.0),
        (-x, -y, phi, 15, -1.0),
    ] {
        if let Some((t, u, v)) = lp_sp_lp(nx, ny, nphi) {
            paths.push(RsPath::new(RS_PATH_TYPE[pidx], sign*t, sign*u, sign*v, 0.0, 0.0));
        }
    }
    for &(nx, ny, nphi, pidx, sign) in &[
        (x, y, phi, 12usize, 1.0_f64),
        (-x, y, -phi, 12, -1.0),
        (x, -y, -phi, 13, 1.0),
        (-x, -y, phi, 13, -1.0),
    ] {
        if let Some((t, u, v)) = lp_sp_rp(nx, ny, nphi) {
            paths.push(RsPath::new(RS_PATH_TYPE[pidx], sign*t, sign*u, sign*v, 0.0, 0.0));
        }
    }
    // CCC
    for &(nx, ny, nphi, pidx, sign, rev) in &[
        (x, y, phi, 0usize, 1.0_f64, false),
        (-x, y, -phi, 0, -1.0, false),
        (x, -y, -phi, 1, 1.0, false),
        (-x, -y, phi, 1, -1.0, false),
        (xb, yb, phi, 0, 1.0, true),
        (-xb, yb, -phi, 0, -1.0, true),
        (xb, -yb, -phi, 1, 1.0, true),
        (-xb, -yb, phi, 1, -1.0, true),
    ] {
        if let Some((t, u, v)) = lp_rm_l(nx, ny, nphi) {
            let (a, b, c) = if rev { (v, u, t) } else { (t, u, v) };
            paths.push(RsPath::new(RS_PATH_TYPE[pidx], sign*a, sign*b, sign*c, 0.0, 0.0));
        }
    }
    // CCCC
    for &(nx, ny, nphi, pidx, sign, neg_mid) in &[
        (x, y, phi, 2usize, 1.0_f64, true),
        (-x, y, -phi, 2, -1.0, true),
        (x, -y, -phi, 3, 1.0, true),
        (-x, -y, phi, 3, -1.0, true),
    ] {
        if let Some((t, u, v)) = lp_rup_lum_rm(nx, ny, nphi) {
            let u2 = if neg_mid { -u } else { u };
            paths.push(RsPath::new(RS_PATH_TYPE[pidx], sign*t, sign*u, sign*u2, sign*v, 0.0));
        }
    }
    for &(nx, ny, nphi, pidx, sign) in &[
        (x, y, phi, 2usize, 1.0_f64),
        (-x, y, -phi, 2, -1.0),
        (x, -y, -phi, 3, 1.0),
        (-x, -y, phi, 3, -1.0),
    ] {
        if let Some((t, u, v)) = lp_rum_lum_rp(nx, ny, nphi) {
            paths.push(RsPath::new(RS_PATH_TYPE[pidx], sign*t, sign*u, sign*u, sign*v, 0.0));
        }
    }
    // CCSC forward
    for &(nx, ny, nphi, pidx, sign) in &[
        (x, y, phi, 4usize, 1.0_f64),
        (-x, y, -phi, 4, -1.0),
        (x, -y, -phi, 5, 1.0),
        (-x, -y, phi, 5, -1.0),
    ] {
        if let Some((t, u, v)) = lp_rm_sm_lm(nx, ny, nphi) {
            paths.push(RsPath::new(RS_PATH_TYPE[pidx], sign*t, sign*(-0.5*PI), sign*u, sign*v, 0.0));
        }
    }
    for &(nx, ny, nphi, pidx, sign) in &[
        (x, y, phi, 8usize, 1.0_f64),
        (-x, y, -phi, 8, -1.0),
        (x, -y, -phi, 9, 1.0),
        (-x, -y, phi, 9, -1.0),
    ] {
        if let Some((t, u, v)) = lp_rm_sm_rm(nx, ny, nphi) {
            paths.push(RsPath::new(RS_PATH_TYPE[pidx], sign*t, sign*(-0.5*PI), sign*u, sign*v, 0.0));
        }
    }
    // CCSC backwards
    for &(nx, ny, nphi, pidx, sign) in &[
        (xb, yb, phi, 6usize, 1.0_f64),
        (-xb, yb, -phi, 6, -1.0),
        (xb, -yb, -phi, 7, 1.0),
        (-xb, -yb, phi, 7, -1.0),
    ] {
        if let Some((t, u, v)) = lp_rm_sm_lm(nx, ny, nphi) {
            paths.push(RsPath::new(RS_PATH_TYPE[pidx], sign*v, sign*u, sign*(-0.5*PI), sign*t, 0.0));
        }
    }
    for &(nx, ny, nphi, pidx, sign) in &[
        (xb, yb, phi, 10usize, 1.0_f64),
        (-xb, yb, -phi, 10, -1.0),
        (xb, -yb, -phi, 11, 1.0),
        (-xb, -yb, phi, 11, -1.0),
    ] {
        if let Some((t, u, v)) = lp_rm_sm_rm(nx, ny, nphi) {
            paths.push(RsPath::new(RS_PATH_TYPE[pidx], sign*v, sign*u, sign*(-0.5*PI), sign*t, 0.0));
        }
    }
    // CCSCC
    for &(nx, ny, nphi, pidx, sign) in &[
        (x, y, phi, 16usize, 1.0_f64),
        (-x, y, -phi, 16, -1.0),
        (x, -y, -phi, 17, 1.0),
        (-x, -y, phi, 17, -1.0),
    ] {
        if let Some((t, u, v)) = lp_rm_s_lm_rp(nx, ny, nphi) {
            paths.push(RsPath::new(RS_PATH_TYPE[pidx], sign*t, sign*(-0.5*PI), sign*u, sign*(-0.5*PI), sign*v));
        }
    }

    paths
}

fn shortest_rs_path(x: f64, y: f64, phi: f64) -> RsPath {
    let paths = collect_all_paths(x, y, phi);
    let mut best = RsPath::invalid();
    for p in paths {
        if p.length() < best.length() { best = p; }
    }
    best
}

fn normalize(s1: &State, s2: &State, kappa: f64) -> (f64, f64, f64) {
    let dx = s2.x - s1.x;
    let dy = s2.y - s1.y;
    let c = s1.theta.cos();
    let s = s1.theta.sin();
    let x = (c * dx + s * dy) * kappa;
    let y = (-s * dx + c * dy) * kappa;
    let phi = s2.theta - s1.theta;
    (x, y, phi)
}

fn controls_from_rs(path: &RsPath, kappa_inv: f64) -> Vec<Control> {
    let mut controls = Vec::new();
    for i in 0..5 {
        let seg = path.type_[i];
        if seg == RS_NOP { break; }
        let delta_s = kappa_inv * path.length_[i];
        if delta_s.abs() <= RS_EPS { continue; }
        let kappa = match seg {
            RS_LEFT  =>  1.0 / kappa_inv,
            RS_RIGHT => -1.0 / kappa_inv,
            _        =>  0.0,
        };
        controls.push(Control { delta_s, kappa, sigma: 0.0 });
    }
    controls
}

/// Reeds-Shepp state space.
pub struct ReedsSheppStateSpace {
    kappa_: f64,
    discretization_: f64,
}

impl ReedsSheppStateSpace {
    pub fn new(kappa: f64, discretization: f64) -> Self {
        assert!(kappa > 0.0 && discretization > 0.0);
        Self { kappa_: kappa, discretization_: discretization }
    }

    pub fn get_distance(&self, s1: &State, s2: &State) -> f64 {
        let (x, y, phi) = normalize(s1, s2, self.kappa_);
        shortest_rs_path(x, y, phi).length() / self.kappa_
    }
}

impl StateSpace for ReedsSheppStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let (x, y, phi) = normalize(s1, s2, self.kappa_);
        let path = shortest_rs_path(x, y, phi);
        controls_from_rs(&path, 1.0 / self.kappa_)
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        let (x, y, phi) = normalize(s1, s2, self.kappa_);
        collect_all_paths(x, y, phi)
            .iter()
            .map(|p| controls_from_rs(p, 1.0 / self.kappa_))
            .filter(|c| !c.is_empty())
            .collect()
    }

    fn discretization(&self) -> f64 { self.discretization_ }
}
