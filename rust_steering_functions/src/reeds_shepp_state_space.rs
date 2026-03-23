use crate::utilities::{end_of_circular_arc, end_of_straight_line, pify, polar, sgn, PI};
use crate::{Control, State, StateSpace};
use std::f64;

const RS_EPS: f64 = 1e-6;
const RS_ZERO: f64 = 10.0 * std::f64::EPSILON;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReedsSheppPathSegmentType {
    Nop = 0,
    Left = 1,
    Straight = 2,
    Right = 3,
}

const REEDS_SHEPP_PATH_TYPE: [[ReedsSheppPathSegmentType; 5]; 18] = [
    [
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Nop,
        ReedsSheppPathSegmentType::Nop,
    ], // 0
    [
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Nop,
        ReedsSheppPathSegmentType::Nop,
    ], // 1
    [
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Nop,
    ], // 2
    [
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Nop,
    ], // 3
    [
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Nop,
    ], // 4
    [
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Nop,
    ], // 5
    [
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Nop,
    ], // 6
    [
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Nop,
    ], // 7
    [
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Nop,
    ], // 8
    [
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Nop,
    ], // 9
    [
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Nop,
    ], // 10
    [
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Nop,
    ], // 11
    [
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Nop,
        ReedsSheppPathSegmentType::Nop,
    ], // 12
    [
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Nop,
        ReedsSheppPathSegmentType::Nop,
    ], // 13
    [
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Nop,
        ReedsSheppPathSegmentType::Nop,
    ], // 14
    [
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Nop,
        ReedsSheppPathSegmentType::Nop,
    ], // 15
    [
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Right,
    ], // 16
    [
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
        ReedsSheppPathSegmentType::Straight,
        ReedsSheppPathSegmentType::Right,
        ReedsSheppPathSegmentType::Left,
    ], // 17
];

#[derive(Debug, Clone, Copy)]
pub struct ReedsSheppPath {
    pub types: [ReedsSheppPathSegmentType; 5],
    pub length: [f64; 5],
    pub total_length: f64,
}

impl ReedsSheppPath {
    pub fn new(
        types: [ReedsSheppPathSegmentType; 5],
        t: f64,
        u: f64,
        v: f64,
        w: f64,
        x: f64,
    ) -> Self {
        Self {
            types,
            length: [t, u, v, w, x],
            total_length: t.abs() + u.abs() + v.abs() + w.abs() + x.abs(),
        }
    }

    pub fn empty() -> Self {
        Self::new(REEDS_SHEPP_PATH_TYPE[0], f64::MAX, 0.0, 0.0, 0.0, 0.0)
    }

    pub fn length(&self) -> f64 {
        self.total_length
    }
}

pub struct ReedsSheppStateSpace {
    pub kappa: f64,
    pub kappa_inv: f64,
    pub discretization: f64,
}

impl ReedsSheppStateSpace {
    pub fn new(kappa: f64, discretization: f64) -> Self {
        assert!(kappa > 0.0 && discretization > 0.0);
        Self {
            kappa,
            kappa_inv: 1.0 / kappa,
            discretization,
        }
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
        let (u, t) = polar(x - phi.sin(), y - 1.0 + phi.cos());
        if t >= -RS_ZERO {
            let v = pify(phi - t);
            if v >= -RS_ZERO {
                return Some((t, u, v));
            }
        }
        None
    }

    fn lp_sp_rp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
        let (mut u1, t1) = polar(x + phi.sin(), y - 1.0 - phi.cos());
        u1 *= u1;
        if u1 >= 4.0 {
            let u = (u1 - 4.0).sqrt();
            let theta = (2.0_f64).atan2(u);
            let t = pify(t1 + theta);
            let v = pify(t - phi);
            if t >= -RS_ZERO && v >= -RS_ZERO {
                return Some((t, u, v));
            }
        }
        None
    }

    fn csc(x: f64, y: f64, phi: f64, paths: &mut Vec<ReedsSheppPath>) {
        if let Some((t, u, v)) = Self::lp_sp_lp(x, y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[14],
                t,
                u,
                v,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_sp_lp(-x, y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[14],
                -t,
                -u,
                -v,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_sp_lp(x, -y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[15],
                t,
                u,
                v,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_sp_lp(-x, -y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[15],
                -t,
                -u,
                -v,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_sp_rp(x, y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[12],
                t,
                u,
                v,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_sp_rp(-x, y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[12],
                -t,
                -u,
                -v,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_sp_rp(x, -y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[13],
                t,
                u,
                v,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_sp_rp(-x, -y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[13],
                -t,
                -u,
                -v,
                0.0,
                0.0,
            ));
        }
    }

    fn lp_rm_l(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
        let xi = x - phi.sin();
        let eta = y - 1.0 + phi.cos();
        let (u1, theta) = polar(xi, eta);
        if u1 <= 4.0 {
            let u = -2.0 * (0.25 * u1).asin();
            let t = pify(theta + 0.5 * u + PI);
            let v = pify(phi - t + u);
            if t >= -RS_ZERO && u <= RS_ZERO {
                return Some((t, u, v));
            }
        }
        None
    }

    fn ccc(x: f64, y: f64, phi: f64, paths: &mut Vec<ReedsSheppPath>) {
        if let Some((t, u, v)) = Self::lp_rm_l(x, y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[0],
                t,
                u,
                v,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_l(-x, y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[0],
                -t,
                -u,
                -v,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_l(x, -y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[1],
                t,
                u,
                v,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_l(-x, -y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[1],
                -t,
                -u,
                -v,
                0.0,
                0.0,
            ));
        }

        let xb = x * phi.cos() + y * phi.sin();
        let yb = x * phi.sin() - y * phi.cos();

        if let Some((t, u, v)) = Self::lp_rm_l(xb, yb, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[0],
                v,
                u,
                t,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_l(-xb, yb, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[0],
                -v,
                -u,
                -t,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_l(xb, -yb, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[1],
                v,
                u,
                t,
                0.0,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_l(-xb, -yb, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[1],
                -v,
                -u,
                -t,
                0.0,
                0.0,
            ));
        }
    }

    fn lp_rup_lum_rm(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
        let xi = x + phi.sin();
        let eta = y - 1.0 - phi.cos();
        let rho = 0.25 * (2.0 + (xi * xi + eta * eta).sqrt());
        if rho <= 1.0 {
            let u = rho.acos();
            let (t, v) = Self::tau_omega(u, -u, xi, eta, phi);
            if t >= -RS_ZERO && v <= RS_ZERO {
                return Some((t, u, v));
            }
        }
        None
    }

    fn lp_rum_lum_rp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
        let xi = x + phi.sin();
        let eta = y - 1.0 - phi.cos();
        let rho = (20.0 - xi * xi - eta * eta) / 16.0;
        if (0.0..=1.0).contains(&rho) {
            let u = -rho.acos();
            if u >= -0.5 * PI {
                let (t, v) = Self::tau_omega(u, u, xi, eta, phi);
                if t >= -RS_ZERO && v >= -RS_ZERO {
                    return Some((t, u, v));
                }
            }
        }
        None
    }

    fn cccc(x: f64, y: f64, phi: f64, paths: &mut Vec<ReedsSheppPath>) {
        if let Some((t, u, v)) = Self::lp_rup_lum_rm(x, y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[2],
                t,
                u,
                -u,
                v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rup_lum_rm(-x, y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[2],
                -t,
                -u,
                u,
                -v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rup_lum_rm(x, -y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[3],
                t,
                u,
                -u,
                v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rup_lum_rm(-x, -y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[3],
                -t,
                -u,
                u,
                -v,
                0.0,
            ));
        }

        if let Some((t, u, v)) = Self::lp_rum_lum_rp(x, y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[2],
                t,
                u,
                u,
                v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rum_lum_rp(-x, y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[2],
                -t,
                -u,
                -u,
                -v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rum_lum_rp(x, -y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[3],
                t,
                u,
                u,
                v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rum_lum_rp(-x, -y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[3],
                -t,
                -u,
                -u,
                -v,
                0.0,
            ));
        }
    }

    fn lp_rm_sm_lm(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
        let xi = x - phi.sin();
        let eta = y - 1.0 + phi.cos();
        let (rho, theta) = polar(xi, eta);
        if rho >= 2.0 {
            let r = (rho * rho - 4.0).sqrt();
            let u = 2.0 - r;
            let t = pify(theta + r.atan2(-2.0));
            let v = pify(phi - 0.5 * PI - t);
            if t >= -RS_ZERO && u <= RS_ZERO && v <= RS_ZERO {
                return Some((t, u, v));
            }
        }
        None
    }

    fn lp_rm_sm_rm(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
        let xi = x + phi.sin();
        let eta = y - 1.0 - phi.cos();
        let (rho, theta) = polar(-eta, xi);
        if rho >= 2.0 {
            let t = theta;
            let u = 2.0 - rho;
            let v = pify(t + 0.5 * PI - phi);
            if t >= -RS_ZERO && u <= RS_ZERO && v <= RS_ZERO {
                return Some((t, u, v));
            }
        }
        None
    }

    fn ccsc(x: f64, y: f64, phi: f64, paths: &mut Vec<ReedsSheppPath>) {
        if let Some((t, u, v)) = Self::lp_rm_sm_lm(x, y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[4],
                t,
                -0.5 * PI,
                u,
                v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_lm(-x, y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[4],
                -t,
                0.5 * PI,
                -u,
                -v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_lm(x, -y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[5],
                t,
                -0.5 * PI,
                u,
                v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_lm(-x, -y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[5],
                -t,
                0.5 * PI,
                -u,
                -v,
                0.0,
            ));
        }

        if let Some((t, u, v)) = Self::lp_rm_sm_rm(x, y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[8],
                t,
                -0.5 * PI,
                u,
                v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_rm(-x, y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[8],
                -t,
                0.5 * PI,
                -u,
                -v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_rm(x, -y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[9],
                t,
                -0.5 * PI,
                u,
                v,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_rm(-x, -y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[9],
                -t,
                0.5 * PI,
                -u,
                -v,
                0.0,
            ));
        }

        let xb = x * phi.cos() + y * phi.sin();
        let yb = x * phi.sin() - y * phi.cos();

        if let Some((t, u, v)) = Self::lp_rm_sm_lm(xb, yb, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[6],
                v,
                u,
                -0.5 * PI,
                t,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_lm(-xb, yb, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[6],
                -v,
                -u,
                0.5 * PI,
                -t,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_lm(xb, -yb, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[7],
                v,
                u,
                -0.5 * PI,
                t,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_lm(-xb, -yb, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[7],
                -v,
                -u,
                0.5 * PI,
                -t,
                0.0,
            ));
        }

        if let Some((t, u, v)) = Self::lp_rm_sm_rm(xb, yb, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[10],
                v,
                u,
                -0.5 * PI,
                t,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_rm(-xb, yb, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[10],
                -v,
                -u,
                0.5 * PI,
                -t,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_rm(xb, -yb, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[11],
                v,
                u,
                -0.5 * PI,
                t,
                0.0,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_sm_rm(-xb, -yb, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[11],
                -v,
                -u,
                0.5 * PI,
                -t,
                0.0,
            ));
        }
    }

    fn lp_rm_s_lm_rp(x: f64, y: f64, phi: f64) -> Option<(f64, f64, f64)> {
        let xi = x + phi.sin();
        let eta = y - 1.0 - phi.cos();
        let (rho, _theta) = polar(xi, eta);
        if rho >= 2.0 {
            let u = 4.0 - (rho * rho - 4.0).sqrt();
            if u <= RS_ZERO {
                let t = pify(((4.0 - u) * xi - 2.0 * eta).atan2(-2.0 * xi + (u - 4.0) * eta));
                let v = pify(t - phi);
                if t >= -RS_ZERO && v >= -RS_ZERO {
                    return Some((t, u, v));
                }
            }
        }
        None
    }

    fn ccscc(x: f64, y: f64, phi: f64, paths: &mut Vec<ReedsSheppPath>) {
        if let Some((t, u, v)) = Self::lp_rm_s_lm_rp(x, y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[16],
                t,
                -0.5 * PI,
                u,
                -0.5 * PI,
                v,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_s_lm_rp(-x, y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[16],
                -t,
                0.5 * PI,
                -u,
                0.5 * PI,
                -v,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_s_lm_rp(x, -y, -phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[17],
                t,
                -0.5 * PI,
                u,
                -0.5 * PI,
                v,
            ));
        }
        if let Some((t, u, v)) = Self::lp_rm_s_lm_rp(-x, -y, phi) {
            paths.push(ReedsSheppPath::new(
                REEDS_SHEPP_PATH_TYPE[17],
                -t,
                0.5 * PI,
                -u,
                0.5 * PI,
                -v,
            ));
        }
    }

    fn get_all_rs_paths_core(x: f64, y: f64, phi: f64) -> Vec<ReedsSheppPath> {
        let mut paths = Vec::with_capacity(48);
        Self::csc(x, y, phi, &mut paths);
        Self::ccc(x, y, phi, &mut paths);
        Self::cccc(x, y, phi, &mut paths);
        Self::ccsc(x, y, phi, &mut paths);
        Self::ccscc(x, y, phi, &mut paths);
        paths
    }

    pub fn reeds_shepp(&self, state1: &State, state2: &State) -> ReedsSheppPath {
        let dx = state2.x - state1.x;
        let dy = state2.y - state1.y;
        let dth = state2.theta - state1.theta;
        let c = state1.theta.cos();
        let s = state1.theta.sin();
        let x = c * dx + s * dy;
        let y = -s * dx + c * dy;

        let paths = Self::get_all_rs_paths_core(x * self.kappa, y * self.kappa, dth);

        let mut best_path = ReedsSheppPath::empty();
        let mut min_len = f64::MAX;

        for path in paths {
            let len = path.length();
            if len < min_len {
                min_len = len;
                best_path = path;
            }
        }
        best_path
    }

    pub fn get_all_rs_paths(&self, state1: &State, state2: &State) -> Vec<ReedsSheppPath> {
        let dx = state2.x - state1.x;
        let dy = state2.y - state1.y;
        let dth = state2.theta - state1.theta;
        let c = state1.theta.cos();
        let s = state1.theta.sin();
        let x = c * dx + s * dy;
        let y = -s * dx + c * dy;

        Self::get_all_rs_paths_core(x * self.kappa, y * self.kappa, dth)
    }

    fn extract_controls_from_path(&self, path: ReedsSheppPath) -> Vec<Control> {
        let mut controls = Vec::with_capacity(5);

        for i in 0..5 {
            let mut control = Control::default();
            match path.types[i] {
                ReedsSheppPathSegmentType::Nop => {
                    return controls;
                }
                ReedsSheppPathSegmentType::Left => {
                    control.delta_s = self.kappa_inv * path.length[i];
                    control.kappa = self.kappa;
                    control.sigma = 0.0;
                }
                ReedsSheppPathSegmentType::Right => {
                    control.delta_s = self.kappa_inv * path.length[i];
                    control.kappa = -self.kappa;
                    control.sigma = 0.0;
                }
                ReedsSheppPathSegmentType::Straight => {
                    control.delta_s = self.kappa_inv * path.length[i];
                    control.kappa = 0.0;
                    control.sigma = 0.0;
                }
            }
            if control.delta_s.abs() > RS_EPS {
                controls.push(control);
            }
        }
        controls
    }

    pub fn get_distance(&self, state1: &State, state2: &State) -> f64 {
        self.kappa_inv * self.reeds_shepp(state1, state2).length()
    }

    #[inline]
    pub fn integrate_ode(&self, state: &State, control: &Control, integration_step: f64) -> State {
        let mut state_next = *state;
        let kappa = control.kappa;
        let d = sgn(control.delta_s);

        if control.kappa.abs() > crate::EPSILON {
            let (x_f, y_f, theta_f) =
                end_of_circular_arc(state.x, state.y, state.theta, kappa, d, integration_step);
            state_next.x = x_f;
            state_next.y = y_f;
            state_next.theta = theta_f;
            state_next.kappa = kappa;
            state_next.d = d;
        } else {
            let (x_f, y_f) =
                end_of_straight_line(state.x, state.y, state.theta, d, integration_step);
            state_next.x = x_f;
            state_next.y = y_f;
            state_next.theta = state.theta;
            state_next.kappa = kappa;
            state_next.d = d;
        }
        state_next.s = state.s + integration_step;
        state_next
    }
}

impl StateSpace for ReedsSheppStateSpace {
    fn get_controls(&self, state1: &State, state2: &State) -> Vec<Control> {
        let path = self.reeds_shepp(state1, state2);
        self.extract_controls_from_path(path)
    }

    fn get_all_controls(&self, state1: &State, state2: &State) -> Vec<Vec<Control>> {
        let mut all_controls = Vec::new();
        let all_rs_paths = self.get_all_rs_paths(state1, state2);

        if (state1.x - state2.x).abs() < 1e-6
            && (state1.y - state2.y).abs() < 1e-6
            && (state1.theta - state2.theta).abs() < 1e-6
        {
            return vec![vec![Control::new(0.0, self.kappa, 0.0)]];
        }

        for rs_path in all_rs_paths {
            let control = self.extract_controls_from_path(rs_path);
            if !control.is_empty() {
                all_controls.push(control);
            }
        }
        all_controls
    }

    fn get_path(&self, state1: &State, state2: &State) -> Vec<State> {
        let controls = self.get_controls(state1, state2);
        self.integrate(state1, &controls)
    }

    fn integrate(&self, state: &State, controls: &[Control]) -> Vec<State> {
        let mut path = Vec::new();
        if controls.is_empty() {
            return path;
        }

        let mut n_states = 0;
        for control in controls {
            n_states += (control.delta_s.abs() / self.discretization).ceil() as usize;
        }
        path.reserve(n_states + 10);

        let mut state_curr = *state;
        state_curr.kappa = controls[0].kappa;
        state_curr.d = sgn(controls[0].delta_s);
        path.push(state_curr);

        for control in controls {
            let abs_delta_s = control.delta_s.abs();
            let n = f64::max(2.0, (abs_delta_s / self.discretization).ceil()) as usize;
            let discretization = abs_delta_s / (n as f64);

            for _ in 0..n {
                let state_next = self.integrate_ode(&state_curr, control, discretization);
                path.push(state_next);
                state_curr = state_next;
            }
        }
        path
    }

    fn interpolate(&self, state: &State, controls: &[Control], mut t: f64) -> State {
        let mut state_curr = *state;

        let mut s_path = 0.0;
        for control in controls {
            s_path += control.delta_s.abs();
        }
        t = t.clamp(0.0, 1.0);
        let s_inter = t * s_path;

        let mut s = 0.0;
        for control in controls {
            let abs_delta_s = control.delta_s.abs();

            if s_inter - s > abs_delta_s {
                state_curr = self.integrate_ode(&state_curr, control, abs_delta_s);
                s += abs_delta_s;
            } else {
                return self.integrate_ode(&state_curr, control, s_inter - s);
            }
        }
        state_curr
    }
}
