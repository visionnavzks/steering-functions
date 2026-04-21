use std::f64::consts::PI;
use crate::utilities::{get_epsilon, global_frame_change, HALF_PI, sgn, end_of_clothoid};
use crate::configuration::{Configuration, configuration_distance, configuration_equal, configuration_aligned};
use crate::hc_cc_circle::{HcCcCircle, HcCcCircleParam, center_distance, configuration_on_hc_cc_circle};
use crate::paths::{
    HcCcRsPath, HcCcRsPathType, NB_HC_CC_RS_PATHS,
    empty_controls, straight_controls, cc_turn_controls, hc_turn_controls, rs_turn_controls,
    state_equal, reverse_control, subtract_control,
};
use crate::hc_cc_state_space::HcCcStateSpaceParams;
use crate::state::{State, Control};
use crate::base_state_space::StateSpace;

const CC_REGULAR: bool = false;
const HC_REGULAR: bool = false;

// ---------------------------------------------------------------------------
// Per-candidate accumulator
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct PathSlot {
    length: f64,
    cstart: Option<HcCcCircle>,
    cend:   Option<HcCcCircle>,
    qi1: Option<Configuration>,
    qi2: Option<Configuration>,
    qi3: Option<Configuration>,
    qi4: Option<Configuration>,
    ci1: Option<HcCcCircle>,
    ci2: Option<HcCcCircle>,
}

impl PathSlot {
    fn infinite() -> Self {
        Self {
            length: f64::INFINITY,
            cstart: None, cend: None,
            qi1: None, qi2: None, qi3: None, qi4: None,
            ci1: None, ci2: None,
        }
    }
}

fn path_type_from_usize(i: usize) -> HcCcRsPathType {
    match i {
        0  => HcCcRsPathType::E,
        1  => HcCcRsPathType::S,
        2  => HcCcRsPathType::T,
        3  => HcCcRsPathType::TT,
        4  => HcCcRsPathType::TcT,
        5  => HcCcRsPathType::TcTcT,
        6  => HcCcRsPathType::TcTT,
        7  => HcCcRsPathType::TTcT,
        8  => HcCcRsPathType::TST,
        9  => HcCcRsPathType::TSTcT,
        10 => HcCcRsPathType::TcTST,
        11 => HcCcRsPathType::TcTSTcT,
        12 => HcCcRsPathType::TTcTT,
        13 => HcCcRsPathType::TcTTcT,
        14 => HcCcRsPathType::TTT,
        15 => HcCcRsPathType::TcST,
        16 => HcCcRsPathType::TScT,
        17 => HcCcRsPathType::TcScT,
        _  => panic!("invalid path type index"),
    }
}

// ---------------------------------------------------------------------------
// HC00 geometry helper
// ---------------------------------------------------------------------------

struct Hc00RsHelper<'a> {
    param:    &'a HcCcCircleParam,
    rs_param: &'a HcCcCircleParam,
    distance: f64,
    angle:    f64,
}

#[allow(non_snake_case)]
impl<'a> Hc00RsHelper<'a> {
    fn new(param: &'a HcCcCircleParam, rs_param: &'a HcCcCircleParam) -> Self {
        Self { param, rs_param, distance: 0.0, angle: 0.0 }
    }

    // ---- TT ---------------------------------------------------------------
    fn tt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left
            && c1.forward != c2.forward
            && (self.distance - 2.0 * c1.radius).abs() < get_epsilon()
    }

    fn tt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Configuration {
        let x = (c1.xc + c2.xc) / 2.0;
        let y = (c1.yc + c2.yc) / 2.0;
        let angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        let theta = if c1.left {
            if c1.forward { angle + HALF_PI - c1.mu } else { angle + HALF_PI + c1.mu }
        } else {
            if c1.forward { angle - HALF_PI + c1.mu } else { angle - HALF_PI - c1.mu }
        };
        Configuration::new(x, y, theta, 0.0)
    }

    fn tt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let q = self.tt_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.cc_turn_length(&q) + cend.cc_turn_length(&q);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q), ..PathSlot::infinite() }
    }

    // ---- TcT --------------------------------------------------------------
    fn tct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left
            && c1.forward == c2.forward
            && (self.distance - 2.0 * c1.kappa_inv.abs()).abs() < get_epsilon()
    }

    fn tct_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Configuration {
        let dist    = center_distance(c1, c2);
        let delta_x = 0.5 * dist;
        let delta_y = 0.0;
        let angle   = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        let theta;
        let (x, y) = if c1.left {
            theta = angle + HALF_PI;
            if c1.forward { global_frame_change(c1.xc, c1.yc, angle, delta_x,  delta_y) }
            else          { global_frame_change(c1.xc, c1.yc, angle, delta_x, -delta_y) }
        } else {
            theta = angle - HALF_PI;
            if c1.forward { global_frame_change(c1.xc, c1.yc, angle, delta_x, -delta_y) }
            else          { global_frame_change(c1.xc, c1.yc, angle, delta_x,  delta_y) }
        };
        Configuration::new(x, y, theta, c1.kappa)
    }

    fn tct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let q = self.tct_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let length = cstart.hc_turn_length(&q) + cend.hc_turn_length(&q);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q), ..PathSlot::infinite() }
    }

    // ---- TcTcT ------------------------------------------------------------
    fn tctct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left
            && c1.forward != c2.forward
            && self.distance <= 4.0 * c1.kappa_inv.abs()
    }

    fn tctct_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration)
    {
        let theta   = self.angle;
        let r       = 2.0 * c1.kappa_inv.abs();
        let delta_x = 0.5 * self.distance;
        let delta_y = (r * r - delta_x * delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.param);
        let q1 = self.tct_tangent_circles(c1,    &tgt1);
        let q2 = self.tct_tangent_circles(&tgt1,  c2);
        let q3 = self.tct_tangent_circles(c1,    &tgt2);
        let q4 = self.tct_tangent_circles(&tgt2,  c2);
        (q1, q2, q3, q4)
    }

    fn tctct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.tctct_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, !c1.forward, true, self.rs_param);
        let mid2   = HcCcCircle::from_configuration(qc, !c1.left, !c1.forward, true, self.rs_param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let l1 = cstart.hc_turn_length(&qa) + mid1.rs_turn_length(&qb) + cend.hc_turn_length(&qb);
        let l2 = cstart.hc_turn_length(&qc) + mid2.rs_turn_length(&qd) + cend.hc_turn_length(&qd);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qa), qi2: Some(qb), ci1: Some(mid1), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qc), qi2: Some(qd), ci1: Some(mid2), ..PathSlot::infinite() }
        }
    }

    // ---- TcTT -------------------------------------------------------------
    fn tctt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward == c2.forward
            && self.distance <= 2.0 * c1.radius + 2.0 * c1.kappa_inv.abs()
            && self.distance >= (2.0 * c1.radius - 2.0 * c1.kappa_inv.abs()).abs()
    }

    fn tctt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration)
    {
        let theta   = self.angle;
        let r1      = 2.0 * c1.kappa_inv.abs();
        let r2      = 2.0 * c1.radius;
        let delta_x = (r1*r1 + self.distance*self.distance - r2*r2) / (2.0 * self.distance);
        let delta_y = (r1*r1 - delta_x*delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.param);
        let q1 = self.tct_tangent_circles(c1,   &tgt1);
        let q2 = self.tt_tangent_circles(&tgt1,  c2);
        let q3 = self.tct_tangent_circles(c1,   &tgt2);
        let q4 = self.tt_tangent_circles(&tgt2,  c2);
        (q1, q2, q3, q4)
    }

    fn tctt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.tctt_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qb, !c1.left, c1.forward, true, self.param);
        let mid2   = HcCcCircle::from_configuration(qd, !c1.left, c1.forward, true, self.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let l1 = cstart.hc_turn_length(&qa) + mid1.hc_turn_length(&qa) + cend.cc_turn_length(&qb);
        let l2 = cstart.hc_turn_length(&qc) + mid2.hc_turn_length(&qc) + cend.cc_turn_length(&qd);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qa), qi2: Some(qb), ci1: Some(mid1), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qc), qi2: Some(qd), ci1: Some(mid2), ..PathSlot::infinite() }
        }
    }

    // ---- TTcT -------------------------------------------------------------
    fn ttct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward == c2.forward
            && self.distance <= 2.0 * c1.radius + 2.0 * c1.kappa_inv.abs()
            && self.distance >= (2.0 * c1.radius - 2.0 * c1.kappa_inv.abs()).abs()
    }

    fn ttct_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration)
    {
        let theta   = self.angle;
        let r1      = 2.0 * c1.radius;
        let r2      = 2.0 * c1.kappa_inv.abs();
        let delta_x = (r1*r1 + self.distance*self.distance - r2*r2) / (2.0 * self.distance);
        let delta_y = (r1*r1 - delta_x*delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.param);
        let q1 = self.tt_tangent_circles(c1,   &tgt1);
        let q2 = self.tct_tangent_circles(&tgt1, c2);
        let q3 = self.tt_tangent_circles(c1,   &tgt2);
        let q4 = self.tct_tangent_circles(&tgt2, c2);
        (q1, q2, q3, q4)
    }

    fn ttct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.ttct_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, true, self.param);
        let mid2   = HcCcCircle::from_configuration(qc, !c1.left, c1.forward, true, self.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let l1 = cstart.cc_turn_length(&qa) + mid1.hc_turn_length(&qb) + cend.hc_turn_length(&qb);
        let l2 = cstart.cc_turn_length(&qc) + mid2.hc_turn_length(&qd) + cend.hc_turn_length(&qd);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qa), qi2: Some(qb), ci1: Some(mid1), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qc), qi2: Some(qd), ci1: Some(mid2), ..PathSlot::infinite() }
        }
    }

    // ---- TST (TiST / TeST) -----------------------------------------------
    fn tist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward != c2.forward && self.distance >= 2.0 * c1.radius
    }
    fn test_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward != c2.forward
            && self.distance >= 2.0 * c1.radius * c1.sin_mu
    }
    fn tst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tist_exists(c1, c2) || self.test_exists(c1, c2)
    }

    fn tist_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration)
    {
        let dist  = center_distance(c1, c2);
        let angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        let alpha = (2.0 * c1.radius * c1.cos_mu / dist).asin();
        let dx = c1.radius * c1.sin_mu;
        let dy = c1.radius * c1.cos_mu;
        if c1.left && c1.forward {
            let theta = angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx,  dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx,  dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else {
            let theta = angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        }
    }

    fn test_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration)
    {
        let dx    = c1.radius * c1.sin_mu;
        let dy    = c1.radius * c1.cos_mu;
        let theta = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        if c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx,  dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx,  dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        }
    }

    fn tist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (q1, q2) = self.tist_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.cc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn test_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (q1, q2) = self.test_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.cc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tist_exists(c1, c2) { Some(self.tist_path(c1, c2)) }
        else if self.test_exists(c1, c2) { Some(self.test_path(c1, c2)) }
        else { None }
    }

    // ---- TSTcT (TiSTcT / TeSTcT) -----------------------------------------
    fn tistct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left == c2.left && c1.forward == c2.forward
            && self.distance >= ((2.0 * c1.radius * c1.sin_mu + 2.0 * kappa_inv_abs).powi(2)
                + (2.0 * c1.radius * c1.cos_mu).powi(2)).sqrt()
    }
    fn testct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left != c2.left && c1.forward == c2.forward
            && self.distance >= 2.0 * (kappa_inv_abs + c1.radius * c1.sin_mu)
    }
    fn tstct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tistct_exists(c1, c2) || self.testct_exists(c1, c2)
    }

    fn tistct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let delta_y = (4.0 * c2.radius * c2.cos_mu) / (c2.kappa.abs() * self.distance);
        let delta_x = ((2.0 * c2.kappa_inv).powi(2) - delta_y * delta_y).max(0.0).sqrt();
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.param);
        let (q1, q2) = self.tist_tangent_circles(c1, &tgt1);
        let q3 = self.tct_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let ci = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2)
            + ci.hc_turn_length(&q3) + cend.hc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn testct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let delta_x = 2.0 * c2.kappa_inv.abs();
        let delta_y = 0.0;
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.param);
        let (q1, q2) = self.test_tangent_circles(c1, &tgt1);
        let q3 = self.tct_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let ci = HcCcCircle::from_configuration(q2, c1.left, c1.forward, true, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2)
            + ci.hc_turn_length(&q3) + cend.hc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tstct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tistct_exists(c1, c2) { Some(self.tistct_path(c1, c2)) }
        else if self.testct_exists(c1, c2) { Some(self.testct_path(c1, c2)) }
        else { None }
    }

    // ---- TcTST (TcTiST / TcTeST) -----------------------------------------
    fn tctist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left == c2.left && c1.forward == c2.forward
            && self.distance >= ((2.0 * c1.radius * c1.sin_mu + 2.0 * kappa_inv_abs).powi(2)
                + (2.0 * c1.radius * c1.cos_mu).powi(2)).sqrt()
    }
    fn tctest_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left != c2.left && c1.forward == c2.forward
            && self.distance >= 2.0 * (kappa_inv_abs + c1.radius * c1.sin_mu)
    }
    fn tctst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tctist_exists(c1, c2) || self.tctest_exists(c1, c2)
    }

    fn tctist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let delta_y = (4.0 * c2.radius * c2.cos_mu) / (c2.kappa.abs() * self.distance);
        let delta_x = ((2.0 * c2.kappa_inv).powi(2) - delta_y * delta_y).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward, c2.regular, self.param);
        let q1 = self.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.tist_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let ci = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.param);
        let length = cstart.hc_turn_length(&q1) + ci.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3) + cend.cc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tctest_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let delta_x = 2.0 * c2.kappa_inv.abs();
        let (x, y)  = global_frame_change(c1.xc, c1.yc, theta, delta_x, 0.0);
        let tgt1 = HcCcCircle::from_center(x, y, c2.left, !c2.forward, c2.regular, self.param);
        let q1 = self.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.test_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let ci = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.param);
        let length = cstart.hc_turn_length(&q1) + ci.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3) + cend.cc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tctst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tctist_exists(c1, c2) { Some(self.tctist_path(c1, c2)) }
        else if self.tctest_exists(c1, c2) { Some(self.tctest_path(c1, c2)) }
        else { None }
    }

    // ---- TcTSTcT (TcTiSTcT / TcTeSTcT) -----------------------------------
    fn tctistct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left != c2.left && c1.forward != c2.forward
            && self.distance >= ((2.0 * c1.radius).powi(2)
                + 16.0 * c1.radius * c1.sin_mu * kappa_inv_abs
                + (4.0 * c1.kappa_inv).powi(2)).sqrt()
    }
    fn tctestct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left == c2.left && c1.forward != c2.forward
            && self.distance >= 4.0 * kappa_inv_abs + 2.0 * c1.radius * c1.sin_mu
    }
    fn tctcstct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tctistct_exists(c1, c2) || self.tctestct_exists(c1, c2)
    }

    fn tctistct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let delta_y = (4.0 * c1.radius * c1.cos_mu) / (self.distance * c1.kappa.abs());
        let delta_x = ((2.0 * c1.kappa_inv).powi(2) - delta_y * delta_y).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left,  c2.forward, c2.regular, self.param);
        let q1 = self.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.tist_tangent_circles(&tgt1, &tgt2);
        let q4 = self.tct_tangent_circles(&tgt2, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let ci1 = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.param);
        let ci2 = HcCcCircle::from_configuration(q3, !c2.left, c2.forward, true, self.param);
        let length = cstart.hc_turn_length(&q1) + ci1.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3) + ci2.hc_turn_length(&q4)
            + cend.hc_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4),
                   ci1: Some(ci1), ci2: Some(ci2) }
    }

    fn tctestct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let delta_x = 2.0 * c1.kappa_inv.abs();
        let delta_y = 0.0;
        let (x, y)  = global_frame_change(c1.xc, c1.yc, theta,  delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.param);
        let (x, y)  = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left,  c2.forward, c2.regular, self.param);
        let q1 = self.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.test_tangent_circles(&tgt1, &tgt2);
        let q4 = self.tct_tangent_circles(&tgt2, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let ci1 = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.param);
        let ci2 = HcCcCircle::from_configuration(q3, !c2.left, c2.forward, true, self.param);
        let length = cstart.hc_turn_length(&q1) + ci1.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3) + ci2.hc_turn_length(&q4)
            + cend.hc_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4),
                   ci1: Some(ci1), ci2: Some(ci2) }
    }

    fn tctcstct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tctistct_exists(c1, c2) { Some(self.tctistct_path(c1, c2)) }
        else if self.tctestct_exists(c1, c2) { Some(self.tctestct_path(c1, c2)) }
        else { None }
    }

    // ---- TTcTT ------------------------------------------------------------
    fn ttctt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left != c2.left && c1.forward == c2.forward
            && self.distance <= 4.0 * c1.radius + 2.0 * kappa_inv_abs
    }

    fn ttctt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration, Configuration, Configuration)
    {
        let theta        = self.angle;
        let kappa_inv_abs = c1.kappa_inv.abs();
        let r1           = 2.0 * kappa_inv_abs;
        let r2           = 2.0 * c1.radius;
        let delta_x = if self.distance < 4.0 * c1.radius - 2.0 * kappa_inv_abs {
            (self.distance + r1) / 2.0
        } else {
            (self.distance - r1) / 2.0
        };
        let delta_y = (r2*r2 - delta_x*delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x,  delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward, c2.regular, self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x, -delta_y);
        let tgt3 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y);
        let tgt4 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward, c2.regular, self.param);
        let q1 = self.tt_tangent_circles(c1,    &tgt1);
        let q2 = self.tct_tangent_circles(&tgt1, &tgt2);
        let q3 = self.tt_tangent_circles(&tgt2,  c2);
        let q4 = self.tt_tangent_circles(c1,    &tgt3);
        let q5 = self.tct_tangent_circles(&tgt3, &tgt4);
        let q6 = self.tt_tangent_circles(&tgt4,  c2);
        (q1, q2, q3, q4, q5, q6)
    }

    fn ttctt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd, qe, qf) = self.ttctt_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, true, self.param);
        let mid2   = HcCcCircle::from_configuration(qc, !c2.left, c2.forward, true, self.param);
        let mid3   = HcCcCircle::from_configuration(qd, !c1.left, c1.forward, true, self.param);
        let mid4   = HcCcCircle::from_configuration(qf, !c2.left, c2.forward, true, self.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let l1 = cstart.cc_turn_length(&qa) + mid1.hc_turn_length(&qb)
            + mid2.hc_turn_length(&qb) + cend.cc_turn_length(&qc);
        let l2 = cstart.cc_turn_length(&qd) + mid3.hc_turn_length(&qe)
            + mid4.hc_turn_length(&qe) + cend.cc_turn_length(&qf);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qa), qi2: Some(qb), qi3: Some(qc),
                       ci1: Some(mid1), ci2: Some(mid2), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qd), qi2: Some(qe), qi3: Some(qf),
                       ci1: Some(mid3), ci2: Some(mid4), ..PathSlot::infinite() }
        }
    }

    // ---- TcTTcT -----------------------------------------------------------
    fn tctTct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left != c2.left && c1.forward != c2.forward
            && self.distance <= 4.0 * kappa_inv_abs + 2.0 * c1.radius
            && self.distance >= (4.0 * kappa_inv_abs - 2.0 * c1.radius).abs()
    }

    fn tctTct_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration, Configuration, Configuration)
    {
        let theta        = self.angle;
        let kappa_inv_abs = c1.kappa_inv.abs();
        let r1           = 2.0 * kappa_inv_abs;
        let r2           = c1.radius;
        let half_d       = self.distance / 2.0;
        let delta_x      = (r1*r1 + half_d*half_d - r2*r2) / self.distance;
        let delta_y      = (r1*r1 - delta_x*delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left,  c2.forward, c2.regular, self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x, -delta_y);
        let tgt3 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x,  delta_y);
        let tgt4 = HcCcCircle::from_center(x, y, !c2.left,  c2.forward, c2.regular, self.param);
        let q1 = self.tct_tangent_circles(c1,    &tgt1);
        let q2 = self.tt_tangent_circles(&tgt1,  &tgt2);
        let q3 = self.tct_tangent_circles(&tgt2,  c2);
        let q4 = self.tct_tangent_circles(c1,    &tgt3);
        let q5 = self.tt_tangent_circles(&tgt3,  &tgt4);
        let q6 = self.tct_tangent_circles(&tgt4,  c2);
        (q1, q2, q3, q4, q5, q6)
    }

    fn tctTct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd, qe, qf) = self.tctTct_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qb, !c1.left,  c1.forward, true, self.param);
        let mid2   = HcCcCircle::from_configuration(qb,  c1.left, !c1.forward, true, self.param);
        let mid3   = HcCcCircle::from_configuration(qe, !c1.left,  c1.forward, true, self.param);
        let mid4   = HcCcCircle::from_configuration(qe,  c1.left, !c1.forward, true, self.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let l1 = cstart.hc_turn_length(&qa) + mid1.hc_turn_length(&qa)
            + mid2.hc_turn_length(&qc) + cend.hc_turn_length(&qc);
        let l2 = cstart.hc_turn_length(&qd) + mid3.hc_turn_length(&qd)
            + mid4.hc_turn_length(&qf) + cend.hc_turn_length(&qf);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qa), qi2: Some(qb), qi3: Some(qc),
                       ci1: Some(mid1), ci2: Some(mid2), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qd), qi2: Some(qe), qi3: Some(qf),
                       ci1: Some(mid3), ci2: Some(mid4), ..PathSlot::infinite() }
        }
    }

    // ---- TTT --------------------------------------------------------------
    fn ttt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward != c2.forward
            && self.distance <= 4.0 * c1.radius
    }

    fn ttt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration)
    {
        let theta   = self.angle;
        let r       = 2.0 * c1.radius;
        let delta_x = 0.5 * self.distance;
        let delta_y = (r*r - delta_x*delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.param);
        let q1 = self.tt_tangent_circles(c1,   &tgt1);
        let q2 = self.tt_tangent_circles(&tgt1,  c2);
        let q3 = self.tt_tangent_circles(c1,   &tgt2);
        let q4 = self.tt_tangent_circles(&tgt2,  c2);
        (q1, q2, q3, q4)
    }

    fn ttt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.ttt_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, true, self.param);
        let mid2   = HcCcCircle::from_configuration(qc, !c1.left, c1.forward, true, self.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let l1 = cstart.cc_turn_length(&qa) + mid1.cc_turn_length(&qb) + cend.cc_turn_length(&qb);
        let l2 = cstart.cc_turn_length(&qc) + mid2.cc_turn_length(&qd) + cend.cc_turn_length(&qd);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qa), qi2: Some(qb), ci1: Some(mid1), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qc), qi2: Some(qd), ci1: Some(mid2), ..PathSlot::infinite() }
        }
    }

    // ---- TcST (TciST / TceST) --------------------------------------------
    fn tcist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left != c2.left && c1.forward == c2.forward
            && self.distance >= ((c1.radius * c1.sin_mu).powi(2)
                + (c1.radius * c1.cos_mu + kappa_inv_abs).powi(2)).sqrt()
    }
    fn tcest_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left == c2.left && c1.forward == c2.forward
            && self.distance >= ((c1.radius * c1.sin_mu).powi(2)
                + (c1.radius * c1.cos_mu - kappa_inv_abs).powi(2)).sqrt()
    }
    fn tcst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tcist_exists(c1, c2) || self.tcest_exists(c1, c2)
    }

    fn tcist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let kappa_inv_abs = c1.kappa_inv.abs();
        let alpha = ((c1.radius * c1.cos_mu + kappa_inv_abs) / self.distance).asin();
        let dx1 = 0.0;
        let dy1 = kappa_inv_abs;
        let dx2 = c1.radius * c1.sin_mu;
        let dy2 = c1.radius * c1.cos_mu;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1,  dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2,  dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2,  dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1,  dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.cc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcest_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let kappa_inv_abs = c1.kappa_inv.abs();
        let alpha = ((c1.radius * c1.cos_mu - kappa_inv_abs) / self.distance).asin();
        let dx1 = 0.0;
        let dy1 = kappa_inv_abs;
        let dx2 = c1.radius * c1.sin_mu;
        let dy2 = c1.radius * c1.cos_mu;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1,  dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2,  dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1,  dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2,  dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.cc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tcist_exists(c1, c2) { Some(self.tcist_path(c1, c2)) }
        else if self.tcest_exists(c1, c2) { Some(self.tcest_path(c1, c2)) }
        else { None }
    }

    // ---- TScT (TiScT / TeScT) --------------------------------------------
    fn tisct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left != c2.left && c1.forward == c2.forward
            && self.distance >= ((c1.radius * c1.sin_mu).powi(2)
                + (c1.radius * c1.cos_mu + kappa_inv_abs).powi(2)).sqrt()
    }
    fn tesct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left == c2.left && c1.forward == c2.forward
            && self.distance >= ((c1.radius * c1.sin_mu).powi(2)
                + (c1.radius * c1.cos_mu - kappa_inv_abs).powi(2)).sqrt()
    }
    fn tsct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tisct_exists(c1, c2) || self.tesct_exists(c1, c2)
    }

    fn tisct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let kappa_inv_abs = c1.kappa_inv.abs();
        let alpha = ((c1.radius * c1.cos_mu + kappa_inv_abs) / self.distance).asin();
        let dx1 = c1.radius * c1.sin_mu;
        let dy1 = c1.radius * c1.cos_mu;
        let dx2 = 0.0;
        let dy2 = kappa_inv_abs;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1, -dy1);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2,  dy2);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1,  dy1);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2, -dy2);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1,  dy1);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2, -dy2);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1, -dy1);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2,  dy2);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.hc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tesct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let kappa_inv_abs = c1.kappa_inv.abs();
        let alpha = ((c1.radius * c1.cos_mu - kappa_inv_abs) / self.distance).asin();
        let dx1 = c1.radius * c1.sin_mu;
        let dy1 = c1.radius * c1.cos_mu;
        let dx2 = 0.0;
        let dy2 = kappa_inv_abs;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1, -dy1);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2, -dy2);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1,  dy1);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2,  dy2);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1,  dy1);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2,  dy2);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1, -dy1);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2, -dy2);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.hc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tsct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tisct_exists(c1, c2) { Some(self.tisct_path(c1, c2)) }
        else if self.tesct_exists(c1, c2) { Some(self.tesct_path(c1, c2)) }
        else { None }
    }

    // ---- TcScT (TciScT / TceScT) -----------------------------------------
    fn tcisct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        let kappa_inv_abs = c1.kappa_inv.abs();
        c1.left != c2.left && c1.forward != c2.forward
            && self.distance > 2.0 * kappa_inv_abs
    }
    fn tcesct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward != c2.forward
            && self.distance >= get_epsilon()
    }
    fn tcsct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tcisct_exists(c1, c2) || self.tcesct_exists(c1, c2)
    }

    fn tcisct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let kappa_inv_abs = c1.kappa_inv.abs();
        let alpha = (2.0 / (c1.kappa.abs() * self.distance)).asin();
        let dx = 0.0;
        let dy = kappa_inv_abs;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.hc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcesct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.angle;
        let dx = 0.0;
        let dy = c1.kappa_inv.abs();
        let (q1, q2) = if c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.param);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.hc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcsct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tcisct_exists(c1, c2) { Some(self.tcisct_path(c1, c2)) }
        else if self.tcesct_exists(c1, c2) { Some(self.tcesct_path(c1, c2)) }
        else { None }
    }
}

// ---------------------------------------------------------------------------
// HC00 circles → shortest RS path
// ---------------------------------------------------------------------------

fn hc00_circles_rs_path(
    c1: &HcCcCircle,
    c2: &HcCcCircle,
    param: &HcCcCircleParam,
    rs_param: &HcCcCircleParam,
) -> HcCcRsPath {
    use HcCcRsPathType::*;

    let mut h = Hc00RsHelper::new(param, rs_param);
    h.distance = center_distance(c1, c2);
    h.angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);

    let mut slots: Vec<PathSlot> = (0..NB_HC_CC_RS_PATHS).map(|_| PathSlot::infinite()).collect();

    let mut skip = false;

    // case E
    if configuration_equal(&c1.start, &c2.start) {
        slots[E as usize].length = 0.0;
        skip = true;
    }

    // case S
    if !skip {
        if configuration_aligned(&c1.start, &c2.start) {
            slots[S as usize].length = configuration_distance(&c1.start, &c2.start);
            skip = true;
        } else if configuration_aligned(&c2.start, &c1.start) {
            slots[S as usize].length = configuration_distance(&c2.start, &c1.start);
            skip = true;
        }
    }

    // case T
    if !skip && configuration_on_hc_cc_circle(c1, &c2.start) {
        let cs = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, param);
        slots[T as usize].length = cs.cc_turn_length(&c2.start);
        slots[T as usize].cstart = Some(cs);
        skip = true;
    }

    if !skip {
        if h.tt_exists(c1, c2) {
            slots[TT as usize] = h.tt_path(c1, c2);
        }
        if h.tct_exists(c1, c2) {
            slots[TcT as usize] = h.tct_path(c1, c2);
        }
        if h.tctct_exists(c1, c2) {
            slots[TcTcT as usize] = h.tctct_path(c1, c2);
        }
        if h.tctt_exists(c1, c2) {
            slots[TcTT as usize] = h.tctt_path(c1, c2);
        }
        if h.ttct_exists(c1, c2) {
            slots[TTcT as usize] = h.ttct_path(c1, c2);
        }
        if h.tst_exists(c1, c2) {
            if let Some(s) = h.tst_path(c1, c2) {
                slots[TST as usize] = s;
            }
        }
        if h.tstct_exists(c1, c2) {
            if let Some(s) = h.tstct_path(c1, c2) {
                slots[TSTcT as usize] = s;
            }
        }
        if h.tctst_exists(c1, c2) {
            if let Some(s) = h.tctst_path(c1, c2) {
                slots[TcTST as usize] = s;
            }
        }
        if h.tctcstct_exists(c1, c2) {
            if let Some(s) = h.tctcstct_path(c1, c2) {
                slots[TcTSTcT as usize] = s;
            }
        }
        if h.ttctt_exists(c1, c2) {
            slots[TTcTT as usize] = h.ttctt_path(c1, c2);
        }
        if h.tctTct_exists(c1, c2) {
            slots[TcTTcT as usize] = h.tctTct_path(c1, c2);
        }
        if h.ttt_exists(c1, c2) {
            slots[TTT as usize] = h.ttt_path(c1, c2);
        }
        if h.tcst_exists(c1, c2) {
            if let Some(s) = h.tcst_path(c1, c2) {
                slots[TcST as usize] = s;
            }
        }
        if h.tsct_exists(c1, c2) {
            if let Some(s) = h.tsct_path(c1, c2) {
                slots[TScT as usize] = s;
            }
        }
        if h.tcsct_exists(c1, c2) {
            if let Some(s) = h.tcsct_path(c1, c2) {
                slots[TcScT as usize] = s;
            }
        }
    }

    let best = (0..NB_HC_CC_RS_PATHS)
        .min_by(|&a, &b| slots[a].length.partial_cmp(&slots[b].length).unwrap())
        .unwrap();

    let s = &slots[best];
    HcCcRsPath::new(
        c1.start,
        c2.start,
        path_type_from_usize(best),
        param.kappa,
        param.sigma,
        s.qi1,
        s.qi2,
        s.qi3,
        s.qi4,
        s.cstart.clone().map(Box::new),
        s.cend.clone().map(Box::new),
        s.ci1.clone().map(Box::new),
        s.ci2.clone().map(Box::new),
        s.length,
    )
}

// ---------------------------------------------------------------------------
// HC00 Reeds-Shepp state space
// ---------------------------------------------------------------------------

/// HC00 Reeds-Shepp state space — zero curvature at both start and end.
pub struct Hc00RsStateSpace {
    params_: HcCcStateSpaceParams,
    rs_circle_param_: HcCcCircleParam,
    discretization_: f64,
}

impl Hc00RsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64) -> Self {
        let params = HcCcStateSpaceParams::new(kappa, sigma);
        let mut rs_param = HcCcCircleParam::default();
        rs_param.set_param(kappa, f64::MAX, 1.0 / kappa, 0.0, 0.0, 1.0, 0.0);
        Self {
            params_: params,
            rs_circle_param_: rs_param,
            discretization_: discretization,
        }
    }

    fn hc00_reeds_shepp(&self, state1: &State, state2: &State) -> HcCcRsPath {
        let start = Configuration::new(state1.x, state1.y, state1.theta, 0.0);
        let end = Configuration::new(state2.x, state2.y, state2.theta, 0.0);
        let p = &self.params_.hc_cc_circle_param_;

        let start_circles = [
            HcCcCircle::from_configuration(start, true, true, true, p),
            HcCcCircle::from_configuration(start, false, true, true, p),
            HcCcCircle::from_configuration(start, true, false, true, p),
            HcCcCircle::from_configuration(start, false, false, true, p),
        ];
        let end_circles = [
            HcCcCircle::from_configuration(end, true, true, true, p),
            HcCcCircle::from_configuration(end, false, true, true, p),
            HcCcCircle::from_configuration(end, true, false, true, p),
            HcCcCircle::from_configuration(end, false, false, true, p),
        ];

        let mut best: Option<HcCcRsPath> = None;
        for sc in &start_circles {
            for ec in &end_circles {
                let path = hc00_circles_rs_path(sc, ec, p, &self.rs_circle_param_);
                if best.as_ref().map_or(true, |b| path.length < b.length) {
                    best = Some(path);
                }
            }
        }
        best.unwrap()
    }
}

impl StateSpace for Hc00RsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let path = self.hc00_reeds_shepp(s1, s2);
        let mut controls: Vec<Control> = Vec::new();

        match path.path_type {
            HcCcRsPathType::E => {
                empty_controls(&mut controls);
            }
            HcCcRsPathType::S => {
                straight_controls(&path.start, &path.end, &mut controls);
            }
            HcCcRsPathType::T => {
                let cs = path.cstart.as_ref().unwrap();
                cc_turn_controls(cs, &path.end, true, &mut controls);
            }
            HcCcRsPathType::TT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                cc_turn_controls(ce, q1, false, &mut controls);
            }
            HcCcRsPathType::TcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ce, q1, false, &mut controls);
            }
            HcCcRsPathType::TcTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                rs_turn_controls(ci, q2, true, &mut controls);
                hc_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TcTT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci, q1, false, &mut controls);
                cc_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci, q2, true, &mut controls);
                hc_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TST => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                straight_controls(q1, q2, &mut controls);
                cc_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TSTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                straight_controls(q1, q2, &mut controls);
                hc_turn_controls(ci, q3, true, &mut controls);
                hc_turn_controls(ce, q3, false, &mut controls);
            }
            HcCcRsPathType::TcTST => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci, q1, false, &mut controls);
                straight_controls(q2, q3, &mut controls);
                cc_turn_controls(ce, q3, false, &mut controls);
            }
            HcCcRsPathType::TcTSTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci1 = path.ci1.as_ref().unwrap();
                let ci2 = path.ci2.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                let q4 = path.qi4.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci1, q1, false, &mut controls);
                straight_controls(q2, q3, &mut controls);
                hc_turn_controls(ci2, q4, true, &mut controls);
                hc_turn_controls(ce, q4, false, &mut controls);
            }
            HcCcRsPathType::TTcTT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci1 = path.ci1.as_ref().unwrap();
                let ci2 = path.ci2.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci1, q2, true, &mut controls);
                hc_turn_controls(ci2, q2, false, &mut controls);
                cc_turn_controls(ce, q3, false, &mut controls);
            }
            HcCcRsPathType::TcTTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci1 = path.ci1.as_ref().unwrap();
                let ci2 = path.ci2.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci1, q1, false, &mut controls);
                hc_turn_controls(ci2, q3, true, &mut controls);
                hc_turn_controls(ce, q3, false, &mut controls);
            }
            HcCcRsPathType::TTT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                cc_turn_controls(ci, q2, true, &mut controls);
                cc_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TcST => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                straight_controls(q1, q2, &mut controls);
                cc_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TScT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                straight_controls(q1, q2, &mut controls);
                hc_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TcScT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                straight_controls(q1, q2, &mut controls);
                hc_turn_controls(ce, q2, false, &mut controls);
            }
        }

        controls
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        vec![self.get_controls(s1, s2)]
    }

    fn discretization(&self) -> f64 {
        self.discretization_
    }
}

// ---------------------------------------------------------------------------
// HC0pm geometry helper  (delegates _exists / tangent methods to HC00)
// ---------------------------------------------------------------------------

struct Hc0pmRsHelper<'a> {
    base: Hc00RsHelper<'a>,
}

#[allow(non_snake_case)]
impl<'a> Hc0pmRsHelper<'a> {
    fn new(param: &'a HcCcCircleParam, rs_param: &'a HcCcCircleParam) -> Self {
        Self { base: Hc00RsHelper::new(param, rs_param) }
    }

    // ---- TT ---------------------------------------------------------------
    fn tt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let q1 = self.base.tt_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.base.param);
        let cend = HcCcCircle::from_configuration(q1, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q2 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let length = cstart.cc_turn_length(&q1) + cend.hc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    // ---- TcT --------------------------------------------------------------
    fn tct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let q = self.base.tct_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let length = cstart.hc_turn_length(&q) + cend.rs_turn_length(&q);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q), ..PathSlot::infinite() }
    }

    // ---- TcTcT ------------------------------------------------------------
    fn tctct_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration)
    {
        let theta   = self.base.angle;
        let r       = 2.0 * c1.kappa_inv.abs();
        let delta_x = 0.5 * self.base.distance;
        let delta_y = (r * r - delta_x * delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.base.rs_param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.base.rs_param);
        let q1 = self.base.tct_tangent_circles(c1,    &tgt1);
        let q2 = self.base.tct_tangent_circles(&tgt1,  c2);
        let q3 = self.base.tct_tangent_circles(c1,    &tgt2);
        let q4 = self.base.tct_tangent_circles(&tgt2,  c2);
        (q1, q2, q3, q4)
    }

    fn tctct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.tctct_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, !c1.forward, true, self.base.rs_param);
        let mid2   = HcCcCircle::from_configuration(qc, !c1.left, !c1.forward, true, self.base.rs_param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let l1 = cstart.hc_turn_length(&qa) + mid1.rs_turn_length(&qb) + cend.rs_turn_length(&qb);
        let l2 = cstart.hc_turn_length(&qc) + mid2.rs_turn_length(&qd) + cend.rs_turn_length(&qd);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qa), qi2: Some(qb), ci1: Some(mid1), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qc), qi2: Some(qd), ci1: Some(mid2), ..PathSlot::infinite() }
        }
    }

    // ---- TcTT -------------------------------------------------------------
    fn tctt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.base.tctt_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qb, !c1.left, c1.forward, true, self.base.param);
        let mid2   = HcCcCircle::from_configuration(qd, !c1.left, c1.forward, true, self.base.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let q2_end = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let l1 = {
            let end_circ = HcCcCircle::from_configuration(qb, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            cstart.hc_turn_length(&qa) + mid1.hc_turn_length(&qa) + end_circ.hc_turn_length(&q2_end)
        };
        let l2 = {
            let end_circ = HcCcCircle::from_configuration(qd, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            cstart.hc_turn_length(&qc) + mid2.hc_turn_length(&qc) + end_circ.hc_turn_length(&q2_end)
        };
        if l1 <= l2 {
            let end_circ = HcCcCircle::from_configuration(qb, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(end_circ),
                       qi1: Some(qa), qi2: Some(q2_end), ci1: Some(mid1), ..PathSlot::infinite() }
        } else {
            let end_circ = HcCcCircle::from_configuration(qd, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(end_circ),
                       qi1: Some(qc), qi2: Some(q2_end), ci1: Some(mid2), ..PathSlot::infinite() }
        }
    }

    // ---- TTcT -------------------------------------------------------------
    fn ttct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.base.ttct_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, true, self.base.param);
        let mid2   = HcCcCircle::from_configuration(qc, !c1.left, c1.forward, true, self.base.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let l1 = cstart.cc_turn_length(&qa) + mid1.hc_turn_length(&qb) + cend.rs_turn_length(&qb);
        let l2 = cstart.cc_turn_length(&qc) + mid2.hc_turn_length(&qd) + cend.rs_turn_length(&qd);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qa), qi2: Some(qb), ci1: Some(mid1), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qc), qi2: Some(qd), ci1: Some(mid2), ..PathSlot::infinite() }
        }
    }

    // ---- TST (TiST / TeST) -----------------------------------------------
    fn tist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (q1, q2) = self.base.tist_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.base.param);
        let cend = HcCcCircle::from_configuration(q2, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q3 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.hc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ..PathSlot::infinite() }
    }

    fn test_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (q1, q2) = self.base.test_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.base.param);
        let cend = HcCcCircle::from_configuration(q2, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q3 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.hc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ..PathSlot::infinite() }
    }

    fn tst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.base.tist_exists(c1, c2) { Some(self.tist_path(c1, c2)) }
        else if self.base.test_exists(c1, c2) { Some(self.test_path(c1, c2)) }
        else { None }
    }

    // ---- TSTcT (TiSTcT / TeSTcT) -----------------------------------------
    fn tistct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.base.angle;
        let delta_y = (4.0 * c2.radius * c2.cos_mu) / (c2.kappa.abs() * self.base.distance);
        let delta_x = ((2.0 * c2.kappa_inv).powi(2) - delta_y * delta_y).max(0.0).sqrt();
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.base.param);
        let (q1, q2) = self.base.tist_tangent_circles(c1, &tgt1);
        let q3 = self.base.tct_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let ci = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.base.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2)
            + ci.hc_turn_length(&q3) + cend.rs_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn testct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.base.angle;
        let delta_x = 2.0 * c2.kappa_inv.abs();
        let delta_y = 0.0;
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.base.param);
        let (q1, q2) = self.base.test_tangent_circles(c1, &tgt1);
        let q3 = self.base.tct_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let ci = HcCcCircle::from_configuration(q2, c1.left, c1.forward, true, self.base.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2)
            + ci.hc_turn_length(&q3) + cend.rs_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tstct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.base.tistct_exists(c1, c2) { Some(self.tistct_path(c1, c2)) }
        else if self.base.testct_exists(c1, c2) { Some(self.testct_path(c1, c2)) }
        else { None }
    }

    // ---- TcTST (TcTiST / TcTeST) -----------------------------------------
    fn tctist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.base.angle;
        let delta_y = (4.0 * c2.radius * c2.cos_mu) / (c2.kappa.abs() * self.base.distance);
        let delta_x = ((2.0 * c2.kappa_inv).powi(2) - delta_y * delta_y).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward, c2.regular, self.base.param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.base.tist_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend_circ = HcCcCircle::from_configuration(q3, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q4 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let ci = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.base.param);
        let length = cstart.hc_turn_length(&q1) + ci.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3) + cend_circ.hc_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend_circ),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4),
                   ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tctest_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.base.angle;
        let delta_x = 2.0 * c2.kappa_inv.abs();
        let (x, y)  = global_frame_change(c1.xc, c1.yc, theta, delta_x, 0.0);
        let tgt1 = HcCcCircle::from_center(x, y, c2.left, !c2.forward, c2.regular, self.base.param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.base.test_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend_circ = HcCcCircle::from_configuration(q3, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q4 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let ci = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.base.param);
        let length = cstart.hc_turn_length(&q1) + ci.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3) + cend_circ.hc_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend_circ),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4),
                   ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tctst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.base.tctist_exists(c1, c2) { Some(self.tctist_path(c1, c2)) }
        else if self.base.tctest_exists(c1, c2) { Some(self.tctest_path(c1, c2)) }
        else { None }
    }

    // ---- TcTSTcT (TcTiSTcT / TcTeSTcT) -----------------------------------
    fn tctistct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.base.angle;
        let delta_y = (4.0 * c1.radius * c1.cos_mu) / (self.base.distance * c1.kappa.abs());
        let delta_x = ((2.0 * c1.kappa_inv).powi(2) - delta_y * delta_y).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.base.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left,  c2.forward, c2.regular, self.base.param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.base.tist_tangent_circles(&tgt1, &tgt2);
        let q4 = self.base.tct_tangent_circles(&tgt2, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let ci1 = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.base.param);
        let ci2 = HcCcCircle::from_configuration(q3, !c2.left, c2.forward, true, self.base.param);
        let length = cstart.hc_turn_length(&q1) + ci1.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3) + ci2.hc_turn_length(&q4)
            + cend.rs_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4),
                   ci1: Some(ci1), ci2: Some(ci2) }
    }

    fn tctestct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.base.angle;
        let delta_x = 2.0 * c1.kappa_inv.abs();
        let delta_y = 0.0;
        let (x, y)  = global_frame_change(c1.xc, c1.yc, theta,  delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.base.param);
        let (x, y)  = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left,  c2.forward, c2.regular, self.base.param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.base.test_tangent_circles(&tgt1, &tgt2);
        let q4 = self.base.tct_tangent_circles(&tgt2, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let ci1 = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.base.param);
        let ci2 = HcCcCircle::from_configuration(q3, !c2.left, c2.forward, true, self.base.param);
        let length = cstart.hc_turn_length(&q1) + ci1.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3) + ci2.hc_turn_length(&q4)
            + cend.rs_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4),
                   ci1: Some(ci1), ci2: Some(ci2) }
    }

    fn tctcstct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.base.tctistct_exists(c1, c2) { Some(self.tctistct_path(c1, c2)) }
        else if self.base.tctestct_exists(c1, c2) { Some(self.tctestct_path(c1, c2)) }
        else { None }
    }

    // ---- TTcTT ------------------------------------------------------------
    fn ttctt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd, qe, qf) = self.base.ttctt_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, true, self.base.param);
        let mid2   = HcCcCircle::from_configuration(qc, !c2.left, c2.forward, true, self.base.param);
        let mid3   = HcCcCircle::from_configuration(qd, !c1.left, c1.forward, true, self.base.param);
        let mid4   = HcCcCircle::from_configuration(qf, !c2.left, c2.forward, true, self.base.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.base.param);
        let q3_end = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let l1 = {
            let end_circ = HcCcCircle::from_configuration(qc, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            cstart.cc_turn_length(&qa) + mid1.hc_turn_length(&qb) + mid2.hc_turn_length(&qb) + end_circ.hc_turn_length(&q3_end)
        };
        let l2 = {
            let end_circ = HcCcCircle::from_configuration(qf, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            cstart.cc_turn_length(&qd) + mid3.hc_turn_length(&qe) + mid4.hc_turn_length(&qe) + end_circ.hc_turn_length(&q3_end)
        };
        if l1 <= l2 {
            let end_circ = HcCcCircle::from_configuration(qc, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(end_circ),
                       qi1: Some(qa), qi2: Some(qb), qi3: Some(q3_end),
                       ci1: Some(mid1), ci2: Some(mid2), ..PathSlot::infinite() }
        } else {
            let end_circ = HcCcCircle::from_configuration(qf, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(end_circ),
                       qi1: Some(qd), qi2: Some(qe), qi3: Some(q3_end),
                       ci1: Some(mid3), ci2: Some(mid4), ..PathSlot::infinite() }
        }
    }

    // ---- TcTTcT -----------------------------------------------------------
    fn tctTct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd, qe, qf) = self.base.tctTct_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qb, !c1.left,  c1.forward, true, self.base.param);
        let mid2   = HcCcCircle::from_configuration(qb,  c1.left, !c1.forward, true, self.base.param);
        let mid3   = HcCcCircle::from_configuration(qe, !c1.left,  c1.forward, true, self.base.param);
        let mid4   = HcCcCircle::from_configuration(qe,  c1.left, !c1.forward, true, self.base.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let l1 = cstart.hc_turn_length(&qa) + mid1.hc_turn_length(&qa)
            + mid2.hc_turn_length(&qc) + cend.rs_turn_length(&qc);
        let l2 = cstart.hc_turn_length(&qd) + mid3.hc_turn_length(&qd)
            + mid4.hc_turn_length(&qf) + cend.rs_turn_length(&qf);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qa), qi2: Some(qc),
                       ci1: Some(mid1), ci2: Some(mid2), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend),
                       qi1: Some(qd), qi2: Some(qf),
                       ci1: Some(mid3), ci2: Some(mid4), ..PathSlot::infinite() }
        }
    }

    // ---- TTT --------------------------------------------------------------
    fn ttt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.base.ttt_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, true, self.base.param);
        let mid2   = HcCcCircle::from_configuration(qc, !c1.left, c1.forward, true, self.base.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.base.param);
        let q3_end1 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let l1 = {
            let end_circ = HcCcCircle::from_configuration(qb, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            cstart.cc_turn_length(&qa) + mid1.cc_turn_length(&qb) + end_circ.hc_turn_length(&q3_end1)
        };
        let l2 = {
            let end_circ = HcCcCircle::from_configuration(qd, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            cstart.cc_turn_length(&qc) + mid2.cc_turn_length(&qd) + end_circ.hc_turn_length(&q3_end1)
        };
        if l1 <= l2 {
            let end_circ = HcCcCircle::from_configuration(qb, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(end_circ),
                       qi1: Some(qa), qi2: Some(qb), qi3: Some(q3_end1),
                       ci1: Some(mid1), ..PathSlot::infinite() }
        } else {
            let end_circ = HcCcCircle::from_configuration(qd, c2.left, !c2.forward, HC_REGULAR, self.base.param);
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(end_circ),
                       qi1: Some(qc), qi2: Some(qd), qi3: Some(q3_end1),
                       ci1: Some(mid2), ..PathSlot::infinite() }
        }
    }

    // ---- TcST (TciST / TceST) --------------------------------------------
    fn tcist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let kappa_inv_abs = c1.kappa_inv.abs();
        let alpha = ((c1.radius * c1.cos_mu + kappa_inv_abs) / self.base.distance).asin();
        let dx1 = 0.0;
        let dy1 = kappa_inv_abs;
        let dx2 = c1.radius * c1.sin_mu;
        let dy2 = c1.radius * c1.cos_mu;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1,  dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2,  dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2,  dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1,  dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend = HcCcCircle::from_configuration(q2, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q3 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.hc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ..PathSlot::infinite() }
    }

    fn tcest_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let kappa_inv_abs = c1.kappa_inv.abs();
        let alpha = ((c1.radius * c1.cos_mu - kappa_inv_abs) / self.base.distance).asin();
        let dx1 = 0.0;
        let dy1 = kappa_inv_abs;
        let dx2 = c1.radius * c1.sin_mu;
        let dy2 = c1.radius * c1.cos_mu;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1,  dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2,  dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1,  dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2,  dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend = HcCcCircle::from_configuration(q2, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q3 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.hc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ..PathSlot::infinite() }
    }

    fn tcst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.base.tcist_exists(c1, c2) { Some(self.tcist_path(c1, c2)) }
        else if self.base.tcest_exists(c1, c2) { Some(self.tcest_path(c1, c2)) }
        else { None }
    }

    // ---- TScT (TiScT / TeScT) --------------------------------------------
    fn tisct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let kappa_inv_abs = c1.kappa_inv.abs();
        let alpha = ((c1.radius * c1.cos_mu + kappa_inv_abs) / self.base.distance).asin();
        let dx1 = c1.radius * c1.sin_mu;
        let dy1 = c1.radius * c1.cos_mu;
        let dx2 = 0.0;
        let dy2 = kappa_inv_abs;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1, -dy1);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2,  dy2);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1,  dy1);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2, -dy2);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1,  dy1);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2, -dy2);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1, -dy1);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2,  dy2);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.rs_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tesct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let kappa_inv_abs = c1.kappa_inv.abs();
        let alpha = ((c1.radius * c1.cos_mu - kappa_inv_abs) / self.base.distance).asin();
        let dx1 = c1.radius * c1.sin_mu;
        let dy1 = c1.radius * c1.cos_mu;
        let dx2 = 0.0;
        let dy2 = kappa_inv_abs;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1, -dy1);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2, -dy2);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1,  dy1);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2,  dy2);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1,  dy1);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2,  dy2);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx1, -dy1);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx2, -dy2);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.rs_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tsct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.base.tisct_exists(c1, c2) { Some(self.tisct_path(c1, c2)) }
        else if self.base.tesct_exists(c1, c2) { Some(self.tesct_path(c1, c2)) }
        else { None }
    }

    // ---- TcScT (TciScT / TceScT) -----------------------------------------
    fn tcisct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let kappa_inv_abs = c1.kappa_inv.abs();
        let alpha = (2.0 / (c1.kappa.abs() * self.base.distance)).asin();
        let dx = 0.0;
        let dy = kappa_inv_abs;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.rs_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcesct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.base.angle;
        let dx = 0.0;
        let dy = c1.kappa_inv.abs();
        let (q1, q2) = if c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.rs_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcsct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.base.tcisct_exists(c1, c2) { Some(self.tcisct_path(c1, c2)) }
        else if self.base.tcesct_exists(c1, c2) { Some(self.tcesct_path(c1, c2)) }
        else { None }
    }
}

// ---------------------------------------------------------------------------
// HC0pm circles → shortest RS path
// ---------------------------------------------------------------------------

fn hc0pm_circles_rs_path(
    c1: &HcCcCircle,
    c2: &HcCcCircle,
    param: &HcCcCircleParam,
    rs_param: &HcCcCircleParam,
) -> HcCcRsPath {
    use HcCcRsPathType::*;

    let mut h = Hc0pmRsHelper::new(param, rs_param);
    h.base.distance = center_distance(c1, c2);
    h.base.angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);

    let mut slots: Vec<PathSlot> = (0..NB_HC_CC_RS_PATHS).map(|_| PathSlot::infinite()).collect();

    let mut skip = false;

    // case E
    if configuration_equal(&c1.start, &c2.start) {
        slots[E as usize].length = 0.0;
        skip = true;
    }

    // case T
    if !skip && h.base.distance < get_epsilon() {
        let cs = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, HC_REGULAR, param);
        slots[T as usize].length = cs.hc_turn_length(&c2.start);
        slots[T as usize].cstart = Some(cs);
        skip = true;
    }

    if !skip {
        if h.base.tt_exists(c1, c2) {
            slots[TT as usize] = h.tt_path(c1, c2);
        }
        if h.base.tct_exists(c1, c2) {
            slots[TcT as usize] = h.tct_path(c1, c2);
        }
        if h.base.tctct_exists(c1, c2) {
            slots[TcTcT as usize] = h.tctct_path(c1, c2);
        }
        if h.base.tctt_exists(c1, c2) {
            slots[TcTT as usize] = h.tctt_path(c1, c2);
        }
        if h.base.ttct_exists(c1, c2) {
            slots[TTcT as usize] = h.ttct_path(c1, c2);
        }
        if h.base.tst_exists(c1, c2) {
            if let Some(s) = h.tst_path(c1, c2) {
                slots[TST as usize] = s;
            }
        }
        if h.base.tstct_exists(c1, c2) {
            if let Some(s) = h.tstct_path(c1, c2) {
                slots[TSTcT as usize] = s;
            }
        }
        if h.base.tctst_exists(c1, c2) {
            if let Some(s) = h.tctst_path(c1, c2) {
                slots[TcTST as usize] = s;
            }
        }
        if h.base.tctcstct_exists(c1, c2) {
            if let Some(s) = h.tctcstct_path(c1, c2) {
                slots[TcTSTcT as usize] = s;
            }
        }
        if h.base.ttctt_exists(c1, c2) {
            slots[TTcTT as usize] = h.ttctt_path(c1, c2);
        }
        if h.base.tctTct_exists(c1, c2) {
            slots[TcTTcT as usize] = h.tctTct_path(c1, c2);
        }
        if h.base.ttt_exists(c1, c2) {
            slots[TTT as usize] = h.ttt_path(c1, c2);
        }
        if h.base.tcst_exists(c1, c2) {
            if let Some(s) = h.tcst_path(c1, c2) {
                slots[TcST as usize] = s;
            }
        }
        if h.base.tsct_exists(c1, c2) {
            if let Some(s) = h.tsct_path(c1, c2) {
                slots[TScT as usize] = s;
            }
        }
        if h.base.tcsct_exists(c1, c2) {
            if let Some(s) = h.tcsct_path(c1, c2) {
                slots[TcScT as usize] = s;
            }
        }
    }

    let best = (0..NB_HC_CC_RS_PATHS)
        .min_by(|&a, &b| slots[a].length.partial_cmp(&slots[b].length).unwrap())
        .unwrap();

    let s = &slots[best];
    HcCcRsPath::new(
        c1.start,
        c2.start,
        path_type_from_usize(best),
        param.kappa,
        param.sigma,
        s.qi1,
        s.qi2,
        s.qi3,
        s.qi4,
        s.cstart.clone().map(Box::new),
        s.cend.clone().map(Box::new),
        s.ci1.clone().map(Box::new),
        s.ci2.clone().map(Box::new),
        s.length,
    )
}

// ---------------------------------------------------------------------------
// HC0pm Reeds-Shepp state space
// ---------------------------------------------------------------------------

/// HC0pm Reeds-Shepp state space — zero curvature at start, ±kappa at end.
pub struct Hc0pmRsStateSpace {
    params_: HcCcStateSpaceParams,
    rs_circle_param_: HcCcCircleParam,
    discretization_: f64,
}

impl Hc0pmRsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64) -> Self {
        let params = HcCcStateSpaceParams::new(kappa, sigma);
        let mut rs_param = HcCcCircleParam::default();
        rs_param.set_param(kappa, f64::MAX, 1.0 / kappa, 0.0, 0.0, 1.0, 0.0);
        Self {
            params_: params,
            rs_circle_param_: rs_param,
            discretization_: discretization,
        }
    }

    fn hc0pm_reeds_shepp(&self, state1: &State, state2: &State) -> HcCcRsPath {
        let kappa = self.params_.kappa_;
        let p   = &self.params_.hc_cc_circle_param_;
        let rsp = &self.rs_circle_param_;

        let start = Configuration::new(state1.x, state1.y, state1.theta, 0.0);
        let end1  = Configuration::new(state2.x, state2.y, state2.theta,  kappa);
        let end2  = Configuration::new(state2.x, state2.y, state2.theta, -kappa);

        let start_circles = [
            HcCcCircle::from_configuration(start,  true,  true, true, p),
            HcCcCircle::from_configuration(start, false,  true, true, p),
            HcCcCircle::from_configuration(start,  true, false, true, p),
            HcCcCircle::from_configuration(start, false, false, true, p),
        ];
        let end_circles = [
            HcCcCircle::from_configuration(end1,  true,  true, true, rsp),
            HcCcCircle::from_configuration(end2, false,  true, true, rsp),
            HcCcCircle::from_configuration(end1,  true, false, true, rsp),
            HcCcCircle::from_configuration(end2, false, false, true, rsp),
        ];

        let mut best: Option<HcCcRsPath> = None;
        for (i, sc) in start_circles.iter().enumerate() {
            for (j, ec) in end_circles.iter().enumerate() {
                // j=0,2 use end1 (kappa>0) → skip if s2.kappa < 0
                // j=1,3 use end2 (kappa<0) → skip if s2.kappa > 0
                if (j == 0 || j == 2) && state2.kappa < 0.0 { continue; }
                if (j == 1 || j == 3) && state2.kappa > 0.0 { continue; }
                let _ = i; // suppress unused warning
                let path = hc0pm_circles_rs_path(sc, ec, p, rsp);
                if best.as_ref().map_or(true, |b| path.length < b.length) {
                    best = Some(path);
                }
            }
        }
        best.unwrap()
    }
}

impl StateSpace for Hc0pmRsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let path = self.hc0pm_reeds_shepp(s1, s2);
        let mut controls: Vec<Control> = Vec::new();

        match path.path_type {
            HcCcRsPathType::E => {
                empty_controls(&mut controls);
            }
            HcCcRsPathType::S => {
                straight_controls(&path.start, &path.end, &mut controls);
            }
            HcCcRsPathType::T => {
                let cs = path.cstart.as_ref().unwrap();
                hc_turn_controls(cs, &path.end, true, &mut controls);
            }
            HcCcRsPathType::TT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ce, q2, true, &mut controls);
            }
            HcCcRsPathType::TcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                rs_turn_controls(ce, q1, false, &mut controls);
            }
            HcCcRsPathType::TcTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                rs_turn_controls(ci, q2, true, &mut controls);
                rs_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TcTT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci, q1, false, &mut controls);
                hc_turn_controls(ce, q2, true, &mut controls);
            }
            HcCcRsPathType::TTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci, q2, true, &mut controls);
                rs_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TST => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                straight_controls(q1, q2, &mut controls);
                hc_turn_controls(ce, q3, true, &mut controls);
            }
            HcCcRsPathType::TSTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                straight_controls(q1, q2, &mut controls);
                hc_turn_controls(ci, q3, true, &mut controls);
                rs_turn_controls(ce, q3, false, &mut controls);
            }
            HcCcRsPathType::TcTST => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                let q4 = path.qi4.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci, q1, false, &mut controls);
                straight_controls(q2, q3, &mut controls);
                hc_turn_controls(ce, q4, true, &mut controls);
            }
            HcCcRsPathType::TcTSTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci1 = path.ci1.as_ref().unwrap();
                let ci2 = path.ci2.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                let q4 = path.qi4.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci1, q1, false, &mut controls);
                straight_controls(q2, q3, &mut controls);
                hc_turn_controls(ci2, q4, true, &mut controls);
                rs_turn_controls(ce, q4, false, &mut controls);
            }
            HcCcRsPathType::TTcTT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci1 = path.ci1.as_ref().unwrap();
                let ci2 = path.ci2.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci1, q2, true, &mut controls);
                hc_turn_controls(ci2, q2, false, &mut controls);
                hc_turn_controls(ce, q3, true, &mut controls);
            }
            HcCcRsPathType::TcTTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci1 = path.ci1.as_ref().unwrap();
                let ci2 = path.ci2.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                hc_turn_controls(ci1, q1, false, &mut controls);
                hc_turn_controls(ci2, q2, true, &mut controls);
                rs_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TTT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                cc_turn_controls(ci, q2, true, &mut controls);
                hc_turn_controls(ce, q3, true, &mut controls);
            }
            HcCcRsPathType::TcST => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                straight_controls(q1, q2, &mut controls);
                hc_turn_controls(ce, q3, true, &mut controls);
            }
            HcCcRsPathType::TScT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true, &mut controls);
                straight_controls(q1, q2, &mut controls);
                rs_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TcScT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                hc_turn_controls(cs, q1, true, &mut controls);
                straight_controls(q1, q2, &mut controls);
                rs_turn_controls(ce, q2, false, &mut controls);
            }
        }

        controls
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        vec![self.get_controls(s1, s2)]
    }

    fn discretization(&self) -> f64 {
        self.discretization_
    }
}

fn state_to_configuration(state: &State) -> Configuration {
    Configuration::new(state.x, state.y, state.theta, state.kappa)
}

fn controls_length(controls: &[Control]) -> f64 {
    controls.iter().map(|control| control.delta_s.abs()).sum()
}

fn reverse_controls(controls: &[Control]) -> Vec<Control> {
    let mut reversed = Vec::with_capacity(controls.len());
    for control in controls.iter().rev() {
        let mut reversed_control = *control;
        reverse_control(&mut reversed_control);
        reversed.push(reversed_control);
    }
    reversed
}

pub struct Hcpm0RsStateSpace {
    inner_: Hc0pmRsStateSpace,
    discretization_: f64,
}

impl Hcpm0RsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64) -> Self {
        Self {
            inner_: Hc0pmRsStateSpace::new(kappa, sigma, discretization),
            discretization_: discretization,
        }
    }
}

impl StateSpace for Hcpm0RsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        reverse_controls(&self.inner_.get_controls(s2, s1))
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        vec![self.get_controls(s1, s2)]
    }

    fn discretization(&self) -> f64 {
        self.discretization_
    }
}

struct HcpmpmRsHelper<'a> {
    base: Hc00RsHelper<'a>,
    radius_: f64,
    mu_: f64,
    sin_mu_: f64,
    cos_mu_: f64,
}

#[allow(non_snake_case)]
impl<'a> HcpmpmRsHelper<'a> {
    fn new(param: &'a HcCcCircleParam, rs_param: &'a HcCcCircleParam) -> Self {
        Self {
            base: Hc00RsHelper::new(param, rs_param),
            radius_: param.radius,
            mu_: param.mu,
            sin_mu_: param.sin_mu,
            cos_mu_: param.cos_mu,
        }
    }

    fn tt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left
            && c1.forward != c2.forward
            && (self.base.distance - 2.0 * self.radius_).abs() < get_epsilon()
    }

    fn tt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Configuration {
        let x = (c1.xc + c2.xc) / 2.0;
        let y = (c1.yc + c2.yc) / 2.0;
        let angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        let theta = if c1.left {
            if c1.forward { angle + HALF_PI - self.mu_ } else { angle + HALF_PI + self.mu_ }
        } else if c1.forward {
            angle - HALF_PI + self.mu_
        } else {
            angle - HALF_PI - self.mu_
        };
        Configuration::new(x, y, theta, 0.0)
    }

    fn tt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let q2 = self.tt_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let cend = HcCcCircle::from_configuration(q2, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q1 = Configuration::new(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
        let q3 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let length = cstart.hc_turn_length(&q1) + cend.hc_turn_length(&q3);
        PathSlot {
            length,
            cstart: Some(cstart),
            cend: Some(cend),
            qi1: Some(q1),
            qi2: Some(q2),
            qi3: Some(q3),
            ..PathSlot::infinite()
        }
    }

    fn tct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let q = self.base.tct_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let length = cstart.rs_turn_length(&q) + cend.rs_turn_length(&q);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q), ..PathSlot::infinite() }
    }

    fn tctct_tangent_circles(
        &self,
        c1: &HcCcCircle,
        c2: &HcCcCircle,
    ) -> (Configuration, Configuration, Configuration, Configuration) {
        let theta = self.base.angle;
        let r = 2.0 * c1.kappa_inv.abs();
        let delta_x = 0.5 * self.base.distance;
        let delta_y = (r * r - delta_x * delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, true, self.base.rs_param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, true, self.base.rs_param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let q2 = self.base.tct_tangent_circles(&tgt1, c2);
        let q3 = self.base.tct_tangent_circles(c1, &tgt2);
        let q4 = self.base.tct_tangent_circles(&tgt2, c2);
        (q1, q2, q3, q4)
    }

    fn tctct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.tctct_tangent_circles(c1, c2);
        let middle1 = HcCcCircle::from_configuration(qa, !c1.left, !c1.forward, true, self.base.rs_param);
        let middle2 = HcCcCircle::from_configuration(qc, !c1.left, !c1.forward, true, self.base.rs_param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let l1 = cstart.rs_turn_length(&qa) + middle1.rs_turn_length(&qb) + cend.rs_turn_length(&qb);
        let l2 = cstart.rs_turn_length(&qc) + middle2.rs_turn_length(&qd) + cend.rs_turn_length(&qd);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend), qi1: Some(qa), qi2: Some(qb), ci1: Some(middle1), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend), qi1: Some(qc), qi2: Some(qd), ci1: Some(middle2), ..PathSlot::infinite() }
        }
    }

    fn tctt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left
            && c1.forward == c2.forward
            && self.base.distance <= 2.0 * self.radius_ + 2.0 * c1.kappa_inv.abs()
            && self.base.distance >= 2.0 * self.radius_ - 2.0 * c1.kappa_inv.abs()
    }

    fn tctt_tangent_circles(
        &self,
        c1: &HcCcCircle,
        c2: &HcCcCircle,
    ) -> (Configuration, Configuration, Configuration, Configuration) {
        let theta = self.base.angle;
        let r1 = 2.0 * c1.kappa_inv.abs();
        let r2 = 2.0 * self.radius_;
        let delta_x = (r1 * r1 + self.base.distance * self.base.distance - r2 * r2) / (2.0 * self.base.distance);
        let delta_y = (r1 * r1 - delta_x * delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.base.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.base.param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let q2 = self.tt_tangent_circles(&tgt1, c2);
        let q3 = self.base.tct_tangent_circles(c1, &tgt2);
        let q4 = self.tt_tangent_circles(&tgt2, c2);
        (q1, q2, q3, q4)
    }

    fn tctt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.tctt_tangent_circles(c1, c2);
        let end1 = HcCcCircle::from_configuration(qb, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let end2 = HcCcCircle::from_configuration(qd, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let middle1 = HcCcCircle::from_configuration(qb, !c1.left, c1.forward, true, self.base.param);
        let middle2 = HcCcCircle::from_configuration(qd, !c1.left, c1.forward, true, self.base.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let q2 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let l1 = cstart.rs_turn_length(&qa) + middle1.hc_turn_length(&qa) + end1.hc_turn_length(&q2);
        let l2 = cstart.rs_turn_length(&qc) + middle2.hc_turn_length(&qc) + end2.hc_turn_length(&q2);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(end1), qi1: Some(qa), qi2: Some(q2), ci1: Some(middle1), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(end2), qi1: Some(qc), qi2: Some(q2), ci1: Some(middle2), ..PathSlot::infinite() }
        }
    }

    fn ttct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left
            && c1.forward != c2.forward
            && self.base.distance <= 2.0 * self.radius_ + 2.0 * c1.kappa_inv.abs()
            && self.base.distance >= 2.0 * self.radius_ - 2.0 * c1.kappa_inv.abs()
    }

    fn ttct_tangent_circles(
        &self,
        c1: &HcCcCircle,
        c2: &HcCcCircle,
    ) -> (Configuration, Configuration, Configuration, Configuration) {
        let theta = self.base.angle;
        let r1 = 2.0 * self.radius_;
        let r2 = 2.0 * c1.kappa_inv.abs();
        let delta_x = (r1 * r1 + self.base.distance * self.base.distance - r2 * r2) / (2.0 * self.base.distance);
        let delta_y = (r1 * r1 - delta_x * delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.base.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.base.param);
        let q1 = self.tt_tangent_circles(c1, &tgt1);
        let q2 = self.base.tct_tangent_circles(&tgt1, c2);
        let q3 = self.tt_tangent_circles(c1, &tgt2);
        let q4 = self.base.tct_tangent_circles(&tgt2, c2);
        (q1, q2, q3, q4)
    }

    fn ttct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.ttct_tangent_circles(c1, c2);
        let start1 = HcCcCircle::from_configuration(qa, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let start2 = HcCcCircle::from_configuration(qc, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let middle1 = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, true, self.base.param);
        let middle2 = HcCcCircle::from_configuration(qc, !c1.left, c1.forward, true, self.base.param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let q1 = Configuration::new(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
        let l1 = start1.hc_turn_length(&q1) + middle1.hc_turn_length(&qb) + cend.rs_turn_length(&qb);
        let l2 = start2.hc_turn_length(&q1) + middle2.hc_turn_length(&qd) + cend.rs_turn_length(&qd);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(start1), cend: Some(cend), qi1: Some(q1), qi2: Some(qb), ci1: Some(middle1), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(start2), cend: Some(cend), qi1: Some(q1), qi2: Some(qd), ci1: Some(middle2), ..PathSlot::infinite() }
        }
    }

    fn tist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward != c2.forward && self.base.distance >= 2.0 * self.radius_
    }

    fn tist_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration) {
        let dist = center_distance(c1, c2);
        let angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        let alpha = (2.0 * self.radius_ * self.cos_mu_ / dist).asin();
        let dx = self.radius_ * self.sin_mu_;
        let dy = self.radius_ * self.cos_mu_;
        if c1.left && c1.forward {
            let theta = angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else {
            let theta = angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        }
    }

    fn test_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward != c2.forward && self.base.distance >= 2.0 * self.radius_ * self.sin_mu_
    }

    fn test_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration) {
        let dx = self.radius_ * self.sin_mu_;
        let dy = self.radius_ * self.cos_mu_;
        let theta = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        if c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        }
    }

    fn tst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tist_exists(c1, c2) || self.test_exists(c1, c2)
    }

    fn tist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (q2, q3) = self.tist_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let cend = HcCcCircle::from_configuration(q3, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q1 = Configuration::new(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
        let q4 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q2, &q3) + cend.hc_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4), ..PathSlot::infinite() }
    }

    fn test_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (q2, q3) = self.test_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let cend = HcCcCircle::from_configuration(q3, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q1 = Configuration::new(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
        let q4 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q2, &q3) + cend.hc_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4), ..PathSlot::infinite() }
    }

    fn tst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tist_exists(c1, c2) { Some(self.tist_path(c1, c2)) } else if self.test_exists(c1, c2) { Some(self.test_path(c1, c2)) } else { None }
    }

    fn tistct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward == c2.forward
            && self.base.distance >= ((2.0 * self.radius_ * self.sin_mu_ + 2.0 * c1.kappa_inv.abs()).powi(2) + (2.0 * self.radius_ * self.cos_mu_).powi(2)).sqrt()
    }

    fn testct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward == c2.forward && self.base.distance >= 2.0 * (c1.kappa_inv.abs() + self.radius_ * self.sin_mu_)
    }

    fn tistct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.base.angle;
        let delta_y = (4.0 * self.radius_ * self.cos_mu_) / (c2.kappa.abs() * self.base.distance);
        let delta_x = ((2.0 * c2.kappa_inv).powi(2) - delta_y * delta_y).max(0.0).sqrt();
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.base.param);
        let (q2, q3) = self.tist_tangent_circles(c1, &tgt1);
        let q4 = self.base.tct_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let q1 = Configuration::new(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
        let ci = HcCcCircle::from_configuration(q3, !c1.left, c1.forward, true, self.base.param);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q2, &q3) + ci.hc_turn_length(&q4) + cend.rs_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn testct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.base.angle;
        let delta_x = 2.0 * c2.kappa_inv.abs();
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, 0.0);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.base.param);
        let (q2, q3) = self.test_tangent_circles(c1, &tgt1);
        let q4 = self.base.tct_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let q1 = Configuration::new(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
        let ci = HcCcCircle::from_configuration(q3, c1.left, c1.forward, true, self.base.param);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q2, &q3) + ci.hc_turn_length(&q4) + cend.rs_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tstct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tistct_exists(c1, c2) { Some(self.tistct_path(c1, c2)) } else if self.testct_exists(c1, c2) { Some(self.testct_path(c1, c2)) } else { None }
    }

    fn tctist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left
            && c1.forward == c2.forward
            && self.base.distance >= ((2.0 * self.radius_ * self.sin_mu_ + 2.0 * c1.kappa_inv.abs()).powi(2) + (2.0 * self.radius_ * self.cos_mu_).powi(2)).sqrt()
    }

    fn tctest_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward == c2.forward && self.base.distance >= 2.0 * (c1.kappa_inv.abs() + self.radius_ * self.sin_mu_)
    }

    fn tctist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.base.angle;
        let delta_y = (4.0 * self.radius_ * self.cos_mu_) / (c2.kappa.abs() * self.base.distance);
        let delta_x = ((2.0 * c2.kappa_inv).powi(2) - delta_y * delta_y).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward, c2.regular, self.base.param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.tist_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(q3, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q4 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let ci = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.base.param);
        let length = cstart.rs_turn_length(&q1) + ci.hc_turn_length(&q1) + configuration_distance(&q2, &q3) + cend.hc_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tctest_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.base.angle;
        let delta_x = 2.0 * c2.kappa_inv.abs();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, 0.0);
        let tgt1 = HcCcCircle::from_center(x, y, c2.left, !c2.forward, c2.regular, self.base.param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.test_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(q3, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q4 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let ci = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.base.param);
        let length = cstart.rs_turn_length(&q1) + ci.hc_turn_length(&q1) + configuration_distance(&q2, &q3) + cend.hc_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tctst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tctist_exists(c1, c2) { Some(self.tctist_path(c1, c2)) } else if self.tctest_exists(c1, c2) { Some(self.tctest_path(c1, c2)) } else { None }
    }

    fn tctistct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left
            && c1.forward != c2.forward
            && self.base.distance >= ((2.0 * self.radius_).powi(2) + 16.0 * self.radius_ * self.sin_mu_ * c1.kappa_inv.abs() + (4.0 * c1.kappa_inv).powi(2)).sqrt()
    }

    fn tctestct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward != c2.forward && self.base.distance >= 4.0 * c1.kappa_inv.abs() + 2.0 * self.radius_ * self.sin_mu_
    }

    fn tctistct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.base.angle;
        let delta_y = (4.0 * self.radius_ * self.cos_mu_) / (self.base.distance * c1.kappa.abs());
        let delta_x = ((2.0 * c1.kappa_inv).powi(2) - delta_y * delta_y).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.base.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.base.param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.tist_tangent_circles(&tgt1, &tgt2);
        let q4 = self.base.tct_tangent_circles(&tgt2, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let ci1 = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.base.param);
        let ci2 = HcCcCircle::from_configuration(q3, !c2.left, c2.forward, true, self.base.param);
        let length = cstart.rs_turn_length(&q1) + ci1.hc_turn_length(&q1) + configuration_distance(&q2, &q3) + ci2.hc_turn_length(&q4) + cend.rs_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4), ci1: Some(ci1), ci2: Some(ci2) }
    }

    fn tctestct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.base.angle;
        let delta_x = 2.0 * c1.kappa_inv.abs();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, 0.0);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.base.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, 0.0);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.base.param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.test_tangent_circles(&tgt1, &tgt2);
        let q4 = self.base.tct_tangent_circles(&tgt2, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let ci1 = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, true, self.base.param);
        let ci2 = HcCcCircle::from_configuration(q3, !c2.left, c2.forward, true, self.base.param);
        let length = cstart.rs_turn_length(&q1) + ci1.hc_turn_length(&q1) + configuration_distance(&q2, &q3) + ci2.hc_turn_length(&q4) + cend.rs_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4), ci1: Some(ci1), ci2: Some(ci2) }
    }

    fn tctcstct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tctistct_exists(c1, c2) { Some(self.tctistct_path(c1, c2)) } else if self.tctestct_exists(c1, c2) { Some(self.tctestct_path(c1, c2)) } else { None }
    }

    fn ttctt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward == c2.forward && self.base.distance <= 4.0 * self.radius_ + 2.0 * c1.kappa_inv.abs()
    }

    fn ttctt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration, Configuration, Configuration, Configuration, Configuration) {
        let theta = self.base.angle;
        let r1 = 2.0 * c1.kappa_inv.abs();
        let r2 = 2.0 * self.radius_;
        let delta_x = if self.base.distance < 4.0 * self.radius_ - 2.0 * c1.kappa_inv.abs() {
            (self.base.distance + r1) / 2.0
        } else {
            (self.base.distance - r1) / 2.0
        };
        let delta_y = (r2 * r2 - delta_x * delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.base.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward, c2.regular, self.base.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt3 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.base.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y);
        let tgt4 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward, c2.regular, self.base.param);
        let q1 = self.tt_tangent_circles(c1, &tgt1);
        let q2 = self.base.tct_tangent_circles(&tgt1, &tgt2);
        let q3 = self.tt_tangent_circles(&tgt2, c2);
        let q4 = self.tt_tangent_circles(c1, &tgt3);
        let q5 = self.base.tct_tangent_circles(&tgt3, &tgt4);
        let q6 = self.tt_tangent_circles(&tgt4, c2);
        (q1, q2, q3, q4, q5, q6)
    }

    fn ttctt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd, qe, qf) = self.ttctt_tangent_circles(c1, c2);
        let start1 = HcCcCircle::from_configuration(qa, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let middle1 = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, true, self.base.param);
        let middle2 = HcCcCircle::from_configuration(qc, !c2.left, c2.forward, true, self.base.param);
        let end1 = HcCcCircle::from_configuration(qc, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let start2 = HcCcCircle::from_configuration(qd, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let middle3 = HcCcCircle::from_configuration(qd, !c1.left, c1.forward, true, self.base.param);
        let middle4 = HcCcCircle::from_configuration(qf, !c2.left, c2.forward, true, self.base.param);
        let end2 = HcCcCircle::from_configuration(qf, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q1 = Configuration::new(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
        let q3 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let l1 = start1.hc_turn_length(&q1) + middle1.hc_turn_length(&qb) + middle2.hc_turn_length(&qb) + end1.hc_turn_length(&q3);
        let l2 = start2.hc_turn_length(&q1) + middle3.hc_turn_length(&qe) + middle4.hc_turn_length(&qe) + end2.hc_turn_length(&q3);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(start1), cend: Some(end1), qi1: Some(q1), qi2: Some(qb), qi3: Some(q3), ci1: Some(middle1), ci2: Some(middle2), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(start2), cend: Some(end2), qi1: Some(q1), qi2: Some(qe), qi3: Some(q3), ci1: Some(middle3), ci2: Some(middle4), ..PathSlot::infinite() }
        }
    }

    fn tctTct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left
            && c1.forward != c2.forward
            && self.base.distance <= 4.0 * c1.kappa_inv.abs() + 2.0 * self.radius_
            && self.base.distance >= 4.0 * c1.kappa_inv.abs() - 2.0 * self.radius_
    }

    fn tctTct_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration, Configuration, Configuration, Configuration, Configuration) {
        let theta = self.base.angle;
        let r1 = 2.0 * c1.kappa_inv.abs();
        let r2 = self.radius_;
        let delta_x = (r1 * r1 + (self.base.distance / 2.0).powi(2) - r2 * r2) / self.base.distance;
        let delta_y = (r1 * r1 - delta_x * delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.base.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.base.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt3 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.base.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y);
        let tgt4 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.base.param);
        let q1 = self.base.tct_tangent_circles(c1, &tgt1);
        let q2 = self.tt_tangent_circles(&tgt1, &tgt2);
        let q3 = self.base.tct_tangent_circles(&tgt2, c2);
        let q4 = self.base.tct_tangent_circles(c1, &tgt3);
        let q5 = self.tt_tangent_circles(&tgt3, &tgt4);
        let q6 = self.base.tct_tangent_circles(&tgt4, c2);
        (q1, q2, q3, q4, q5, q6)
    }

    fn tctTct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd, qe, qf) = self.tctTct_tangent_circles(c1, c2);
        let middle1 = HcCcCircle::from_configuration(qb, !c1.left, c1.forward, true, self.base.param);
        let middle2 = HcCcCircle::from_configuration(qb, c1.left, !c1.forward, true, self.base.param);
        let middle3 = HcCcCircle::from_configuration(qe, !c1.left, c1.forward, true, self.base.param);
        let middle4 = HcCcCircle::from_configuration(qe, c1.left, !c1.forward, true, self.base.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let l1 = cstart.rs_turn_length(&qa) + middle1.hc_turn_length(&qa) + middle2.hc_turn_length(&qc) + cend.rs_turn_length(&qc);
        let l2 = cstart.rs_turn_length(&qd) + middle3.hc_turn_length(&qd) + middle4.hc_turn_length(&qf) + cend.rs_turn_length(&qf);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(cstart), cend: Some(cend), qi1: Some(qa), qi2: Some(qc), ci1: Some(middle1), ci2: Some(middle2), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(cstart), cend: Some(cend), qi1: Some(qd), qi2: Some(qf), ci1: Some(middle3), ci2: Some(middle4), ..PathSlot::infinite() }
        }
    }

    fn ttt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward != c2.forward && self.base.distance <= 4.0 * self.radius_
    }

    fn ttt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration, Configuration, Configuration) {
        let theta = self.base.angle;
        let r = 2.0 * self.radius_;
        let delta_x = 0.5 * self.base.distance;
        let delta_y = (r * r - delta_x * delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.base.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, self.base.param);
        let q1 = self.tt_tangent_circles(c1, &tgt1);
        let q2 = self.tt_tangent_circles(&tgt1, c2);
        let q3 = self.tt_tangent_circles(c1, &tgt2);
        let q4 = self.tt_tangent_circles(&tgt2, c2);
        (q1, q2, q3, q4)
    }

    fn ttt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.ttt_tangent_circles(c1, c2);
        let start1 = HcCcCircle::from_configuration(qa, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let middle1 = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, CC_REGULAR, self.base.param);
        let end1 = HcCcCircle::from_configuration(qb, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let start2 = HcCcCircle::from_configuration(qc, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let middle2 = HcCcCircle::from_configuration(qc, !c1.left, c1.forward, CC_REGULAR, self.base.param);
        let end2 = HcCcCircle::from_configuration(qd, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let q1 = Configuration::new(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
        let q3 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let l1 = start1.hc_turn_length(&q1) + middle1.cc_turn_length(&qb) + end1.hc_turn_length(&q3);
        let l2 = start2.hc_turn_length(&q1) + middle2.cc_turn_length(&qd) + end2.hc_turn_length(&q3);
        if l1 <= l2 {
            PathSlot { length: l1, cstart: Some(start1), cend: Some(end1), qi1: Some(q1), qi2: Some(qb), qi3: Some(q3), ci1: Some(middle1), ..PathSlot::infinite() }
        } else {
            PathSlot { length: l2, cstart: Some(start2), cend: Some(end2), qi1: Some(q1), qi2: Some(qd), qi3: Some(q3), ci1: Some(middle2), ..PathSlot::infinite() }
        }
    }

    fn tcist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left
            && c1.forward == c2.forward
            && self.base.distance >= ((self.radius_ * self.sin_mu_).powi(2) + (self.radius_ * self.cos_mu_ + c1.kappa_inv.abs()).powi(2)).sqrt()
    }

    fn tcest_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left
            && c1.forward == c2.forward
            && self.base.distance >= ((self.radius_ * self.sin_mu_).powi(2) + (self.radius_ * self.cos_mu_ - c1.kappa_inv.abs()).powi(2)).sqrt()
    }

    fn tcist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let alpha = ((self.radius_ * self.cos_mu_ + c1.kappa_inv.abs()) / self.base.distance).asin();
        let dx1 = 0.0;
        let dy1 = c1.kappa_inv.abs();
        let dx2 = self.radius_ * self.sin_mu_;
        let dy2 = self.radius_ * self.cos_mu_;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        };
        let q3 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(q2, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let length = cstart.rs_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.hc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ..PathSlot::infinite() }
    }

    fn tcest_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let alpha = ((self.radius_ * self.cos_mu_ - c1.kappa_inv.abs()) / self.base.distance).asin();
        let dx1 = 0.0;
        let dy1 = c1.kappa_inv.abs();
        let dx2 = self.radius_ * self.sin_mu_;
        let dy2 = self.radius_ * self.cos_mu_;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, -dy1);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, -dy2);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx1, dy1);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx2, dy2);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        };
        let q3 = Configuration::new(c2.start.x, c2.start.y, c2.start.theta, c2.kappa);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(q2, c2.left, !c2.forward, HC_REGULAR, self.base.param);
        let length = cstart.rs_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.hc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ..PathSlot::infinite() }
    }

    fn tcst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tcist_exists(c1, c2) || self.tcest_exists(c1, c2)
    }

    fn tcst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tcist_exists(c1, c2) { Some(self.tcist_path(c1, c2)) } else if self.tcest_exists(c1, c2) { Some(self.tcest_path(c1, c2)) } else { None }
    }

    fn tisct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left
            && c1.forward == c2.forward
            && self.base.distance >= ((self.radius_ * self.sin_mu_).powi(2) + (self.radius_ * self.cos_mu_ + c1.kappa_inv.abs()).powi(2)).sqrt()
    }

    fn tesct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left
            && c1.forward == c2.forward
            && self.base.distance >= ((self.radius_ * self.sin_mu_).powi(2) + (self.radius_ * self.cos_mu_ - c1.kappa_inv.abs()).powi(2)).sqrt()
    }

    fn tisct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let alpha = ((self.radius_ * self.cos_mu_ + c1.kappa_inv.abs()) / self.base.distance).asin();
        let dx1 = self.radius_ * self.sin_mu_;
        let dy1 = self.radius_ * self.cos_mu_;
        let dx2 = 0.0;
        let dy2 = c1.kappa_inv.abs();
        let (q2, q3) = if c1.left && c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1);
            let q2 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2);
            let q3 = Configuration::new(x, y, theta, c2.kappa);
            (q2, q3)
        } else if c1.left && !c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2);
            let q3 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q2, q3)
        } else if !c1.left && c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1);
            let q2 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2);
            let q3 = Configuration::new(x, y, theta, c2.kappa);
            (q2, q3)
        } else {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2);
            let q3 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q2, q3)
        };
        let q1 = Configuration::new(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q2, &q3) + cend.rs_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ..PathSlot::infinite() }
    }

    fn tesct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let alpha = ((self.radius_ * self.cos_mu_ - c1.kappa_inv.abs()) / self.base.distance).asin();
        let dx1 = self.radius_ * self.sin_mu_;
        let dy1 = self.radius_ * self.cos_mu_;
        let dx2 = 0.0;
        let dy2 = c1.kappa_inv.abs();
        let (q2, q3) = if c1.left && c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1);
            let q2 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2);
            let q3 = Configuration::new(x, y, theta, c2.kappa);
            (q2, q3)
        } else if c1.left && !c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2);
            let q3 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q2, q3)
        } else if !c1.left && c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx1, dy1);
            let q2 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx2, dy2);
            let q3 = Configuration::new(x, y, theta, c2.kappa);
            (q2, q3)
        } else {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx1, -dy1);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx2, -dy2);
            let q3 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q2, q3)
        };
        let q1 = Configuration::new(c1.start.x, c1.start.y, c1.start.theta, c1.kappa);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, HC_REGULAR, self.base.param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let length = cstart.hc_turn_length(&q1) + configuration_distance(&q2, &q3) + cend.rs_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ..PathSlot::infinite() }
    }

    fn tsct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tisct_exists(c1, c2) || self.tesct_exists(c1, c2)
    }

    fn tsct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tisct_exists(c1, c2) { Some(self.tisct_path(c1, c2)) } else if self.tesct_exists(c1, c2) { Some(self.tesct_path(c1, c2)) } else { None }
    }

    fn tcisct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward != c2.forward && self.base.distance > 2.0 * c1.kappa_inv.abs()
    }

    fn tcesct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward != c2.forward && self.base.distance >= get_epsilon()
    }

    fn tcisct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let alpha = (2.0 / (c1.kappa.abs() * self.base.distance)).asin();
        let dx = 0.0;
        let dy = c1.kappa_inv.abs();
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx, dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.base.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx, dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else {
            let theta = self.base.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx, -dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let length = cstart.rs_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.rs_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcesct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.base.angle;
        let dx = 0.0;
        let dy = c1.kappa_inv.abs();
        let (q1, q2) = if c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx, dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx, -dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, c2.kappa);
            (q1, q2)
        } else {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, dy);
            let q1 = Configuration::new(x, y, theta, c1.kappa);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, dx, dy);
            let q2 = Configuration::new(x, y, theta, c2.kappa);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, self.base.rs_param);
        let cend = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, self.base.rs_param);
        let length = cstart.rs_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.rs_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcsct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tcisct_exists(c1, c2) || self.tcesct_exists(c1, c2)
    }

    fn tcsct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tcisct_exists(c1, c2) { Some(self.tcisct_path(c1, c2)) } else if self.tcesct_exists(c1, c2) { Some(self.tcesct_path(c1, c2)) } else { None }
    }
}

fn hcpmpm_circles_rs_path(
    c1: &HcCcCircle,
    c2: &HcCcCircle,
    param: &HcCcCircleParam,
    rs_param: &HcCcCircleParam,
) -> HcCcRsPath {
    use HcCcRsPathType::*;

    let mut h = HcpmpmRsHelper::new(param, rs_param);
    h.base.distance = center_distance(c1, c2);
    h.base.angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);

    let mut slots: Vec<PathSlot> = (0..NB_HC_CC_RS_PATHS).map(|_| PathSlot::infinite()).collect();
    let mut skip = false;

    if configuration_equal(&c1.start, &c2.start) {
        slots[E as usize].length = 0.0;
        skip = true;
    }

    if !skip && configuration_on_hc_cc_circle(c1, &c2.start) {
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, true, rs_param);
        slots[T as usize].length = cstart.rs_turn_length(&c2.start);
        slots[T as usize].cstart = Some(cstart);
        skip = true;
    }

    if !skip {
        if h.tt_exists(c1, c2) { slots[TT as usize] = h.tt_path(c1, c2); }
        if h.base.tct_exists(c1, c2) { slots[TcT as usize] = h.tct_path(c1, c2); }
        if h.base.tctct_exists(c1, c2) { slots[TcTcT as usize] = h.tctct_path(c1, c2); }
        if h.tctt_exists(c1, c2) { slots[TcTT as usize] = h.tctt_path(c1, c2); }
        if h.ttct_exists(c1, c2) { slots[TTcT as usize] = h.ttct_path(c1, c2); }
        if h.tst_exists(c1, c2) { if let Some(s) = h.tst_path(c1, c2) { slots[TST as usize] = s; } }
        if h.tistct_exists(c1, c2) || h.testct_exists(c1, c2) {
            if let Some(s) = h.tstct_path(c1, c2) { slots[TSTcT as usize] = s; }
        }
        if h.tctist_exists(c1, c2) || h.tctest_exists(c1, c2) {
            if let Some(s) = h.tctst_path(c1, c2) { slots[TcTST as usize] = s; }
        }
        if h.tctistct_exists(c1, c2) || h.tctestct_exists(c1, c2) {
            if let Some(s) = h.tctcstct_path(c1, c2) { slots[TcTSTcT as usize] = s; }
        }
        if h.ttctt_exists(c1, c2) { slots[TTcTT as usize] = h.ttctt_path(c1, c2); }
        if h.tctTct_exists(c1, c2) { slots[TcTTcT as usize] = h.tctTct_path(c1, c2); }
        if h.ttt_exists(c1, c2) { slots[TTT as usize] = h.ttt_path(c1, c2); }
        if h.tcst_exists(c1, c2) { if let Some(s) = h.tcst_path(c1, c2) { slots[TcST as usize] = s; } }
        if h.tsct_exists(c1, c2) { if let Some(s) = h.tsct_path(c1, c2) { slots[TScT as usize] = s; } }
        if h.tcsct_exists(c1, c2) { if let Some(s) = h.tcsct_path(c1, c2) { slots[TcScT as usize] = s; } }
    }

    let best = (0..NB_HC_CC_RS_PATHS)
        .min_by(|&a, &b| slots[a].length.partial_cmp(&slots[b].length).unwrap())
        .unwrap();

    let s = &slots[best];
    HcCcRsPath::new(
        c1.start,
        c2.start,
        path_type_from_usize(best),
        param.kappa,
        param.sigma,
        s.qi1,
        s.qi2,
        s.qi3,
        s.qi4,
        s.cstart.clone().map(Box::new),
        s.cend.clone().map(Box::new),
        s.ci1.clone().map(Box::new),
        s.ci2.clone().map(Box::new),
        s.length,
    )
}

pub struct HcpmpmRsStateSpace {
    params_: HcCcStateSpaceParams,
    rs_circle_param_: HcCcCircleParam,
    discretization_: f64,
}

impl HcpmpmRsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64) -> Self {
        let params = HcCcStateSpaceParams::new(kappa, sigma);
        let mut rs_param = HcCcCircleParam::default();
        rs_param.set_param(kappa, f64::MAX, 1.0 / kappa, 0.0, 0.0, 1.0, 0.0);
        Self {
            params_: params,
            rs_circle_param_: rs_param,
            discretization_: discretization,
        }
    }

    fn hcpmpm_reeds_shepp(&self, state1: &State, state2: &State) -> HcCcRsPath {
        let kappa = self.params_.kappa_;
        let rsp = &self.rs_circle_param_;
        let start1 = Configuration::new(state1.x, state1.y, state1.theta, kappa);
        let start2 = Configuration::new(state1.x, state1.y, state1.theta, -kappa);
        let end1 = Configuration::new(state2.x, state2.y, state2.theta, kappa);
        let end2 = Configuration::new(state2.x, state2.y, state2.theta, -kappa);
        let start_circles = [
            HcCcCircle::from_configuration(start1, true, true, true, rsp),
            HcCcCircle::from_configuration(start2, false, true, true, rsp),
            HcCcCircle::from_configuration(start1, true, false, true, rsp),
            HcCcCircle::from_configuration(start2, false, false, true, rsp),
        ];
        let end_circles = [
            HcCcCircle::from_configuration(end1, true, true, true, rsp),
            HcCcCircle::from_configuration(end2, false, true, true, rsp),
            HcCcCircle::from_configuration(end1, true, false, true, rsp),
            HcCcCircle::from_configuration(end2, false, false, true, rsp),
        ];

        let mut best: Option<HcCcRsPath> = None;
        for (i, sc) in start_circles.iter().enumerate() {
            if (i == 0 || i == 2) && state1.kappa < 0.0 { continue; }
            if (i == 1 || i == 3) && state1.kappa > 0.0 { continue; }
            for (j, ec) in end_circles.iter().enumerate() {
                if (j == 0 || j == 2) && state2.kappa < 0.0 { continue; }
                if (j == 1 || j == 3) && state2.kappa > 0.0 { continue; }
                let path = hcpmpm_circles_rs_path(sc, ec, &self.params_.hc_cc_circle_param_, &self.rs_circle_param_);
                if best.as_ref().map_or(true, |b| path.length < b.length) {
                    best = Some(path);
                }
            }
        }
        best.unwrap()
    }
}

impl StateSpace for HcpmpmRsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let path = self.hcpmpm_reeds_shepp(s1, s2);
        let mut controls = Vec::new();

        match path.path_type {
            HcCcRsPathType::E => empty_controls(&mut controls),
            HcCcRsPathType::T => {
                rs_turn_controls(path.cstart.as_ref().unwrap(), &path.end, true, &mut controls);
            }
            HcCcRsPathType::TT => {
                hc_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                hc_turn_controls(path.cend.as_ref().unwrap(), path.qi3.as_ref().unwrap(), true, &mut controls);
            }
            HcCcRsPathType::TcT => {
                rs_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), true, &mut controls);
                rs_turn_controls(path.cend.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
            }
            HcCcRsPathType::TcTcT => {
                rs_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), true, &mut controls);
                rs_turn_controls(path.ci1.as_ref().unwrap(), path.qi2.as_ref().unwrap(), true, &mut controls);
                rs_turn_controls(path.cend.as_ref().unwrap(), path.qi2.as_ref().unwrap(), false, &mut controls);
            }
            HcCcRsPathType::TcTT => {
                rs_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), true, &mut controls);
                hc_turn_controls(path.ci1.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                hc_turn_controls(path.cend.as_ref().unwrap(), path.qi2.as_ref().unwrap(), true, &mut controls);
            }
            HcCcRsPathType::TTcT => {
                hc_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                hc_turn_controls(path.ci1.as_ref().unwrap(), path.qi2.as_ref().unwrap(), true, &mut controls);
                rs_turn_controls(path.cend.as_ref().unwrap(), path.qi2.as_ref().unwrap(), false, &mut controls);
            }
            HcCcRsPathType::TST => {
                hc_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                straight_controls(path.qi2.as_ref().unwrap(), path.qi3.as_ref().unwrap(), &mut controls);
                hc_turn_controls(path.cend.as_ref().unwrap(), path.qi4.as_ref().unwrap(), true, &mut controls);
            }
            HcCcRsPathType::TSTcT => {
                hc_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                straight_controls(path.qi2.as_ref().unwrap(), path.qi3.as_ref().unwrap(), &mut controls);
                hc_turn_controls(path.ci1.as_ref().unwrap(), path.qi4.as_ref().unwrap(), true, &mut controls);
                rs_turn_controls(path.cend.as_ref().unwrap(), path.qi4.as_ref().unwrap(), false, &mut controls);
            }
            HcCcRsPathType::TcTST => {
                rs_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), true, &mut controls);
                hc_turn_controls(path.ci1.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                straight_controls(path.qi2.as_ref().unwrap(), path.qi3.as_ref().unwrap(), &mut controls);
                hc_turn_controls(path.cend.as_ref().unwrap(), path.qi4.as_ref().unwrap(), true, &mut controls);
            }
            HcCcRsPathType::TcTSTcT => {
                rs_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), true, &mut controls);
                hc_turn_controls(path.ci1.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                straight_controls(path.qi2.as_ref().unwrap(), path.qi3.as_ref().unwrap(), &mut controls);
                hc_turn_controls(path.ci2.as_ref().unwrap(), path.qi4.as_ref().unwrap(), true, &mut controls);
                rs_turn_controls(path.cend.as_ref().unwrap(), path.qi4.as_ref().unwrap(), false, &mut controls);
            }
            HcCcRsPathType::TTcTT => {
                hc_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                hc_turn_controls(path.ci1.as_ref().unwrap(), path.qi2.as_ref().unwrap(), true, &mut controls);
                hc_turn_controls(path.ci2.as_ref().unwrap(), path.qi2.as_ref().unwrap(), false, &mut controls);
                hc_turn_controls(path.cend.as_ref().unwrap(), path.qi3.as_ref().unwrap(), true, &mut controls);
            }
            HcCcRsPathType::TcTTcT => {
                rs_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), true, &mut controls);
                hc_turn_controls(path.ci1.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                hc_turn_controls(path.ci2.as_ref().unwrap(), path.qi2.as_ref().unwrap(), true, &mut controls);
                rs_turn_controls(path.cend.as_ref().unwrap(), path.qi2.as_ref().unwrap(), false, &mut controls);
            }
            HcCcRsPathType::TTT => {
                hc_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                cc_turn_controls(path.ci1.as_ref().unwrap(), path.qi2.as_ref().unwrap(), true, &mut controls);
                hc_turn_controls(path.cend.as_ref().unwrap(), path.qi3.as_ref().unwrap(), true, &mut controls);
            }
            HcCcRsPathType::TcST => {
                rs_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), true, &mut controls);
                straight_controls(path.qi1.as_ref().unwrap(), path.qi2.as_ref().unwrap(), &mut controls);
                hc_turn_controls(path.cend.as_ref().unwrap(), path.qi3.as_ref().unwrap(), true, &mut controls);
            }
            HcCcRsPathType::TScT => {
                hc_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), false, &mut controls);
                straight_controls(path.qi2.as_ref().unwrap(), path.qi3.as_ref().unwrap(), &mut controls);
                rs_turn_controls(path.cend.as_ref().unwrap(), path.qi3.as_ref().unwrap(), false, &mut controls);
            }
            HcCcRsPathType::TcScT => {
                rs_turn_controls(path.cstart.as_ref().unwrap(), path.qi1.as_ref().unwrap(), true, &mut controls);
                straight_controls(path.qi1.as_ref().unwrap(), path.qi2.as_ref().unwrap(), &mut controls);
                rs_turn_controls(path.cend.as_ref().unwrap(), path.qi2.as_ref().unwrap(), false, &mut controls);
            }
            HcCcRsPathType::S => {
                straight_controls(&path.start, &path.end, &mut controls);
            }
        }

        controls
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        vec![self.get_controls(s1, s2)]
    }

    fn discretization(&self) -> f64 {
        self.discretization_
    }
}

pub struct HcRsStateSpace {
    kappa_: f64,
    sigma_: f64,
    discretization_: f64,
    hc00_: Hc00RsStateSpace,
    hc0pm_: Hc0pmRsStateSpace,
    hcpm0_: Hcpm0RsStateSpace,
    hcpmpm_: HcpmpmRsStateSpace,
}

impl HcRsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64) -> Self {
        Self {
            kappa_: kappa,
            sigma_: sigma,
            discretization_: discretization,
            hc00_: Hc00RsStateSpace::new(kappa, sigma, discretization),
            hc0pm_: Hc0pmRsStateSpace::new(kappa, sigma, discretization),
            hcpm0_: Hcpm0RsStateSpace::new(kappa, sigma, discretization),
            hcpmpm_: HcpmpmRsStateSpace::new(kappa, sigma, discretization),
        }
    }

    fn predict_state(&self, state: &State) -> Vec<(State, Control)> {
        let eps = get_epsilon();
        if state.kappa.abs() < eps || (self.kappa_ - state.kappa.abs()).abs() < eps {
            return vec![(
                *state,
                Control {
                    delta_s: 0.0,
                    kappa: state.kappa,
                    sigma: 0.0,
                },
            )];
        }

        let mut results = Vec::with_capacity(4);
        let sgn_kappa = sgn(state.kappa);

        let mut c1 = Control {
            delta_s: (self.kappa_ - sgn_kappa * state.kappa) / self.sigma_,
            kappa: state.kappa,
            sigma: sgn_kappa * self.sigma_,
        };
        if state.kappa.abs() > self.kappa_ {
            c1.sigma = -sgn_kappa * self.sigma_;
        }

        let controls = [
            c1,
            Control {
                delta_s: -c1.delta_s,
                kappa: state.kappa,
                sigma: c1.sigma,
            },
            Control {
                delta_s: sgn_kappa * state.kappa / self.sigma_,
                kappa: state.kappa,
                sigma: -sgn_kappa * self.sigma_,
            },
            Control {
                delta_s: -(sgn_kappa * state.kappa / self.sigma_),
                kappa: state.kappa,
                sigma: -sgn_kappa * self.sigma_,
            },
        ];

        for ctrl in controls {
            let d = sgn(ctrl.delta_s);
            let abs_ds = ctrl.delta_s.abs();
            let (x, y, theta, kappa) = end_of_clothoid(
                state.x,
                state.y,
                state.theta,
                state.kappa,
                ctrl.sigma,
                d,
                abs_ds,
            );
            results.push((
                State {
                    x,
                    y,
                    theta,
                    kappa,
                    ..State::default()
                },
                ctrl,
            ));
        }

        results
    }

    fn candidate_controls(&self, state1: &State, state2: &State) -> Vec<Vec<Control>> {
        let eps = get_epsilon();
        let start_predictions = self.predict_state(state1);
        let end_predictions = self.predict_state(state2);
        let mut candidates = Vec::new();

        for (start_state, start_control) in &start_predictions {
            for (end_state, end_control) in &end_predictions {
                let mut controls = Vec::new();
                if state_equal(&state_to_configuration(start_state), &state_to_configuration(end_state)) {
                    controls.push(subtract_control(start_control, end_control));
                } else {
                    if start_state.kappa.abs() < eps {
                        if end_state.kappa.abs() < eps {
                            controls = self.hc00_.get_controls(start_state, end_state);
                        } else {
                            controls = self.hc0pm_.get_controls(start_state, end_state);
                        }
                    } else if end_state.kappa.abs() < eps {
                        controls = self.hcpm0_.get_controls(start_state, end_state);
                    } else {
                        controls = self.hcpmpm_.get_controls(start_state, end_state);
                    }

                    if start_control.delta_s.abs() > eps {
                        controls.insert(0, *start_control);
                    }
                    if end_control.delta_s.abs() > eps {
                        let mut reversed_end = *end_control;
                        reverse_control(&mut reversed_end);
                        controls.push(reversed_end);
                    }
                }
                candidates.push(controls);
            }
        }

        candidates
    }
}

impl StateSpace for HcRsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let mut candidates = self.candidate_controls(s1, s2);
        candidates.sort_by(|left, right| {
            controls_length(left)
                .partial_cmp(&controls_length(right))
                .unwrap()
        });
        candidates.into_iter().next().unwrap_or_default()
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        let mut candidates = self.candidate_controls(s1, s2);
        candidates.sort_by(|left, right| {
            controls_length(left)
                .partial_cmp(&controls_length(right))
                .unwrap()
        });
        candidates
    }

    fn discretization(&self) -> f64 {
        self.discretization_
    }
}
