use std::f64::consts::PI;
use crate::utilities::{get_epsilon, global_frame_change, HALF_PI, sgn, end_of_clothoid, point_distance};
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
