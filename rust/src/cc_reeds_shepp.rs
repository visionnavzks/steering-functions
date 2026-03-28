use std::f64::consts::PI;
use crate::utilities::{get_epsilon, global_frame_change, HALF_PI};
use crate::configuration::{Configuration, configuration_distance, configuration_equal, configuration_aligned};
use crate::hc_cc_circle::{HcCcCircle, HcCcCircleParam, center_distance, configuration_on_hc_cc_circle};
use crate::paths::{
    HcCcRsPath, HcCcRsPathType, NB_HC_CC_RS_PATHS,
    empty_controls, straight_controls, cc_turn_controls,
};
use crate::hc_cc_state_space::HcCcStateSpaceParams;
use crate::state::{State, Control};
use crate::base_state_space::StateSpace;

const CC_REGULAR: bool = false;

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
// Geometry helper
// ---------------------------------------------------------------------------

struct Cc00RsHelper<'a> {
    param:    &'a HcCcCircleParam,
    distance: f64,
    angle:    f64,
}

impl<'a> Cc00RsHelper<'a> {
    fn new(param: &'a HcCcCircleParam) -> Self {
        Self { param, distance: 0.0, angle: 0.0 }
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

    // ---- TcT ---------------------------------------------------------------
    fn tct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left
            && c1.forward == c2.forward
            && (self.distance - 2.0 * c1.radius * c1.cos_mu).abs() < get_epsilon()
    }

    fn tct_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Configuration {
        let dist    = center_distance(c1, c2);
        let delta_x = 0.5 * dist;
        let delta_y = (c1.radius * c1.radius - delta_x * delta_x).max(0.0).sqrt();
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
        Configuration::new(x, y, theta, 0.0)
    }

    fn tct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let q = self.tct_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.cc_turn_length(&q) + cend.cc_turn_length(&q);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend), qi1: Some(q), ..PathSlot::infinite() }
    }

    // ---- TcTcT -------------------------------------------------------------
    fn tctct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left
            && c1.forward != c2.forward
            && self.distance <= 4.0 * c1.radius * c1.cos_mu
    }

    fn tctct_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration)
    {
        let theta   = self.angle;
        let r       = 2.0 * c1.radius * c1.cos_mu;
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
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, !c1.forward, c1.regular, self.param);
        let mid2   = HcCcCircle::from_configuration(qc, !c1.left, !c1.forward, c1.regular, self.param);
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

    // ---- TcTT --------------------------------------------------------------
    fn tctt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward == c2.forward
            && self.distance <= 2.0 * c1.radius * (1.0 + c1.cos_mu)
            && self.distance >= 2.0 * c1.radius * (1.0 - c1.cos_mu)
    }

    fn tctt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration)
    {
        let theta   = self.angle;
        let r1      = 2.0 * c1.radius * c1.cos_mu;
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
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, !c1.forward, c1.regular, self.param);
        let mid2   = HcCcCircle::from_configuration(qc, !c1.left, !c1.forward, c1.regular, self.param);
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

    // ---- TTcT --------------------------------------------------------------
    fn ttct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward == c2.forward
            && self.distance <= 2.0 * c1.radius * (1.0 + c1.cos_mu)
            && self.distance >= 2.0 * c1.radius * (1.0 - c1.cos_mu)
    }

    fn ttct_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration)
    {
        let theta   = self.angle;
        let r1      = 2.0 * c1.radius;
        let r2      = 2.0 * c1.radius * c1.cos_mu;
        let delta_x = (r1*r1 + self.distance*self.distance - r2*r2) / (2.0 * self.distance);
        let delta_y = (r1*r1 - delta_x*delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left,  c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left,  c1.forward, c1.regular, self.param);
        let q1 = self.tt_tangent_circles(c1,   &tgt1);
        let q2 = self.tct_tangent_circles(&tgt1, c2);
        let q3 = self.tt_tangent_circles(c1,   &tgt2);
        let q4 = self.tct_tangent_circles(&tgt2, c2);
        (q1, q2, q3, q4)
    }

    fn ttct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let (qa, qb, qc, qd) = self.ttct_tangent_circles(c1, c2);
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left,  c1.forward, c1.regular, self.param);
        let mid2   = HcCcCircle::from_configuration(qc, !c1.left,  c1.forward, c1.regular, self.param);
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

    // ---- TST ---------------------------------------------------------------
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

    // ---- TSTcT -------------------------------------------------------------
    fn tistct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward == c2.forward
            && self.distance >= 2.0 * c1.radius
                * (1.0 + 2.0*c1.sin_mu*c1.cos_mu + c1.cos_mu*c1.cos_mu).sqrt()
    }
    fn testct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward == c2.forward
            && self.distance >= 2.0 * c1.radius * (c1.cos_mu + c1.sin_mu)
    }
    fn tstct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tistct_exists(c1, c2) || self.testct_exists(c1, c2)
    }

    fn tistct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let r       = c2.radius * c2.cos_mu;
        let delta_y = (2.0*r)*(2.0*r) / self.distance;
        let delta_x = 2.0*r * (1.0 - delta_y/self.distance).max(0.0).sqrt();
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.param);
        let (q1, q2) = self.tist_tangent_circles(c1, &tgt1);
        let q3 = self.tct_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let ci = HcCcCircle::from_configuration(q2, !c1.left, c1.forward, c1.regular, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2)
            + ci.cc_turn_length(&q3) + cend.cc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn testct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let delta_x = 2.0 * c2.radius * c2.cos_mu;
        let (x, y)  = global_frame_change(c2.xc, c2.yc, theta, -delta_x, 0.0);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, c2.forward, c2.regular, self.param);
        let (q1, q2) = self.test_tangent_circles(c1, &tgt1);
        let q3 = self.tct_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let ci = HcCcCircle::from_configuration(q2, c1.left, c1.forward, c1.regular, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2)
            + ci.cc_turn_length(&q3) + cend.cc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tstct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tistct_exists(c1, c2) { Some(self.tistct_path(c1, c2)) }
        else if self.testct_exists(c1, c2) { Some(self.testct_path(c1, c2)) }
        else { None }
    }

    // ---- TcTST -------------------------------------------------------------
    fn tctist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward == c2.forward
            && self.distance >= 2.0 * c1.radius
                * (1.0 + 2.0*c1.sin_mu*c1.cos_mu + c1.cos_mu*c1.cos_mu).sqrt()
    }
    fn tctest_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward == c2.forward
            && self.distance >= 2.0 * c1.radius * (c1.cos_mu + c1.sin_mu)
    }
    fn tctst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tctist_exists(c1, c2) || self.tctest_exists(c1, c2)
    }

    fn tctist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let r       = c1.radius * c1.cos_mu;
        let delta_y = (2.0*r)*(2.0*r) / self.distance;
        let delta_x = 2.0*r * (1.0 - delta_y/self.distance).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward, c2.regular, self.param);
        let q1 = self.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.tist_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let ci = HcCcCircle::from_configuration(q1, !c1.left, !c1.forward, c1.regular, self.param);
        let length = cstart.cc_turn_length(&q1) + ci.cc_turn_length(&q2)
            + configuration_distance(&q2, &q3) + cend.cc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tctest_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let delta_x = 2.0 * c2.radius * c2.cos_mu;
        let (x, y)  = global_frame_change(c1.xc, c1.yc, theta, delta_x, 0.0);
        let tgt1 = HcCcCircle::from_center(x, y, c2.left, !c2.forward, c2.regular, self.param);
        let q1 = self.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.test_tangent_circles(&tgt1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let ci = HcCcCircle::from_configuration(q1, !c1.left, !c1.forward, c1.regular, self.param);
        let length = cstart.cc_turn_length(&q1) + ci.cc_turn_length(&q2)
            + configuration_distance(&q2, &q3) + cend.cc_turn_length(&q3);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), ci1: Some(ci), ..PathSlot::infinite() }
    }

    fn tctst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tctist_exists(c1, c2) { Some(self.tctist_path(c1, c2)) }
        else if self.tctest_exists(c1, c2) { Some(self.tctest_path(c1, c2)) }
        else { None }
    }

    // ---- TcTSTcT -----------------------------------------------------------
    fn tctistct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward != c2.forward
            && self.distance >= 2.0 * c1.radius
                * (1.0 + 4.0*c1.cos_mu*c1.sin_mu + 4.0*c1.cos_mu*c1.cos_mu).sqrt()
    }
    fn tctestct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward != c2.forward
            && self.distance >= 2.0 * c1.radius * (2.0*c1.cos_mu + c1.sin_mu)
    }
    fn tctcstct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tctistct_exists(c1, c2) || self.tctestct_exists(c1, c2)
    }

    fn tctistct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let r       = c1.radius * c1.cos_mu;
        let delta_y = (2.0*r)*(2.0*r) / self.distance;
        let delta_x = 2.0*r * (1.0 - delta_y/self.distance).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left,  c2.forward, c2.regular, self.param);
        let q1 = self.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.tist_tangent_circles(&tgt1, &tgt2);
        let q4 = self.tct_tangent_circles(&tgt2, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let ci1 = HcCcCircle::from_configuration(q1, !c1.left, !c1.forward, c1.regular, self.param);
        let ci2 = HcCcCircle::from_configuration(q3, !c2.left,  c2.forward, c2.regular, self.param);
        let length = cstart.cc_turn_length(&q1) + ci1.cc_turn_length(&q2)
            + configuration_distance(&q2, &q3) + ci2.cc_turn_length(&q4)
            + cend.cc_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4),
                   ci1: Some(ci1), ci2: Some(ci2) }
    }

    fn tctestct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta   = self.angle;
        let delta_x = 2.0 * c1.radius * c1.cos_mu;
        let (x, y)  = global_frame_change(c1.xc, c1.yc, theta,  delta_x, 0.0);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, !c1.forward, c1.regular, self.param);
        let (x, y)  = global_frame_change(c2.xc, c2.yc, theta, -delta_x, 0.0);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left,  c2.forward, c2.regular, self.param);
        let q1 = self.tct_tangent_circles(c1, &tgt1);
        let (q2, q3) = self.test_tangent_circles(&tgt1, &tgt2);
        let q4 = self.tct_tangent_circles(&tgt2, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let ci1 = HcCcCircle::from_configuration(q1, !c1.left, !c1.forward, c1.regular, self.param);
        let ci2 = HcCcCircle::from_configuration(q3, !c2.left,  c2.forward, c2.regular, self.param);
        let length = cstart.cc_turn_length(&q1) + ci1.cc_turn_length(&q2)
            + configuration_distance(&q2, &q3) + ci2.cc_turn_length(&q4)
            + cend.cc_turn_length(&q4);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), qi3: Some(q3), qi4: Some(q4),
                   ci1: Some(ci1), ci2: Some(ci2) }
    }

    fn tctcstct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tctistct_exists(c1, c2) { Some(self.tctistct_path(c1, c2)) }
        else if self.tctestct_exists(c1, c2) { Some(self.tctestct_path(c1, c2)) }
        else { None }
    }

    // ---- TTcTT -------------------------------------------------------------
    fn ttctt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward == c2.forward
            && self.distance <= 2.0 * c1.radius * (c1.cos_mu + 2.0)
    }

    fn ttctt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration, Configuration, Configuration)
    {
        let theta   = self.angle;
        let r1      = 2.0 * c1.radius * c1.cos_mu;
        let r2      = 2.0 * c1.radius;
        let delta_x = if self.distance < 2.0 * c1.radius * (2.0 - c1.cos_mu) {
            (self.distance + r1) / 2.0
        } else {
            (self.distance - r1) / 2.0
        };
        let delta_y = (r2*r2 - delta_x*delta_x).max(0.0).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left,  c1.forward,  c1.regular, self.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x,  delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward,  c2.regular, self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x, -delta_y);
        let tgt3 = HcCcCircle::from_center(x, y, !c1.left,  c1.forward,  c1.regular, self.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y);
        let tgt4 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward,  c2.regular, self.param);
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
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left,  c1.forward,  c1.regular, self.param);
        let mid2   = HcCcCircle::from_configuration(qb, !c2.left, !c2.forward,  c2.regular, self.param);
        let mid3   = HcCcCircle::from_configuration(qd, !c1.left,  c1.forward,  c1.regular, self.param);
        let mid4   = HcCcCircle::from_configuration(qe, !c2.left, !c2.forward,  c2.regular, self.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let l1 = cstart.cc_turn_length(&qa) + mid1.cc_turn_length(&qb)
            + mid2.cc_turn_length(&qc) + cend.cc_turn_length(&qc);
        let l2 = cstart.cc_turn_length(&qd) + mid3.cc_turn_length(&qe)
            + mid4.cc_turn_length(&qf) + cend.cc_turn_length(&qf);
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

    // ---- TcTTcT ------------------------------------------------------------
    fn tctTct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward != c2.forward
            && self.distance <= 2.0 * c1.radius * (2.0*c1.cos_mu + 1.0)
            && self.distance >= 2.0 * c1.radius * (2.0*c1.cos_mu - 1.0)
    }

    fn tctTct_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (Configuration, Configuration, Configuration, Configuration, Configuration, Configuration)
    {
        let theta   = self.angle;
        let r1      = 2.0 * c1.radius * c1.cos_mu;
        let r2      = c1.radius;
        let half_d  = self.distance / 2.0;
        let delta_x = (r1*r1 + half_d*half_d - r2*r2) / self.distance;
        let delta_y = (r1*r1 - delta_x*delta_x).max(0.0).sqrt();
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
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, !c1.forward, c1.regular, self.param);
        let mid2   = HcCcCircle::from_configuration(qb,  c1.left, !c1.forward, c1.regular, self.param);
        let mid3   = HcCcCircle::from_configuration(qd, !c1.left, !c1.forward, c1.regular, self.param);
        let mid4   = HcCcCircle::from_configuration(qe,  c1.left, !c1.forward, c1.regular, self.param);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let l1 = cstart.cc_turn_length(&qa) + mid1.cc_turn_length(&qb)
            + mid2.cc_turn_length(&qc) + cend.cc_turn_length(&qc);
        let l2 = cstart.cc_turn_length(&qd) + mid3.cc_turn_length(&qe)
            + mid4.cc_turn_length(&qf) + cend.cc_turn_length(&qf);
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

    // ---- TTT ---------------------------------------------------------------
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
        let mid1   = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, c1.regular, self.param);
        let mid2   = HcCcCircle::from_configuration(qc, !c1.left, c1.forward, c1.regular, self.param);
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

    // ---- TcST --------------------------------------------------------------
    fn tcist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward == c2.forward
            && self.distance >= 2.0 * c1.radius * c1.cos_mu
    }
    fn tcest_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward == c2.forward && self.distance >= get_epsilon()
    }
    fn tcst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tcist_exists(c1, c2) || self.tcest_exists(c1, c2)
    }

    fn tcist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let alpha = (2.0 * c1.radius * c1.cos_mu / self.distance).asin();
        let dx = c1.radius * c1.sin_mu;
        let dy = c1.radius * c1.cos_mu;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx,  dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.cc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcest_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.angle;
        let dx = c1.radius * c1.sin_mu;
        let dy = c1.radius * c1.cos_mu;
        let (q1, q2) = if c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx,  dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.cc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tcist_exists(c1, c2) { Some(self.tcist_path(c1, c2)) }
        else if self.tcest_exists(c1, c2) { Some(self.tcest_path(c1, c2)) }
        else { None }
    }

    // ---- TScT --------------------------------------------------------------
    fn tisce_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward == c2.forward
            && self.distance >= 2.0 * c1.radius * c1.cos_mu
    }
    fn tesc_t_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward == c2.forward && self.distance >= get_epsilon()
    }
    fn tsct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tisce_exists(c1, c2) || self.tesc_t_exists(c1, c2)
    }

    fn tiscT_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let alpha = (2.0 * c1.radius * c1.cos_mu / self.distance).asin();
        let dx = c1.radius * c1.sin_mu;
        let dy = c1.radius * c1.cos_mu;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx,  dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.cc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tescT_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.angle;
        let dx = c1.radius * c1.sin_mu;
        let dy = c1.radius * c1.cos_mu;
        let (q1, q2) = if c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx,  dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.cc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tsct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tisce_exists(c1, c2) { Some(self.tiscT_path(c1, c2)) }
        else if self.tesc_t_exists(c1, c2) { Some(self.tescT_path(c1, c2)) }
        else { None }
    }

    // ---- TcScT -------------------------------------------------------------
    fn tcisce_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left != c2.left && c1.forward != c2.forward
            && self.distance >= 2.0 * c1.radius * c1.cos_mu
    }
    fn tcesce_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        c1.left == c2.left && c1.forward != c2.forward && self.distance >= get_epsilon()
    }
    fn tcsct_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tcisce_exists(c1, c2) || self.tcesce_exists(c1, c2)
    }

    fn tciscT_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let alpha = (2.0 * c1.radius * c1.cos_mu / self.distance).asin();
        let dx = c1.radius * c1.sin_mu;
        let dy = c1.radius * c1.cos_mu;
        let (q1, q2) = if c1.left && c1.forward {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let theta = self.angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else {
            let theta = self.angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.cc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcescT_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> PathSlot {
        let theta = self.angle;
        let dx = c1.radius * c1.sin_mu;
        let dy = c1.radius * c1.cos_mu;
        let (q1, q2) = if c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else if c1.left && !c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        } else if !c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx, -dy);
            let q2 = Configuration::new(x, y, theta + PI, 0.0);
            (q1, q2)
        } else {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, -dx,  dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta,  dx,  dy);
            let q2 = Configuration::new(x, y, theta, 0.0);
            (q1, q2)
        };
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, CC_REGULAR, self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, CC_REGULAR, self.param);
        let length = cstart.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + cend.cc_turn_length(&q2);
        PathSlot { length, cstart: Some(cstart), cend: Some(cend),
                   qi1: Some(q1), qi2: Some(q2), ..PathSlot::infinite() }
    }

    fn tcsct_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Option<PathSlot> {
        if self.tcisce_exists(c1, c2) { Some(self.tciscT_path(c1, c2)) }
        else if self.tcesce_exists(c1, c2) { Some(self.tcescT_path(c1, c2)) }
        else { None }
    }
}

// ---------------------------------------------------------------------------
// Main path-evaluation function for one (c1, c2) pair
// ---------------------------------------------------------------------------

fn cc00_circles_rs_path(
    c1: &HcCcCircle,
    c2: &HcCcCircle,
    param: &HcCcCircleParam,
    kappa: f64,
    sigma: f64,
) -> HcCcRsPath {
    use HcCcRsPathType::*;

    let mut slots: Vec<PathSlot> = (0..NB_HC_CC_RS_PATHS).map(|_| PathSlot::infinite()).collect();

    let mut h = Cc00RsHelper::new(param);
    h.distance = center_distance(c1, c2);
    h.angle    = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);

    let mut skip = false;

    // E
    if configuration_equal(&c1.start, &c2.start) {
        slots[E as usize].length = 0.0;
        skip = true;
    }

    // S
    if !skip {
        if configuration_aligned(&c1.start, &c2.start) {
            slots[S as usize].length = configuration_distance(&c1.start, &c2.start);
            skip = true;
        } else if configuration_aligned(&c2.start, &c1.start) {
            slots[S as usize].length = configuration_distance(&c2.start, &c1.start);
            skip = true;
        }
    }

    // T
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
        kappa,
        sigma,
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
// Public state space
// ---------------------------------------------------------------------------

/// CC00 Reeds-Shepp state space — zero curvature at both start and end.
pub struct CcReedsSheppStateSpace {
    params_:         HcCcStateSpaceParams,
    discretization_: f64,
}

impl CcReedsSheppStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64) -> Self {
        Self {
            params_: HcCcStateSpaceParams::new(kappa, sigma),
            discretization_: discretization,
        }
    }

    fn cc00_reeds_shepp(&self, state1: &State, state2: &State) -> HcCcRsPath {
        let start = Configuration::new(state1.x, state1.y, state1.theta, 0.0);
        let end   = Configuration::new(state2.x, state2.y, state2.theta, 0.0);
        let p     = &self.params_.hc_cc_circle_param_;

        let start_circles = [
            HcCcCircle::from_configuration(start, true,  true,  CC_REGULAR, p),
            HcCcCircle::from_configuration(start, false, true,  CC_REGULAR, p),
            HcCcCircle::from_configuration(start, true,  false, CC_REGULAR, p),
            HcCcCircle::from_configuration(start, false, false, CC_REGULAR, p),
        ];
        let end_circles = [
            HcCcCircle::from_configuration(end, true,  true,  CC_REGULAR, p),
            HcCcCircle::from_configuration(end, false, true,  CC_REGULAR, p),
            HcCcCircle::from_configuration(end, true,  false, CC_REGULAR, p),
            HcCcCircle::from_configuration(end, false, false, CC_REGULAR, p),
        ];

        let mut best: Option<HcCcRsPath> = None;
        for sc in &start_circles {
            for ec in &end_circles {
                let path = cc00_circles_rs_path(sc, ec, p, self.params_.kappa_, self.params_.sigma_);
                if best.as_ref().map_or(true, |b| path.length < b.length) {
                    best = Some(path);
                }
            }
        }
        best.unwrap()
    }
}

impl StateSpace for CcReedsSheppStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let path = self.cc00_reeds_shepp(s1, s2);
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
            HcCcRsPathType::TT | HcCcRsPathType::TcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                cc_turn_controls(cs, q1, true,  &mut controls);
                cc_turn_controls(ce, q1, false, &mut controls);
            }
            HcCcRsPathType::TcTcT | HcCcRsPathType::TcTT | HcCcRsPathType::TTcT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true,  &mut controls);
                cc_turn_controls(ci, q2, true,  &mut controls);
                cc_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TST => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true,  &mut controls);
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
                cc_turn_controls(cs, q1, true,  &mut controls);
                straight_controls(q1, q2, &mut controls);
                cc_turn_controls(ci, q3, true,  &mut controls);
                cc_turn_controls(ce, q3, false, &mut controls);
            }
            HcCcRsPathType::TcTST => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                let q3 = path.qi3.as_ref().unwrap();
                cc_turn_controls(cs, q1, true,  &mut controls);
                cc_turn_controls(ci, q2, true,  &mut controls);
                straight_controls(q2, q3, &mut controls);
                cc_turn_controls(ce, q3, false, &mut controls);
            }
            HcCcRsPathType::TcTSTcT => {
                let cs  = path.cstart.as_ref().unwrap();
                let ce  = path.cend.as_ref().unwrap();
                let ci1 = path.ci1.as_ref().unwrap();
                let ci2 = path.ci2.as_ref().unwrap();
                let q1  = path.qi1.as_ref().unwrap();
                let q2  = path.qi2.as_ref().unwrap();
                let q3  = path.qi3.as_ref().unwrap();
                let q4  = path.qi4.as_ref().unwrap();
                cc_turn_controls(cs,  q1, true,  &mut controls);
                cc_turn_controls(ci1, q2, true,  &mut controls);
                straight_controls(q2, q3, &mut controls);
                cc_turn_controls(ci2, q4, true,  &mut controls);
                cc_turn_controls(ce,  q4, false, &mut controls);
            }
            HcCcRsPathType::TTcTT | HcCcRsPathType::TcTTcT => {
                let cs  = path.cstart.as_ref().unwrap();
                let ce  = path.cend.as_ref().unwrap();
                let ci1 = path.ci1.as_ref().unwrap();
                let ci2 = path.ci2.as_ref().unwrap();
                let q1  = path.qi1.as_ref().unwrap();
                let q2  = path.qi2.as_ref().unwrap();
                let q3  = path.qi3.as_ref().unwrap();
                cc_turn_controls(cs,  q1, true,  &mut controls);
                cc_turn_controls(ci1, q2, true,  &mut controls);
                cc_turn_controls(ci2, q3, true,  &mut controls);
                cc_turn_controls(ce,  q3, false, &mut controls);
            }
            HcCcRsPathType::TTT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let ci = path.ci1.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true,  &mut controls);
                cc_turn_controls(ci, q2, true,  &mut controls);
                cc_turn_controls(ce, q2, false, &mut controls);
            }
            HcCcRsPathType::TcST | HcCcRsPathType::TScT | HcCcRsPathType::TcScT => {
                let cs = path.cstart.as_ref().unwrap();
                let ce = path.cend.as_ref().unwrap();
                let q1 = path.qi1.as_ref().unwrap();
                let q2 = path.qi2.as_ref().unwrap();
                cc_turn_controls(cs, q1, true,  &mut controls);
                straight_controls(q1, q2, &mut controls);
                cc_turn_controls(ce, q2, false, &mut controls);
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
