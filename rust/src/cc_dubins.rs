use std::f64::consts::PI;

use crate::base_state_space::StateSpace;
use crate::configuration::{
    configuration_aligned, configuration_distance, configuration_equal, Configuration,
};
use crate::hc_cc_circle::{
    center_distance, configuration_on_hc_cc_circle, HcCcCircle, HcCcCircleParam,
};
use crate::hc_cc_state_space::HcCcStateSpaceParams;
use crate::paths::{
    cc_turn_controls, empty_controls, hc_turn_controls, reverse_control, rs_turn_controls,
    state_equal, straight_controls, subtract_control, CcDubinsPath, CcDubinsPathType,
    NB_CC_DUBINS_PATHS,
};
use crate::state::{Control, State};
use crate::utilities::{end_of_clothoid, get_epsilon, global_frame_change, sgn, HALF_PI};

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn array_index_min(arr: &[f64; NB_CC_DUBINS_PATHS]) -> usize {
    arr.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

fn idx_to_path_type(i: usize) -> CcDubinsPathType {
    match i {
        0 => CcDubinsPathType::E,
        1 => CcDubinsPathType::S,
        2 => CcDubinsPathType::T,
        3 => CcDubinsPathType::TT,
        4 => CcDubinsPathType::TST,
        5 => CcDubinsPathType::TTT,
        6 => CcDubinsPathType::TTTT,
        _ => CcDubinsPathType::E,
    }
}

fn make_rs_circle_param(kappa: f64) -> HcCcCircleParam {
    let mut p = HcCcCircleParam::default();
    p.set_param(kappa, f64::INFINITY, 1.0 / kappa, 0.0, 0.0, 1.0, 0.0);
    p
}

// Shared TiST tangent computation where c1 provides the radius/mu geometry.
fn tist_tangent_from_c1(
    c1: &HcCcCircle,
    c2: &HcCcCircle,
    radius: f64,
    sin_mu: f64,
    cos_mu: f64,
) -> (Configuration, Configuration) {
    let dist = center_distance(c1, c2);
    let angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
    let alpha = (2.0 * radius * cos_mu / dist).asin();
    let dx = radius * sin_mu;
    let dy = radius * cos_mu;
    if c1.left && c1.forward {
        let theta = angle + alpha;
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
        let q1 = Configuration::new(x, y, theta, 0.0);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
        (q1, Configuration::new(x, y, theta, 0.0))
    } else if c1.left && !c1.forward {
        let theta = angle - alpha;
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
        let q1 = Configuration::new(x, y, theta + PI, 0.0);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
        (q1, Configuration::new(x, y, theta + PI, 0.0))
    } else if !c1.left && c1.forward {
        let theta = angle - alpha;
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
        let q1 = Configuration::new(x, y, theta, 0.0);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
        (q1, Configuration::new(x, y, theta, 0.0))
    } else {
        let theta = angle + alpha;
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
        let q1 = Configuration::new(x, y, theta + PI, 0.0);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
        (q1, Configuration::new(x, y, theta + PI, 0.0))
    }
}

// Shared TeST tangent computation where c1 provides the radius/mu geometry.
fn test_tangent_from_c1(
    c1: &HcCcCircle,
    c2: &HcCcCircle,
    radius: f64,
    sin_mu: f64,
    cos_mu: f64,
) -> (Configuration, Configuration) {
    let dx = radius * sin_mu;
    let dy = radius * cos_mu;
    let theta = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
    if c1.left && c1.forward {
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
        let q1 = Configuration::new(x, y, theta, 0.0);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
        (q1, Configuration::new(x, y, theta, 0.0))
    } else if c1.left && !c1.forward {
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
        let q1 = Configuration::new(x, y, theta + PI, 0.0);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
        (q1, Configuration::new(x, y, theta + PI, 0.0))
    } else if !c1.left && c1.forward {
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
        let q1 = Configuration::new(x, y, theta, 0.0);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
        (q1, Configuration::new(x, y, theta, 0.0))
    } else {
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
        let q1 = Configuration::new(x, y, theta + PI, 0.0);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
        (q1, Configuration::new(x, y, theta + PI, 0.0))
    }
}

// ===========================================================================
// CC00 Helper
// ===========================================================================

struct Cc00Helper {
    distance: f64,
    angle: f64,
    param: HcCcCircleParam,
}

impl Cc00Helper {
    fn new(param: HcCcCircleParam) -> Self {
        Self { distance: 0.0, angle: 0.0, param }
    }

    // ----- TT ---------------------------------------------------------------
    fn tt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left == c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        (self.distance - 2.0 * c1.radius).abs() < get_epsilon()
    }

    fn tt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Configuration {
        let x = 0.5 * (c1.xc + c2.xc);
        let y = 0.5 * (c1.yc + c2.yc);
        let angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        let theta = if c1.left {
            if c1.forward { angle + HALF_PI - c1.mu } else { angle + HALF_PI + c1.mu }
        } else if c1.forward { angle - HALF_PI + c1.mu } else { angle - HALF_PI - c1.mu };
        Configuration::new(x, y, theta, 0.0)
    }

    fn tt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (f64, Configuration) {
        let q = self.tt_tangent_circles(c1, c2);
        (c1.cc_turn_length(&q) + c2.cc_turn_length(&q), q)
    }

    // ----- TST --------------------------------------------------------------
    fn tist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left == c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance >= 2.0 * c1.radius
    }

    fn test_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left != c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance >= 2.0 * c1.radius * c1.sin_mu
    }

    fn tst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tist_exists(c1, c2) || self.test_exists(c1, c2)
    }

    fn tist_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration) {
        tist_tangent_from_c1(c1, c2, c1.radius, c1.sin_mu, c1.cos_mu)
    }

    fn test_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration) {
        test_tangent_from_c1(c1, c2, c1.radius, c1.sin_mu, c1.cos_mu)
    }

    fn tist_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (f64, Configuration, Configuration) {
        let (q1, q2) = self.tist_tangent_circles(c1, c2);
        let l = c1.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + c2.cc_turn_length(&q2);
        (l, q1, q2)
    }

    fn test_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (f64, Configuration, Configuration) {
        let (q1, q2) = self.test_tangent_circles(c1, c2);
        let l = c1.cc_turn_length(&q1) + configuration_distance(&q1, &q2) + c2.cc_turn_length(&q2);
        (l, q1, q2)
    }

    fn tst_path(&self, c1: &HcCcCircle, c2: &HcCcCircle)
        -> (f64, Option<Configuration>, Option<Configuration>)
    {
        if self.tist_exists(c1, c2) {
            let (l, q1, q2) = self.tist_path(c1, c2);
            return (l, Some(q1), Some(q2));
        }
        if self.test_exists(c1, c2) {
            let (l, q1, q2) = self.test_path(c1, c2);
            return (l, Some(q1), Some(q2));
        }
        (f64::INFINITY, None, None)
    }

    // ----- TTT --------------------------------------------------------------
    fn ttt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left != c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance <= 4.0 * c1.radius
    }

    fn ttt_tangent_circles(
        &self, c1: &HcCcCircle, c2: &HcCcCircle, radius: f64,
    ) -> (Configuration, Configuration, Configuration, Configuration) {
        let theta = self.angle;
        let r = 2.0 * radius;
        let delta_x = 0.5 * self.distance;
        let delta_y = (r * r - delta_x * delta_x).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, &self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, &self.param);
        let qa = self.tt_tangent_circles(c1, &tgt1);
        let qb = self.tt_tangent_circles(&tgt1, c2);
        let qc = self.tt_tangent_circles(c1, &tgt2);
        let qd = self.tt_tangent_circles(&tgt2, c2);
        (qa, qb, qc, qd)
    }

    fn ttt_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (f64, Configuration, Configuration, HcCcCircle) {
        let (qa, qb, qc, qd) = self.ttt_tangent_circles(c1, c2, c1.radius);
        let middle1 = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, c1.regular, &self.param);
        let middle2 = HcCcCircle::from_configuration(qc, !c1.left, c1.forward, c1.regular, &self.param);
        let len1 = c1.cc_turn_length(&qa) + middle1.cc_turn_length(&qb) + c2.cc_turn_length(&qb);
        let len2 = c1.cc_turn_length(&qc) + middle2.cc_turn_length(&qd) + c2.cc_turn_length(&qd);
        if len1 < len2 { (len1, qa, qb, middle1) } else { (len2, qc, qd, middle2) }
    }
}

// ===========================================================================
// CC00_Dubins_State_Space
// ===========================================================================

pub struct CC00DubinsStateSpace {
    params: HcCcStateSpaceParams,
    discretization: f64,
    forward: bool,
}

impl CC00DubinsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64, forward: bool) -> Self {
        Self {
            params: HcCcStateSpaceParams::new(kappa, sigma),
            discretization,
            forward,
        }
    }

    fn cc00_circles_dubins_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> CcDubinsPath {
        let mut helper = Cc00Helper::new(self.params.hc_cc_circle_param_);
        helper.distance = center_distance(c1, c2);
        helper.angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);

        let mut length = [f64::INFINITY; NB_CC_DUBINS_PATHS];
        let mut qi1 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut qi2 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut cstart: [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] =
            std::array::from_fn(|_| None);
        let mut cend: [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] =
            std::array::from_fn(|_| None);
        let mut ci1: [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] =
            std::array::from_fn(|_| None);

        const E: usize = CcDubinsPathType::E as usize;
        const S: usize = CcDubinsPathType::S as usize;
        const T: usize = CcDubinsPathType::T as usize;
        const TT: usize = CcDubinsPathType::TT as usize;
        const TST: usize = CcDubinsPathType::TST as usize;
        const TTT: usize = CcDubinsPathType::TTT as usize;

        let p = &self.params.hc_cc_circle_param_;

        if configuration_equal(&c1.start, &c2.start) {
            length[E] = 0.0;
        } else if self.forward && configuration_aligned(&c1.start, &c2.start) {
            length[S] = configuration_distance(&c1.start, &c2.start);
        } else if !self.forward && configuration_aligned(&c2.start, &c1.start) {
            length[S] = configuration_distance(&c2.start, &c1.start);
        } else if configuration_on_hc_cc_circle(c1, &c2.start) {
            let cs = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, p);
            length[T] = cs.cc_turn_length(&c2.start);
            cstart[T] = Some(Box::new(cs));
        } else {
            // TT
            if helper.tt_exists(c1, c2) {
                let cs = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, p);
                let ce = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, p);
                let (l, q) = helper.tt_path(&cs, &ce);
                length[TT] = l;
                qi1[TT] = Some(q);
                cstart[TT] = Some(Box::new(cs));
                cend[TT] = Some(Box::new(ce));
            }
            // TST
            if helper.tst_exists(c1, c2) {
                let cs = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, p);
                let ce = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, p);
                let (l, q1, q2) = helper.tst_path(&cs, &ce);
                length[TST] = l;
                qi1[TST] = q1;
                qi2[TST] = q2;
                cstart[TST] = Some(Box::new(cs));
                cend[TST] = Some(Box::new(ce));
            }
            // TTT
            if helper.ttt_exists(c1, c2) {
                let cs = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, c1.regular, p);
                let ce = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, c2.regular, p);
                let (l, q1, q2, mid) = helper.ttt_path(&cs, &ce);
                length[TTT] = l;
                qi1[TTT] = Some(q1);
                qi2[TTT] = Some(q2);
                ci1[TTT] = Some(Box::new(mid));
                cstart[TTT] = Some(Box::new(cs));
                cend[TTT] = Some(Box::new(ce));
            }
        }

        let best = array_index_min(&length);
        CcDubinsPath::new(
            c1.start, c2.start, idx_to_path_type(best),
            self.params.kappa_, self.params.sigma_,
            qi1[best], qi2[best], None, None,
            cstart[best].take(), cend[best].take(), ci1[best].take(), None,
            length[best],
        )
    }

    fn cc00_dubins(&self, state1: &State, state2: &State) -> CcDubinsPath {
        let start = Configuration::new(state1.x, state1.y, state1.theta, 0.0);
        let end = Configuration::new(state2.x, state2.y, state2.theta, 0.0);
        let p = &self.params.hc_cc_circle_param_;

        let (start_circles, end_circles) = if self.forward {
            (
                [
                    HcCcCircle::from_configuration(start, true,  true,  true, p),
                    HcCcCircle::from_configuration(start, false, true,  true, p),
                ],
                [
                    HcCcCircle::from_configuration(end, true,  false, true, p),
                    HcCcCircle::from_configuration(end, false, false, true, p),
                ],
            )
        } else {
            (
                [
                    HcCcCircle::from_configuration(start, true,  false, true, p),
                    HcCcCircle::from_configuration(start, false, false, true, p),
                ],
                [
                    HcCcCircle::from_configuration(end, true,  true, true, p),
                    HcCcCircle::from_configuration(end, false, true, true, p),
                ],
            )
        };

        let mut best_path: Option<CcDubinsPath> = None;
        let mut best_len = f64::INFINITY;
        for i in 0..2 {
            for j in 0..2 {
                let path = self.cc00_circles_dubins_path(&start_circles[i], &end_circles[j]);
                if path.length < best_len {
                    best_len = path.length;
                    best_path = Some(path);
                }
            }
        }
        best_path.unwrap()
    }
}

impl StateSpace for CC00DubinsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let mut controls = Vec::new();
        let p = self.cc00_dubins(s1, s2);
        match p.path_type {
            CcDubinsPathType::E => empty_controls(&mut controls),
            CcDubinsPathType::S => straight_controls(&p.start, &p.end, &mut controls),
            CcDubinsPathType::T => {
                cc_turn_controls(p.cstart.as_ref().unwrap(), &p.end, true, &mut controls);
            }
            CcDubinsPathType::TT => {
                let qi1 = p.qi1.as_ref().unwrap();
                cc_turn_controls(p.cstart.as_ref().unwrap(), qi1, true,  &mut controls);
                cc_turn_controls(p.cend.as_ref().unwrap(),   qi1, false, &mut controls);
            }
            CcDubinsPathType::TST => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                cc_turn_controls(p.cstart.as_ref().unwrap(), qi1, true,  &mut controls);
                straight_controls(qi1, qi2, &mut controls);
                cc_turn_controls(p.cend.as_ref().unwrap(),   qi2, false, &mut controls);
            }
            CcDubinsPathType::TTT => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                cc_turn_controls(p.cstart.as_ref().unwrap(), qi1, true,  &mut controls);
                cc_turn_controls(p.ci1.as_ref().unwrap(),    qi2, true,  &mut controls);
                cc_turn_controls(p.cend.as_ref().unwrap(),   qi2, false, &mut controls);
            }
            _ => {}
        }
        controls
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        vec![self.get_controls(s1, s2)]
    }

    fn discretization(&self) -> f64 { self.discretization }
}

// ===========================================================================
// CC0pm Helper
// ===========================================================================

struct Cc0pmHelper {
    distance: f64,
    angle: f64,
    param: HcCcCircleParam,
}

impl Cc0pmHelper {
    fn new(param: HcCcCircleParam) -> Self {
        Self { distance: 0.0, angle: 0.0, param }
    }

    // ----- TT ---------------------------------------------------------------
    fn tt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left == c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        (self.distance - 2.0 * c1.radius).abs() < get_epsilon()
    }

    // CC0pm uses c1.mu (same as CC00)
    fn tt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Configuration {
        let x = 0.5 * (c1.xc + c2.xc);
        let y = 0.5 * (c1.yc + c2.yc);
        let angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        let theta = if c1.left {
            if c1.forward { angle + HALF_PI - c1.mu } else { angle + HALF_PI + c1.mu }
        } else if c1.forward { angle - HALF_PI + c1.mu } else { angle - HALF_PI - c1.mu };
        Configuration::new(x, y, theta, 0.0)
    }

    fn tt_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration) {
        let q1 = self.tt_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, true, &self.param);
        let cend   = HcCcCircle::from_configuration(q1, c2.left, !c2.forward, true, &self.param);
        let q2 = c2.start; // same kappa as c2
        let length = cstart.cc_turn_length(&q1) + cend.hc_turn_length(&q2);
        (length, cstart, cend, q1, q2)
    }

    // ----- TST --------------------------------------------------------------
    fn tist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left == c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance >= 2.0 * c1.radius
    }

    fn test_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left != c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance >= 2.0 * c1.radius * c1.sin_mu
    }

    fn tst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tist_exists(c1, c2) || self.test_exists(c1, c2)
    }

    fn tist_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration) {
        tist_tangent_from_c1(c1, c2, c1.radius, c1.sin_mu, c1.cos_mu)
    }

    fn test_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration) {
        test_tangent_from_c1(c1, c2, c1.radius, c1.sin_mu, c1.cos_mu)
    }

    fn tist_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration) {
        let (q1, q2) = self.tist_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, true, &self.param);
        let cend   = HcCcCircle::from_configuration(q2, c2.left, !c2.forward, true, &self.param);
        let q3 = c2.start;
        let length = cstart.cc_turn_length(&q1)
            + configuration_distance(&q1, &q2)
            + cend.hc_turn_length(&q3);
        (length, cstart, cend, q1, q2, q3)
    }

    fn test_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration) {
        let (q1, q2) = self.test_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, true, &self.param);
        let cend   = HcCcCircle::from_configuration(q2, c2.left, !c2.forward, true, &self.param);
        let q3 = c2.start;
        let length = cstart.cc_turn_length(&q1)
            + configuration_distance(&q1, &q2)
            + cend.hc_turn_length(&q3);
        (length, cstart, cend, q1, q2, q3)
    }

    #[allow(clippy::type_complexity)]
    fn tst_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> Option<(f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration)> {
        if self.tist_exists(c1, c2) { return Some(self.tist_path(c1, c2)); }
        if self.test_exists(c1, c2)  { return Some(self.test_path(c1, c2)); }
        None
    }

    // ----- TTT --------------------------------------------------------------
    fn ttt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left != c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance <= 4.0 * c1.radius
    }

    fn ttt_tangent_circles(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (Configuration, Configuration, Configuration, Configuration) {
        let theta = self.angle;
        let r = 2.0 * c1.radius;
        let delta_x = 0.5 * self.distance;
        let delta_y = (r * r - delta_x * delta_x).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, &self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, &self.param);
        let qa = self.tt_tangent_circles(c1,    &tgt1);
        let qb = self.tt_tangent_circles(&tgt1, c2);
        let qc = self.tt_tangent_circles(c1,    &tgt2);
        let qd = self.tt_tangent_circles(&tgt2, c2);
        (qa, qb, qc, qd)
    }

    #[allow(clippy::type_complexity)]
    fn ttt_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration, HcCcCircle) {
        let (qa, qb, qc, qd) = self.ttt_tangent_circles(c1, c2);
        let middle1 = HcCcCircle::from_configuration(qa, !c1.left, c1.forward, true, &self.param);
        let end1    = HcCcCircle::from_configuration(qb, c2.left, !c2.forward, true, &self.param);
        let middle2 = HcCcCircle::from_configuration(qc, !c1.left, c1.forward, true, &self.param);
        let end2    = HcCcCircle::from_configuration(qd, c2.left, !c2.forward, true, &self.param);
        let cstart  = HcCcCircle::from_configuration(c1.start, c1.left, c1.forward, true, &self.param);
        let q3 = c2.start;
        let len1 = cstart.cc_turn_length(&qa) + middle1.cc_turn_length(&qb) + end1.hc_turn_length(&q3);
        let len2 = cstart.cc_turn_length(&qc) + middle2.cc_turn_length(&qd) + end2.hc_turn_length(&q3);
        if len1 < len2 {
            (len1, cstart, end1, qa, qb, q3, middle1)
        } else {
            (len2, cstart, end2, qc, qd, q3, middle2)
        }
    }
}

// ===========================================================================
// CC0pm_Dubins_State_Space
// ===========================================================================

pub struct CC0pmDubinsStateSpace {
    params: HcCcStateSpaceParams,
    rs_param: HcCcCircleParam,
    discretization: f64,
    forward: bool,
}

impl CC0pmDubinsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64, forward: bool) -> Self {
        Self {
            params: HcCcStateSpaceParams::new(kappa, sigma),
            rs_param: make_rs_circle_param(kappa),
            discretization,
            forward,
        }
    }

    fn cc0pm_circles_dubins_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> CcDubinsPath {
        let mut helper = Cc0pmHelper::new(self.params.hc_cc_circle_param_);
        helper.distance = center_distance(c1, c2);
        helper.angle    = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);

        let mut length = [f64::INFINITY; NB_CC_DUBINS_PATHS];
        let mut qi1 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut qi2 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut qi3 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut cstart: [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] = std::array::from_fn(|_| None);
        let mut cend:   [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] = std::array::from_fn(|_| None);
        let mut ci1:    [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] = std::array::from_fn(|_| None);

        const E:   usize = CcDubinsPathType::E   as usize;
        const T:   usize = CcDubinsPathType::T   as usize;
        const TT:  usize = CcDubinsPathType::TT  as usize;
        const TST: usize = CcDubinsPathType::TST as usize;
        const TTT: usize = CcDubinsPathType::TTT as usize;

        if configuration_equal(&c1.start, &c2.start) {
            length[E] = 0.0;
        } else if helper.distance < get_epsilon() {
            // centers coincide: single HC turn from c1 to c2.start
            let cs = HcCcCircle::from_configuration(
                c1.start, c1.left, c1.forward, true, &self.params.hc_cc_circle_param_,
            );
            length[T] = cs.hc_turn_length(&c2.start);
            cstart[T] = Some(Box::new(cs));
        } else {
            if helper.tt_exists(c1, c2) {
                let (l, cs, ce, q1, q2) = helper.tt_path(c1, c2);
                length[TT] = l;
                qi1[TT] = Some(q1);
                qi2[TT] = Some(q2);
                cstart[TT] = Some(Box::new(cs));
                cend[TT]   = Some(Box::new(ce));
            }
            if helper.tst_exists(c1, c2) {
                if let Some((l, cs, ce, q1, q2, q3)) = helper.tst_path(c1, c2) {
                    length[TST] = l;
                    qi1[TST] = Some(q1);
                    qi2[TST] = Some(q2);
                    qi3[TST] = Some(q3);
                    cstart[TST] = Some(Box::new(cs));
                    cend[TST]   = Some(Box::new(ce));
                }
            }
            if helper.ttt_exists(c1, c2) {
                let (l, cs, ce, q1, q2, q3, mid) = helper.ttt_path(c1, c2);
                length[TTT] = l;
                qi1[TTT] = Some(q1);
                qi2[TTT] = Some(q2);
                qi3[TTT] = Some(q3);
                ci1[TTT] = Some(Box::new(mid));
                cstart[TTT] = Some(Box::new(cs));
                cend[TTT]   = Some(Box::new(ce));
            }
        }

        let best = array_index_min(&length);
        CcDubinsPath::new(
            c1.start, c2.start, idx_to_path_type(best),
            self.params.kappa_, self.params.sigma_,
            qi1[best], qi2[best], qi3[best], None,
            cstart[best].take(), cend[best].take(), ci1[best].take(), None,
            length[best],
        )
    }

    fn cc0pm_dubins(&self, state1: &State, state2: &State) -> CcDubinsPath {
        let start = Configuration::new(state1.x, state1.y, state1.theta, 0.0);
        let end1  = Configuration::new(state2.x, state2.y, state2.theta,  self.params.kappa_);
        let end2  = Configuration::new(state2.x, state2.y, state2.theta, -self.params.kappa_);
        let hp = &self.params.hc_cc_circle_param_;
        let rp = &self.rs_param;

        let (start_circles, end_circles) = if self.forward {
            (
                [
                    HcCcCircle::from_configuration(start, true,  true, true, hp),
                    HcCcCircle::from_configuration(start, false, true, true, hp),
                ],
                [
                    HcCcCircle::from_configuration(end1, true,  false, true, rp),
                    HcCcCircle::from_configuration(end2, false, false, true, rp),
                ],
            )
        } else {
            (
                [
                    HcCcCircle::from_configuration(start, true,  false, true, hp),
                    HcCcCircle::from_configuration(start, false, false, true, hp),
                ],
                [
                    HcCcCircle::from_configuration(end1, true,  true, true, rp),
                    HcCcCircle::from_configuration(end2, false, true, true, rp),
                ],
            )
        };

        let mut best_path: Option<CcDubinsPath> = None;
        let mut best_len = f64::INFINITY;
        for i in 0..2 {
            for j in 0..2 {
                // skip circle at the end for curvature continuity
                if j == 0 && state2.kappa < 0.0 { continue; }
                if j == 1 && state2.kappa > 0.0 { continue; }
                let path = self.cc0pm_circles_dubins_path(&start_circles[i], &end_circles[j]);
                if path.length < best_len {
                    best_len = path.length;
                    best_path = Some(path);
                }
            }
        }
        best_path.unwrap()
    }
}

impl StateSpace for CC0pmDubinsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let mut controls = Vec::new();
        let p = self.cc0pm_dubins(s1, s2);
        match p.path_type {
            CcDubinsPathType::E => empty_controls(&mut controls),
            CcDubinsPathType::T => {
                hc_turn_controls(p.cstart.as_ref().unwrap(), &p.end, true, &mut controls);
            }
            CcDubinsPathType::TT => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                cc_turn_controls(p.cstart.as_ref().unwrap(), qi1, true, &mut controls);
                hc_turn_controls(p.cend.as_ref().unwrap(),   qi2, true, &mut controls);
            }
            CcDubinsPathType::TST => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                let qi3 = p.qi3.as_ref().unwrap();
                cc_turn_controls(p.cstart.as_ref().unwrap(), qi1, true, &mut controls);
                straight_controls(qi1, qi2, &mut controls);
                hc_turn_controls(p.cend.as_ref().unwrap(),   qi3, true, &mut controls);
            }
            CcDubinsPathType::TTT => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                let qi3 = p.qi3.as_ref().unwrap();
                cc_turn_controls(p.cstart.as_ref().unwrap(), qi1, true, &mut controls);
                cc_turn_controls(p.ci1.as_ref().unwrap(),    qi2, true, &mut controls);
                hc_turn_controls(p.cend.as_ref().unwrap(),   qi3, true, &mut controls);
            }
            _ => {}
        }
        controls
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        vec![self.get_controls(s1, s2)]
    }

    fn discretization(&self) -> f64 { self.discretization }
}

// ===========================================================================
// CCpm0 Helper
// ===========================================================================

struct Ccpm0Helper {
    distance: f64,
    angle: f64,
    param: HcCcCircleParam,
}

impl Ccpm0Helper {
    fn new(param: HcCcCircleParam) -> Self {
        Self { distance: 0.0, angle: 0.0, param }
    }

    // ----- TT ---------------------------------------------------------------
    fn tt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left == c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        (self.distance - 2.0 * c2.radius).abs() < get_epsilon()
    }

    // CCpm0 uses c2.mu
    fn tt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Configuration {
        let x = 0.5 * (c1.xc + c2.xc);
        let y = 0.5 * (c1.yc + c2.yc);
        let angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        let theta = if c1.left {
            if c1.forward { angle + HALF_PI - c2.mu } else { angle + HALF_PI + c2.mu }
        } else if c1.forward { angle - HALF_PI + c2.mu } else { angle - HALF_PI - c2.mu };
        Configuration::new(x, y, theta, 0.0)
    }

    fn tt_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration) {
        let q2 = self.tt_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(q2, !c2.left, c2.forward, true, &self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, true, &self.param);
        let q1 = c1.start; // kappa already set in c1.start
        let length = cstart.hc_turn_length(&q1) + cend.cc_turn_length(&q2);
        (length, cstart, cend, q1, q2)
    }

    // ----- TST --------------------------------------------------------------
    fn tist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left == c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance >= 2.0 * c2.radius
    }

    fn test_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left != c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance >= 2.0 * c2.radius * c2.sin_mu
    }

    fn tst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tist_exists(c1, c2) || self.test_exists(c1, c2)
    }

    // CCpm0 tist uses c2's geometry (passing c2 as source for tist_tangent_from_c1 is wrong;
    // we need to pass c2's geometry but with c1's position context).
    fn tist_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration) {
        // alpha uses c2.radius/c2.cos_mu; delta uses c2.radius/c2.sin_mu/c2.cos_mu
        let dist = center_distance(c1, c2);
        let angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        let alpha = (2.0 * c2.radius * c2.cos_mu / dist).asin();
        let dx = c2.radius * c2.sin_mu;
        let dy = c2.radius * c2.cos_mu;
        if c1.left && c1.forward {
            let theta = angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
            (q1, Configuration::new(x, y, theta, 0.0))
        } else if c1.left && !c1.forward {
            let theta = angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            (q1, Configuration::new(x, y, theta + PI, 0.0))
        } else if !c1.left && c1.forward {
            let theta = angle - alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            (q1, Configuration::new(x, y, theta, 0.0))
        } else {
            let theta = angle + alpha;
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
            (q1, Configuration::new(x, y, theta + PI, 0.0))
        }
    }

    fn test_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration) {
        let dx = c2.radius * c2.sin_mu;
        let dy = c2.radius * c2.cos_mu;
        let theta = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        if c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            (q1, Configuration::new(x, y, theta, 0.0))
        } else if c1.left && !c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
            (q1, Configuration::new(x, y, theta + PI, 0.0))
        } else if !c1.left && c1.forward {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, dy);
            let q1 = Configuration::new(x, y, theta, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, dy);
            (q1, Configuration::new(x, y, theta, 0.0))
        } else {
            let (x, y) = global_frame_change(c1.xc, c1.yc, theta, dx, -dy);
            let q1 = Configuration::new(x, y, theta + PI, 0.0);
            let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -dx, -dy);
            (q1, Configuration::new(x, y, theta + PI, 0.0))
        }
    }

    fn tist_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration) {
        let (q2, q3) = self.tist_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, true, &self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, true, &self.param);
        let q1 = c1.start;
        let length = cstart.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3)
            + cend.cc_turn_length(&q3);
        (length, cstart, cend, q1, q2, q3)
    }

    fn test_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration) {
        let (q2, q3) = self.test_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, true, &self.param);
        let cend   = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, true, &self.param);
        let q1 = c1.start;
        let length = cstart.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3)
            + cend.cc_turn_length(&q3);
        (length, cstart, cend, q1, q2, q3)
    }

    #[allow(clippy::type_complexity)]
    fn tst_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> Option<(f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration)> {
        if self.tist_exists(c1, c2) { return Some(self.tist_path(c1, c2)); }
        if self.test_exists(c1, c2)  { return Some(self.test_path(c1, c2)); }
        None
    }

    // ----- TTT --------------------------------------------------------------
    fn ttt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left != c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance <= 4.0 * c2.radius
    }

    fn ttt_tangent_circles(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (Configuration, Configuration, Configuration, Configuration) {
        let theta = self.angle;
        let r = 2.0 * c2.radius;
        let delta_x = 0.5 * self.distance;
        let delta_y = (r * r - delta_x * delta_x).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, &self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, &self.param);
        let qa = self.tt_tangent_circles(c1,    &tgt1);
        let qb = self.tt_tangent_circles(&tgt1, c2);
        let qc = self.tt_tangent_circles(c1,    &tgt2);
        let qd = self.tt_tangent_circles(&tgt2, c2);
        (qa, qb, qc, qd)
    }

    #[allow(clippy::type_complexity)]
    fn ttt_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, HcCcCircle) {
        let (qa, qb, qc, qd) = self.ttt_tangent_circles(c1, c2);
        let start1  = HcCcCircle::from_configuration(qa, c1.left,  !c1.forward, true, &self.param);
        let middle1 = HcCcCircle::from_configuration(qa, !c1.left,  c1.forward, true, &self.param);
        let start2  = HcCcCircle::from_configuration(qc, c1.left,  !c1.forward, true, &self.param);
        let middle2 = HcCcCircle::from_configuration(qc, !c1.left,  c1.forward, true, &self.param);
        let cend    = HcCcCircle::from_configuration(c2.start, c2.left, c2.forward, true, &self.param);
        let q1 = c1.start;
        let len1 = start1.hc_turn_length(&q1) + middle1.cc_turn_length(&qb) + cend.cc_turn_length(&qb);
        let len2 = start2.hc_turn_length(&q1) + middle2.cc_turn_length(&qd) + cend.cc_turn_length(&qd);
        if len1 < len2 {
            (len1, start1, cend, q1, qb, middle1)
        } else {
            (len2, start2, cend, q1, qd, middle2)
        }
    }
}

// ===========================================================================
// CCpm0_Dubins_State_Space
// ===========================================================================

pub struct CCpm0DubinsStateSpace {
    params: HcCcStateSpaceParams,
    rs_param: HcCcCircleParam,
    discretization: f64,
    forward: bool,
}

impl CCpm0DubinsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64, forward: bool) -> Self {
        Self {
            params: HcCcStateSpaceParams::new(kappa, sigma),
            rs_param: make_rs_circle_param(kappa),
            discretization,
            forward,
        }
    }

    fn ccpm0_circles_dubins_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> CcDubinsPath {
        let mut helper = Ccpm0Helper::new(self.params.hc_cc_circle_param_);
        helper.distance = center_distance(c1, c2);
        helper.angle    = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);

        let mut length = [f64::INFINITY; NB_CC_DUBINS_PATHS];
        let mut qi1 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut qi2 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut qi3 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut cstart: [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] = std::array::from_fn(|_| None);
        let mut cend:   [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] = std::array::from_fn(|_| None);
        let mut ci1:    [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] = std::array::from_fn(|_| None);

        const E:   usize = CcDubinsPathType::E   as usize;
        const T:   usize = CcDubinsPathType::T   as usize;
        const TT:  usize = CcDubinsPathType::TT  as usize;
        const TST: usize = CcDubinsPathType::TST as usize;
        const TTT: usize = CcDubinsPathType::TTT as usize;

        if configuration_equal(&c1.start, &c2.start) {
            length[E] = 0.0;
        } else if helper.distance < get_epsilon() {
            let ce = HcCcCircle::from_configuration(
                c2.start, c2.left, c2.forward, true, &self.params.hc_cc_circle_param_,
            );
            length[T] = ce.hc_turn_length(&c1.start);
            cend[T] = Some(Box::new(ce));
        } else {
            if helper.tt_exists(c1, c2) {
                let (l, cs, ce, q1, q2) = helper.tt_path(c1, c2);
                length[TT] = l;
                qi1[TT] = Some(q1);
                qi2[TT] = Some(q2);
                cstart[TT] = Some(Box::new(cs));
                cend[TT]   = Some(Box::new(ce));
            }
            if helper.tst_exists(c1, c2) {
                if let Some((l, cs, ce, q1, q2, q3)) = helper.tst_path(c1, c2) {
                    length[TST] = l;
                    qi1[TST] = Some(q1);
                    qi2[TST] = Some(q2);
                    qi3[TST] = Some(q3);
                    cstart[TST] = Some(Box::new(cs));
                    cend[TST]   = Some(Box::new(ce));
                }
            }
            if helper.ttt_exists(c1, c2) {
                let (l, cs, ce, q1, q2, mid) = helper.ttt_path(c1, c2);
                length[TTT] = l;
                qi1[TTT] = Some(q1);
                qi2[TTT] = Some(q2);
                ci1[TTT] = Some(Box::new(mid));
                cstart[TTT] = Some(Box::new(cs));
                cend[TTT]   = Some(Box::new(ce));
            }
        }

        let best = array_index_min(&length);
        CcDubinsPath::new(
            c1.start, c2.start, idx_to_path_type(best),
            self.params.kappa_, self.params.sigma_,
            qi1[best], qi2[best], qi3[best], None,
            cstart[best].take(), cend[best].take(), ci1[best].take(), None,
            length[best],
        )
    }

    fn ccpm0_dubins(&self, state1: &State, state2: &State) -> CcDubinsPath {
        let start1 = Configuration::new(state1.x, state1.y, state1.theta,  self.params.kappa_);
        let start2 = Configuration::new(state1.x, state1.y, state1.theta, -self.params.kappa_);
        let end    = Configuration::new(state2.x, state2.y, state2.theta, 0.0);
        let hp = &self.params.hc_cc_circle_param_;
        let rp = &self.rs_param;

        let (start_circles, end_circles) = if self.forward {
            (
                [
                    HcCcCircle::from_configuration(start1, true,  true, true, rp),
                    HcCcCircle::from_configuration(start2, false, true, true, rp),
                ],
                [
                    HcCcCircle::from_configuration(end, true,  false, true, hp),
                    HcCcCircle::from_configuration(end, false, false, true, hp),
                ],
            )
        } else {
            (
                [
                    HcCcCircle::from_configuration(start1, true,  false, true, rp),
                    HcCcCircle::from_configuration(start2, false, false, true, rp),
                ],
                [
                    HcCcCircle::from_configuration(end, true,  true, true, hp),
                    HcCcCircle::from_configuration(end, false, true, true, hp),
                ],
            )
        };

        let mut best_path: Option<CcDubinsPath> = None;
        let mut best_len = f64::INFINITY;
        for i in 0..2 {
            // skip circle at the beginning for curvature continuity
            if i == 0 && state1.kappa < 0.0 { continue; }
            if i == 1 && state1.kappa > 0.0 { continue; }
            for j in 0..2 {
                let path = self.ccpm0_circles_dubins_path(&start_circles[i], &end_circles[j]);
                if path.length < best_len {
                    best_len = path.length;
                    best_path = Some(path);
                }
            }
        }
        best_path.unwrap()
    }
}

impl StateSpace for CCpm0DubinsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let mut controls = Vec::new();
        let p = self.ccpm0_dubins(s1, s2);
        match p.path_type {
            CcDubinsPathType::E => empty_controls(&mut controls),
            CcDubinsPathType::T => {
                hc_turn_controls(p.cend.as_ref().unwrap(), &p.start, false, &mut controls);
            }
            CcDubinsPathType::TT => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                hc_turn_controls(p.cstart.as_ref().unwrap(), qi1, false, &mut controls);
                cc_turn_controls(p.cend.as_ref().unwrap(),   qi2, false, &mut controls);
            }
            CcDubinsPathType::TST => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                let qi3 = p.qi3.as_ref().unwrap();
                hc_turn_controls(p.cstart.as_ref().unwrap(), qi1, false, &mut controls);
                straight_controls(qi2, qi3, &mut controls);
                cc_turn_controls(p.cend.as_ref().unwrap(),   qi3, false, &mut controls);
            }
            CcDubinsPathType::TTT => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                hc_turn_controls(p.cstart.as_ref().unwrap(), qi1, false, &mut controls);
                cc_turn_controls(p.ci1.as_ref().unwrap(),    qi2, true,  &mut controls);
                cc_turn_controls(p.cend.as_ref().unwrap(),   qi2, false, &mut controls);
            }
            _ => {}
        }
        controls
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        vec![self.get_controls(s1, s2)]
    }

    fn discretization(&self) -> f64 { self.discretization }
}

// ===========================================================================
// CCpmpm Helper
// ===========================================================================

struct CcpmpmHelper {
    distance: f64,
    angle: f64,
    /// hc_cc_circle_param_ — used for creating CC circles and for mu/radius
    param: HcCcCircleParam,
}

impl CcpmpmHelper {
    fn new(param: HcCcCircleParam) -> Self {
        Self { distance: 0.0, angle: 0.0, param }
    }

    // ----- TT ---------------------------------------------------------------
    fn tt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left == c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        (self.distance - 2.0 * self.param.radius).abs() < get_epsilon()
    }

    // CCpmpm uses self.param.mu (the CC mu, not the RS circle's mu=0)
    fn tt_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> Configuration {
        let x = 0.5 * (c1.xc + c2.xc);
        let y = 0.5 * (c1.yc + c2.yc);
        let angle = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);
        let mu = self.param.mu;
        let theta = if c1.left {
            if c1.forward { angle + HALF_PI - mu } else { angle + HALF_PI + mu }
        } else if c1.forward { angle - HALF_PI + mu } else { angle - HALF_PI - mu };
        Configuration::new(x, y, theta, 0.0)
    }

    fn tt_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration) {
        let q2 = self.tt_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, true, &self.param);
        let cend   = HcCcCircle::from_configuration(q2, c2.left, !c2.forward, true, &self.param);
        let q1 = c1.start;
        let q3 = c2.start;
        let length = cstart.hc_turn_length(&q1) + cend.hc_turn_length(&q3);
        (length, cstart, cend, q1, q2, q3)
    }

    // ----- TST --------------------------------------------------------------
    fn tist_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left == c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance >= 2.0 * self.param.radius
    }

    fn test_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left != c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance >= 2.0 * self.param.radius * self.param.sin_mu
    }

    fn tst_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        self.tist_exists(c1, c2) || self.test_exists(c1, c2)
    }

    fn tist_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration) {
        tist_tangent_from_c1(c1, c2, self.param.radius, self.param.sin_mu, self.param.cos_mu)
    }

    fn test_tangent_circles(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> (Configuration, Configuration) {
        test_tangent_from_c1(c1, c2, self.param.radius, self.param.sin_mu, self.param.cos_mu)
    }

    fn tist_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration, Configuration) {
        let (q2, q3) = self.tist_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, true, &self.param);
        let cend   = HcCcCircle::from_configuration(q3, c2.left, !c2.forward, true, &self.param);
        let q1 = c1.start;
        let q4 = c2.start;
        let length = cstart.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3)
            + cend.hc_turn_length(&q4);
        (length, cstart, cend, q1, q2, q3, q4)
    }

    fn test_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration, Configuration) {
        let (q2, q3) = self.test_tangent_circles(c1, c2);
        let cstart = HcCcCircle::from_configuration(q2, c1.left, !c1.forward, true, &self.param);
        let cend   = HcCcCircle::from_configuration(q3, c2.left, !c2.forward, true, &self.param);
        let q1 = c1.start;
        let q4 = c2.start;
        let length = cstart.hc_turn_length(&q1)
            + configuration_distance(&q2, &q3)
            + cend.hc_turn_length(&q4);
        (length, cstart, cend, q1, q2, q3, q4)
    }

    #[allow(clippy::type_complexity)]
    fn tst_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> Option<(f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration, Configuration)>
    {
        if self.tist_exists(c1, c2) { return Some(self.tist_path(c1, c2)); }
        if self.test_exists(c1, c2)  { return Some(self.test_path(c1, c2)); }
        None
    }

    // ----- TTT --------------------------------------------------------------
    fn ttt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left != c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance <= 4.0 * self.param.radius
    }

    fn ttt_tangent_circles(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (Configuration, Configuration, Configuration, Configuration) {
        let theta = self.angle;
        let r = 2.0 * self.param.radius;
        let delta_x = 0.5 * self.distance;
        let delta_y = (r * r - delta_x * delta_x).sqrt();
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, &self.param);
        let (x, y) = global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, &self.param);
        let qa = self.tt_tangent_circles(c1,    &tgt1);
        let qb = self.tt_tangent_circles(&tgt1, c2);
        let qc = self.tt_tangent_circles(c1,    &tgt2);
        let qd = self.tt_tangent_circles(&tgt2, c2);
        (qa, qb, qc, qd)
    }

    #[allow(clippy::type_complexity)]
    fn ttt_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration, HcCcCircle) {
        let (qa, qb, qc, qd) = self.ttt_tangent_circles(c1, c2);
        let start1  = HcCcCircle::from_configuration(qa, c1.left,  !c1.forward, true, &self.param);
        let middle1 = HcCcCircle::from_configuration(qa, !c1.left,  c1.forward, true, &self.param);
        let end1    = HcCcCircle::from_configuration(qb, c2.left,  !c2.forward, true, &self.param);
        let start2  = HcCcCircle::from_configuration(qc, c1.left,  !c1.forward, true, &self.param);
        let middle2 = HcCcCircle::from_configuration(qc, !c1.left,  c1.forward, true, &self.param);
        let end2    = HcCcCircle::from_configuration(qd, c2.left,  !c2.forward, true, &self.param);
        let q1 = c1.start;
        let q3 = c2.start;
        let len1 = start1.hc_turn_length(&q1) + middle1.cc_turn_length(&qb) + end1.hc_turn_length(&q3);
        let len2 = start2.hc_turn_length(&q1) + middle2.cc_turn_length(&qd) + end2.hc_turn_length(&q3);
        if len1 < len2 {
            (len1, start1, end1, q1, qb, q3, middle1)
        } else {
            (len2, start2, end2, q1, qd, q3, middle2)
        }
    }

    // ----- TTTT -------------------------------------------------------------
    fn tttt_exists(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> bool {
        if c1.left == c2.left { return false; }
        if c1.forward == c2.forward { return false; }
        self.distance <= 6.0 * self.param.radius
    }

    #[allow(clippy::type_complexity)]
    fn tttt_tangent_circles(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (Configuration, Configuration, Configuration, Configuration, Configuration, Configuration) {
        let theta = self.angle;
        let r1 = 2.0 * self.param.radius;
        let (delta_x, delta_y) = if self.distance < r1 {
            let dx = (self.distance + r1) / 2.0;
            let dy = (r1 * r1 - dx * dx).sqrt();
            (dx, dy)
        } else {
            let dx = (self.distance - r1) / 2.0;
            let dy = (r1 * r1 - dx * dx).sqrt();
            (dx, dy)
        };

        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x,  delta_y);
        let tgt1 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, &self.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x,  delta_y);
        let tgt2 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward, c2.regular, &self.param);

        let (x, y) = global_frame_change(c1.xc, c1.yc, theta,  delta_x, -delta_y);
        let tgt3 = HcCcCircle::from_center(x, y, !c1.left, c1.forward, c1.regular, &self.param);
        let (x, y) = global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y);
        let tgt4 = HcCcCircle::from_center(x, y, !c2.left, !c2.forward, c2.regular, &self.param);

        let qa = self.tt_tangent_circles(c1,    &tgt1);
        let qb = self.tt_tangent_circles(&tgt1, &tgt2);
        let qc = self.tt_tangent_circles(&tgt2, c2);

        let qd = self.tt_tangent_circles(c1,    &tgt3);
        let qe = self.tt_tangent_circles(&tgt3, &tgt4);
        let qf = self.tt_tangent_circles(&tgt4, c2);
        (qa, qb, qc, qd, qe, qf)
    }

    #[allow(clippy::type_complexity)]
    fn tttt_path(
        &self, c1: &HcCcCircle, c2: &HcCcCircle,
    ) -> (f64, HcCcCircle, HcCcCircle, Configuration, Configuration, Configuration, HcCcCircle, HcCcCircle) {
        let (qa, qb, qc, qd, qe, qf) = self.tttt_tangent_circles(c1, c2);

        let start1   = HcCcCircle::from_configuration(qa, c1.left,  !c1.forward, true, &self.param);
        let middle1  = HcCcCircle::from_configuration(qa, !c1.left,  c1.forward, true, &self.param);
        let middle2  = HcCcCircle::from_configuration(qc, !c2.left,  c2.forward, true, &self.param);
        let end1     = HcCcCircle::from_configuration(qc, c2.left,  !c2.forward, true, &self.param);
        let start2   = HcCcCircle::from_configuration(qd, c1.left,  !c1.forward, true, &self.param);
        let middle3  = HcCcCircle::from_configuration(qd, !c1.left,  c1.forward, true, &self.param);
        let middle4  = HcCcCircle::from_configuration(qf, !c2.left,  c2.forward, true, &self.param);
        let end2     = HcCcCircle::from_configuration(qf, c2.left,  !c2.forward, true, &self.param);

        let q1 = c1.start;
        let q3 = c2.start;

        let len1 = start1.hc_turn_length(&q1)
            + middle1.cc_turn_length(&qb)
            + middle2.cc_turn_length(&qb)
            + end1.hc_turn_length(&q3);
        let len2 = start2.hc_turn_length(&q1)
            + middle3.cc_turn_length(&qe)
            + middle4.cc_turn_length(&qe)
            + end2.hc_turn_length(&q3);

        if len1 < len2 {
            (len1, start1, end1, q1, qb, q3, middle1, middle2)
        } else {
            (len2, start2, end2, q1, qe, q3, middle3, middle4)
        }
    }
}

// ===========================================================================
// CCpmpm_Dubins_State_Space
// ===========================================================================

pub struct CCpmpmDubinsStateSpace {
    params: HcCcStateSpaceParams,
    rs_param: HcCcCircleParam,
    discretization: f64,
    forward: bool,
}

impl CCpmpmDubinsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64, forward: bool) -> Self {
        Self {
            params: HcCcStateSpaceParams::new(kappa, sigma),
            rs_param: make_rs_circle_param(kappa),
            discretization,
            forward,
        }
    }

    fn ccpmpm_circles_dubins_path(&self, c1: &HcCcCircle, c2: &HcCcCircle) -> CcDubinsPath {
        let mut helper = CcpmpmHelper::new(self.params.hc_cc_circle_param_);
        helper.distance = center_distance(c1, c2);
        helper.angle    = (c2.yc - c1.yc).atan2(c2.xc - c1.xc);

        let mut length = [f64::INFINITY; NB_CC_DUBINS_PATHS];
        let mut qi1 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut qi2 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut qi3 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut qi4 = [None::<Configuration>; NB_CC_DUBINS_PATHS];
        let mut cstart: [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] = std::array::from_fn(|_| None);
        let mut cend:   [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] = std::array::from_fn(|_| None);
        let mut ci1:    [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] = std::array::from_fn(|_| None);
        let mut ci2:    [Option<Box<HcCcCircle>>; NB_CC_DUBINS_PATHS] = std::array::from_fn(|_| None);

        const E:    usize = CcDubinsPathType::E    as usize;
        const T:    usize = CcDubinsPathType::T    as usize;
        const TT:   usize = CcDubinsPathType::TT   as usize;
        const TST:  usize = CcDubinsPathType::TST  as usize;
        const TTT:  usize = CcDubinsPathType::TTT  as usize;
        const TTTT: usize = CcDubinsPathType::TTTT as usize;

        if configuration_equal(&c1.start, &c2.start) {
            length[E] = 0.0;
        } else if configuration_on_hc_cc_circle(c1, &c2.start) {
            let cs = HcCcCircle::from_configuration(
                c1.start, c1.left, c1.forward, true, &self.rs_param,
            );
            length[T] = cs.rs_turn_length(&c2.start);
            cstart[T] = Some(Box::new(cs));
        } else {
            if helper.tt_exists(c1, c2) {
                let (l, cs, ce, q1, q2, q3) = helper.tt_path(c1, c2);
                length[TT] = l;
                qi1[TT] = Some(q1);
                qi2[TT] = Some(q2);
                qi3[TT] = Some(q3);
                cstart[TT] = Some(Box::new(cs));
                cend[TT]   = Some(Box::new(ce));
            }
            if helper.tst_exists(c1, c2) {
                if let Some((l, cs, ce, q1, q2, q3, q4)) = helper.tst_path(c1, c2) {
                    length[TST] = l;
                    qi1[TST] = Some(q1);
                    qi2[TST] = Some(q2);
                    qi3[TST] = Some(q3);
                    qi4[TST] = Some(q4);
                    cstart[TST] = Some(Box::new(cs));
                    cend[TST]   = Some(Box::new(ce));
                }
            }
            if helper.ttt_exists(c1, c2) {
                let (l, cs, ce, q1, q2, q3, mid) = helper.ttt_path(c1, c2);
                length[TTT] = l;
                qi1[TTT] = Some(q1);
                qi2[TTT] = Some(q2);
                qi3[TTT] = Some(q3);
                ci1[TTT] = Some(Box::new(mid));
                cstart[TTT] = Some(Box::new(cs));
                cend[TTT]   = Some(Box::new(ce));
            }
            if helper.tttt_exists(c1, c2) {
                let (l, cs, ce, q1, q2, q3, m1, m2) = helper.tttt_path(c1, c2);
                length[TTTT] = l;
                qi1[TTTT] = Some(q1);
                qi2[TTTT] = Some(q2);
                qi3[TTTT] = Some(q3);
                ci1[TTTT] = Some(Box::new(m1));
                ci2[TTTT] = Some(Box::new(m2));
                cstart[TTTT] = Some(Box::new(cs));
                cend[TTTT]   = Some(Box::new(ce));
            }
        }

        let best = array_index_min(&length);
        CcDubinsPath::new(
            c1.start, c2.start, idx_to_path_type(best),
            self.params.kappa_, self.params.sigma_,
            qi1[best], qi2[best], qi3[best], qi4[best],
            cstart[best].take(), cend[best].take(), ci1[best].take(), ci2[best].take(),
            length[best],
        )
    }

    fn ccpmpm_dubins(&self, state1: &State, state2: &State) -> CcDubinsPath {
        let start1 = Configuration::new(state1.x, state1.y, state1.theta,  self.params.kappa_);
        let start2 = Configuration::new(state1.x, state1.y, state1.theta, -self.params.kappa_);
        let end1   = Configuration::new(state2.x, state2.y, state2.theta,  self.params.kappa_);
        let end2   = Configuration::new(state2.x, state2.y, state2.theta, -self.params.kappa_);
        let rp = &self.rs_param;

        let (start_circles, end_circles) = if self.forward {
            (
                [
                    HcCcCircle::from_configuration(start1, true,  true, true, rp),
                    HcCcCircle::from_configuration(start2, false, true, true, rp),
                ],
                [
                    HcCcCircle::from_configuration(end1, true,  false, true, rp),
                    HcCcCircle::from_configuration(end2, false, false, true, rp),
                ],
            )
        } else {
            (
                [
                    HcCcCircle::from_configuration(start1, true,  false, true, rp),
                    HcCcCircle::from_configuration(start2, false, false, true, rp),
                ],
                [
                    HcCcCircle::from_configuration(end1, true,  true, true, rp),
                    HcCcCircle::from_configuration(end2, false, true, true, rp),
                ],
            )
        };

        let mut best_path: Option<CcDubinsPath> = None;
        let mut best_len = f64::INFINITY;
        for i in 0..2 {
            if i == 0 && state1.kappa < 0.0 { continue; }
            if i == 1 && state1.kappa > 0.0 { continue; }
            for j in 0..2 {
                if j == 0 && state2.kappa < 0.0 { continue; }
                if j == 1 && state2.kappa > 0.0 { continue; }
                let path = self.ccpmpm_circles_dubins_path(&start_circles[i], &end_circles[j]);
                if path.length < best_len {
                    best_len = path.length;
                    best_path = Some(path);
                }
            }
        }
        best_path.unwrap()
    }
}

impl StateSpace for CCpmpmDubinsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        let mut controls = Vec::new();
        let p = self.ccpmpm_dubins(s1, s2);
        match p.path_type {
            CcDubinsPathType::E => empty_controls(&mut controls),
            CcDubinsPathType::T => {
                rs_turn_controls(p.cstart.as_ref().unwrap(), &p.end, true, &mut controls);
            }
            CcDubinsPathType::TT => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi3 = p.qi3.as_ref().unwrap();
                hc_turn_controls(p.cstart.as_ref().unwrap(), qi1, false, &mut controls);
                hc_turn_controls(p.cend.as_ref().unwrap(),   qi3, true,  &mut controls);
            }
            CcDubinsPathType::TST => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                let qi3 = p.qi3.as_ref().unwrap();
                let qi4 = p.qi4.as_ref().unwrap();
                hc_turn_controls(p.cstart.as_ref().unwrap(), qi1, false, &mut controls);
                straight_controls(qi2, qi3, &mut controls);
                hc_turn_controls(p.cend.as_ref().unwrap(),   qi4, true,  &mut controls);
            }
            CcDubinsPathType::TTT => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                let qi3 = p.qi3.as_ref().unwrap();
                hc_turn_controls(p.cstart.as_ref().unwrap(), qi1, false, &mut controls);
                cc_turn_controls(p.ci1.as_ref().unwrap(),    qi2, true,  &mut controls);
                hc_turn_controls(p.cend.as_ref().unwrap(),   qi3, true,  &mut controls);
            }
            CcDubinsPathType::TTTT => {
                let qi1 = p.qi1.as_ref().unwrap();
                let qi2 = p.qi2.as_ref().unwrap();
                let qi3 = p.qi3.as_ref().unwrap();
                hc_turn_controls(p.cstart.as_ref().unwrap(), qi1, false, &mut controls);
                cc_turn_controls(p.ci1.as_ref().unwrap(),    qi2, true,  &mut controls);
                cc_turn_controls(p.ci2.as_ref().unwrap(),    qi2, false, &mut controls);
                hc_turn_controls(p.cend.as_ref().unwrap(),   qi3, true,  &mut controls);
            }
            _ => {}
        }
        controls
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        vec![self.get_controls(s1, s2)]
    }

    fn discretization(&self) -> f64 { self.discretization }
}

// ===========================================================================
// CC_Dubins_State_Space  (general wrapper)
// ===========================================================================

pub struct CCDubinsStateSpace {
    params: HcCcStateSpaceParams,
    cc00:   CC00DubinsStateSpace,
    cc0pm:  CC0pmDubinsStateSpace,
    ccpm0:  CCpm0DubinsStateSpace,
    ccpmpm: CCpmpmDubinsStateSpace,
    discretization: f64,
    forward: bool,
}

impl CCDubinsStateSpace {
    pub fn new(kappa: f64, sigma: f64, discretization: f64, forward: bool) -> Self {
        Self {
            params: HcCcStateSpaceParams::new(kappa, sigma),
            cc00:   CC00DubinsStateSpace::new(kappa, sigma, discretization, forward),
            cc0pm:  CC0pmDubinsStateSpace::new(kappa, sigma, discretization, forward),
            ccpm0:  CCpm0DubinsStateSpace::new(kappa, sigma, discretization, forward),
            ccpmpm: CCpmpmDubinsStateSpace::new(kappa, sigma, discretization, forward),
            discretization,
            forward,
        }
    }

    /// Predict states from `state` ramping curvature to 0 or to ±kappa_max.
    /// Returns a list of `(predicted_state, control_to_reach_it)`.
    fn predict_state(&self, state: &State, forwards: bool) -> Vec<(State, Control)> {
        let mut result = Vec::new();

        // No prediction needed if kappa is already 0 or at max
        if state.kappa.abs() < get_epsilon()
            || (self.params.kappa_ - state.kappa.abs()).abs() < get_epsilon()
        {
            let mut st = State::default();
            st.x     = state.x;
            st.y     = state.y;
            st.theta = state.theta;
            st.kappa = state.kappa;
            st.d     = state.d;
            let ctrl = Control { delta_s: 0.0, kappa: state.kappa, sigma: 0.0 };
            result.push((st, ctrl));
            return result;
        }

        let sgn_kappa = sgn(state.kappa);
        let kappa   = self.params.kappa_;
        let sigma   = self.params.sigma_;

        let (c1, c2) = if forwards {
            (
                Control {
                    delta_s: (kappa - sgn_kappa * state.kappa) / sigma,
                    kappa:   state.kappa,
                    sigma:   sgn_kappa * sigma,
                },
                Control {
                    delta_s: sgn_kappa * state.kappa / sigma,
                    kappa:   state.kappa,
                    sigma:   -sgn_kappa * sigma,
                },
            )
        } else {
            (
                Control {
                    delta_s: -(kappa - sgn_kappa * state.kappa) / sigma,
                    kappa:    state.kappa,
                    sigma:    sgn_kappa * sigma,
                },
                Control {
                    delta_s: -sgn_kappa * state.kappa / sigma,
                    kappa:    state.kappa,
                    sigma:    -sgn_kappa * sigma,
                },
            )
        };

        for ctrl in [c1, c2] {
            let d = sgn(ctrl.delta_s);
            let abs_ds = ctrl.delta_s.abs();
            let (xf, yf, thf, kf) =
                end_of_clothoid(state.x, state.y, state.theta, state.kappa, ctrl.sigma, d, abs_ds);
            let mut predicted = State::default();
            predicted.x     = xf;
            predicted.y     = yf;
            predicted.theta = thf;
            predicted.kappa = kf;
            result.push((predicted, ctrl));
        }
        result
    }

    fn states_equal(s1: &State, s2: &State) -> bool {
        let c1 = Configuration::new(s1.x, s1.y, s1.theta, s1.kappa);
        let c2 = Configuration::new(s2.x, s2.y, s2.theta, s2.kappa);
        state_equal(&c1, &c2)
    }

    #[allow(dead_code)]
    fn get_distance_inner(&self, state1: &State, state2: &State) -> f64 {
        let start_scs = self.predict_state(state1,  self.forward);
        let end_scs   = self.predict_state(state2, !self.forward);
        let mut distances = Vec::new();

        for (start_state, start_ctrl) in &start_scs {
            for (end_state, end_ctrl) in &end_scs {
                if Self::states_equal(start_state, end_state) {
                    let ctrl = subtract_control(start_ctrl, end_ctrl);
                    distances.push(ctrl.delta_s.abs());
                } else {
                    let mut dist = if start_state.kappa.abs() < get_epsilon() {
                        if end_state.kappa.abs() < get_epsilon() {
                            self.cc00.cc00_dubins(start_state, end_state).length
                        } else {
                            self.cc0pm.cc0pm_dubins(start_state, end_state).length
                        }
                    } else if end_state.kappa.abs() < get_epsilon() {
                        self.ccpm0.ccpm0_dubins(start_state, end_state).length
                    } else {
                        self.ccpmpm.ccpmpm_dubins(start_state, end_state).length
                    };
                    if start_ctrl.delta_s.abs() > get_epsilon() {
                        dist += start_ctrl.delta_s.abs();
                    }
                    if end_ctrl.delta_s.abs() > get_epsilon() {
                        dist += end_ctrl.delta_s.abs();
                    }
                    distances.push(dist);
                }
            }
        }

        distances.into_iter().fold(f64::INFINITY, f64::min)
    }

    fn get_controls_inner(&self, state1: &State, state2: &State) -> Vec<Control> {
        let start_scs = self.predict_state(state1,  self.forward);
        let end_scs   = self.predict_state(state2, !self.forward);
        let mut pairs: Vec<(Vec<Control>, f64)> = Vec::new();

        for (start_state, start_ctrl) in &start_scs {
            for (end_state, end_ctrl) in &end_scs {
                let cc_ctrls: Vec<Control> = if Self::states_equal(start_state, end_state) {
                    vec![subtract_control(start_ctrl, end_ctrl)]
                } else {
                    let inner = if start_state.kappa.abs() < get_epsilon() {
                        if end_state.kappa.abs() < get_epsilon() {
                            self.cc00.get_controls(start_state, end_state)
                        } else {
                            self.cc0pm.get_controls(start_state, end_state)
                        }
                    } else if end_state.kappa.abs() < get_epsilon() {
                        self.ccpm0.get_controls(start_state, end_state)
                    } else {
                        self.ccpmpm.get_controls(start_state, end_state)
                    };
                    let mut v = inner;
                    if start_ctrl.delta_s.abs() > get_epsilon() {
                        v.insert(0, *start_ctrl);
                    }
                    if end_ctrl.delta_s.abs() > get_epsilon() {
                        let mut ec = *end_ctrl;
                        reverse_control(&mut ec);
                        v.push(ec);
                    }
                    v
                };
                let total: f64 = cc_ctrls.iter().map(|c| c.delta_s.abs()).sum();
                pairs.push((cc_ctrls, total));
            }
        }

        pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        pairs.into_iter().next().map(|(v, _)| v).unwrap_or_default()
    }

    fn get_controls_reverse(&self, state1: &State, state2: &State) -> Vec<Control> {
        let mut ctrls = self.get_controls_inner(state2, state1);
        ctrls.reverse();
        for c in &mut ctrls {
            reverse_control(c);
        }
        ctrls
    }
}

impl StateSpace for CCDubinsStateSpace {
    fn get_controls(&self, s1: &State, s2: &State) -> Vec<Control> {
        self.get_controls_inner(s1, s2)
    }

    fn get_all_controls(&self, s1: &State, s2: &State) -> Vec<Vec<Control>> {
        vec![
            self.get_controls_inner(s1, s2),
            self.get_controls_reverse(s1, s2),
        ]
    }

    fn discretization(&self) -> f64 { self.discretization }
}

// ---------------------------------------------------------------------------
// Keep the old stub alias so existing code that uses CcDubinsStateSpace still
// compiles while the full implementation is provided above.
// ---------------------------------------------------------------------------
pub use CCDubinsStateSpace as CcDubinsStateSpace;
