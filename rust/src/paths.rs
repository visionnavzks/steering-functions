use crate::configuration::Configuration;
use crate::hc_cc_circle::HcCcCircle;
use crate::state::Control;
use crate::utilities::{get_epsilon, point_distance, twopify, sgn};

// ---------------------------------------------------------------------------
// Path-type enums
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum CcDubinsPathType {
    E = 0,
    S,
    T,
    TT,
    TST,
    TTT,
    TTTT,
}
pub const NB_CC_DUBINS_PATHS: usize = 7;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum HcCcRsPathType {
    E = 0,
    S,
    T,
    TT,
    TcT,
    TcTcT,
    TcTT,
    TTcT,
    TST,
    TSTcT,
    TcTST,
    TcTSTcT,
    TTcTT,
    TcTTcT,
    TTT,
    TcST,
    TScT,
    TcScT,
}
pub const NB_HC_CC_RS_PATHS: usize = 18;

// ---------------------------------------------------------------------------
// Base path struct
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct Path {
    pub start: Configuration,
    pub end: Configuration,
    pub kappa: f64,
    pub sigma: f64,
    pub length: f64,
}

impl Path {
    pub fn new(
        start: Configuration,
        end: Configuration,
        kappa: f64,
        sigma: f64,
        length: f64,
    ) -> Self {
        Self { start, end, kappa, sigma, length }
    }
}

// ---------------------------------------------------------------------------
// CC-Dubins path
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct CcDubinsPath {
    pub start: Configuration,
    pub end: Configuration,
    pub kappa: f64,
    pub sigma: f64,
    pub length: f64,
    pub path_type: CcDubinsPathType,
    pub qi1: Option<Configuration>,
    pub qi2: Option<Configuration>,
    pub qi3: Option<Configuration>,
    pub qi4: Option<Configuration>,
    pub cstart: Option<Box<HcCcCircle>>,
    pub cend: Option<Box<HcCcCircle>>,
    pub ci1: Option<Box<HcCcCircle>>,
    pub ci2: Option<Box<HcCcCircle>>,
}

impl CcDubinsPath {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        start: Configuration,
        end: Configuration,
        path_type: CcDubinsPathType,
        kappa: f64,
        sigma: f64,
        qi1: Option<Configuration>,
        qi2: Option<Configuration>,
        qi3: Option<Configuration>,
        qi4: Option<Configuration>,
        cstart: Option<Box<HcCcCircle>>,
        cend: Option<Box<HcCcCircle>>,
        ci1: Option<Box<HcCcCircle>>,
        ci2: Option<Box<HcCcCircle>>,
        length: f64,
    ) -> Self {
        Self {
            start,
            end,
            kappa,
            sigma,
            length,
            path_type,
            qi1,
            qi2,
            qi3,
            qi4,
            cstart,
            cend,
            ci1,
            ci2,
        }
    }
}

// ---------------------------------------------------------------------------
// HC/CC Reeds-Shepp path
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct HcCcRsPath {
    pub start: Configuration,
    pub end: Configuration,
    pub kappa: f64,
    pub sigma: f64,
    pub length: f64,
    pub path_type: HcCcRsPathType,
    pub qi1: Option<Configuration>,
    pub qi2: Option<Configuration>,
    pub qi3: Option<Configuration>,
    pub qi4: Option<Configuration>,
    pub cstart: Option<Box<HcCcCircle>>,
    pub cend: Option<Box<HcCcCircle>>,
    pub ci1: Option<Box<HcCcCircle>>,
    pub ci2: Option<Box<HcCcCircle>>,
}

impl HcCcRsPath {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        start: Configuration,
        end: Configuration,
        path_type: HcCcRsPathType,
        kappa: f64,
        sigma: f64,
        qi1: Option<Configuration>,
        qi2: Option<Configuration>,
        qi3: Option<Configuration>,
        qi4: Option<Configuration>,
        cstart: Option<Box<HcCcCircle>>,
        cend: Option<Box<HcCcCircle>>,
        ci1: Option<Box<HcCcCircle>>,
        ci2: Option<Box<HcCcCircle>>,
        length: f64,
    ) -> Self {
        Self {
            start,
            end,
            kappa,
            sigma,
            length,
            path_type,
            qi1,
            qi2,
            qi3,
            qi4,
            cstart,
            cend,
            ci1,
            ci2,
        }
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

pub fn state_equal(s1: &Configuration, s2: &Configuration) -> bool {
    if (s2.kappa - s1.kappa).abs() > get_epsilon() {
        return false;
    }
    if (twopify(s2.theta) - twopify(s1.theta)).abs() > get_epsilon() {
        return false;
    }
    point_distance(s1.x, s1.y, s2.x, s2.y) <= get_epsilon()
}

pub fn reverse_control(control: &mut Control) {
    control.kappa += control.delta_s.abs() * control.sigma;
    control.delta_s = -control.delta_s;
    control.sigma = -control.sigma;
}

pub fn subtract_control(c1: &Control, c2: &Control) -> Control {
    debug_assert!(
        (sgn(c1.delta_s) * c1.sigma - sgn(c2.delta_s) * c2.sigma).abs() < get_epsilon()
    );
    Control {
        delta_s: c1.delta_s - c2.delta_s,
        kappa: c1.kappa,
        sigma: c1.sigma,
    }
}

pub fn empty_controls(controls: &mut Vec<Control>) {
    controls.push(Control { delta_s: 0.0, kappa: 0.0, sigma: 0.0 });
}

pub fn straight_controls(q1: &Configuration, q2: &Configuration, controls: &mut Vec<Control>) {
    let length = point_distance(q1.x, q1.y, q2.x, q2.y);
    let dot = q1.theta.cos() * (q2.x - q1.x) + q1.theta.sin() * (q2.y - q1.y);
    let d = sgn(dot);
    controls.push(Control { delta_s: d * length, kappa: 0.0, sigma: 0.0 });
}

fn direction(forward: bool, order: bool) -> f64 {
    if (forward && order) || (!forward && !order) { 1.0 } else { -1.0 }
}

pub fn rs_turn_controls(c: &HcCcCircle, q: &Configuration, order: bool, controls: &mut Vec<Control>) {
    let delta = c.deflection(q);
    let length_arc = c.kappa_inv.abs() * c.rs_circular_deflection(delta);
    let d = direction(c.forward, order);
    controls.push(Control { delta_s: d * length_arc, kappa: c.kappa, sigma: 0.0 });
}

pub fn hc_turn_controls(c: &HcCcCircle, q: &Configuration, order: bool, controls: &mut Vec<Control>) {
    let delta = c.deflection(q);
    let length_min = (c.kappa / c.sigma).abs();
    let length_arc = c.kappa_inv.abs() * c.hc_circular_deflection(delta);
    let d = direction(c.forward, order);

    if order {
        controls.push(Control { delta_s: d * length_min, kappa: 0.0, sigma: c.sigma });
    }
    controls.push(Control { delta_s: d * length_arc, kappa: c.kappa, sigma: 0.0 });
    if !order {
        controls.push(Control { delta_s: d * length_min, kappa: c.kappa, sigma: -c.sigma });
    }
}

pub fn cc_elementary_controls(
    c: &HcCcCircle,
    q: &Configuration,
    delta: f64,
    order: bool,
    controls: &mut Vec<Control>,
) -> bool {
    let (success, sigma0) = c.cc_elementary_sharpness(q, delta);
    if success {
        let length = (delta / sigma0.abs()).sqrt();
        let d = direction(c.forward, order);
        controls.push(Control { delta_s: d * length, kappa: 0.0, sigma: sigma0 });
        controls.push(Control { delta_s: d * length, kappa: sigma0 * length, sigma: -sigma0 });
    }
    success
}

pub fn cc_default_controls(
    c: &HcCcCircle,
    q: &Configuration,
    delta: f64,
    order: bool,
    controls: &mut Vec<Control>,
) {
    let length_min = (c.kappa / c.sigma).abs();
    let length_arc = c.kappa_inv.abs() * c.cc_circular_deflection(delta);
    let d = direction(c.forward, order);
    controls.push(Control { delta_s: d * length_min, kappa: 0.0, sigma: c.sigma });
    controls.push(Control { delta_s: d * length_arc, kappa: c.kappa, sigma: 0.0 });
    controls.push(Control { delta_s: d * length_min, kappa: c.kappa, sigma: -c.sigma });
}

pub fn cc_turn_controls(c: &HcCcCircle, q: &Configuration, order: bool, controls: &mut Vec<Control>) {
    let delta = c.deflection(q);
    if delta < get_epsilon() {
        if order {
            straight_controls(&c.start, q, controls);
        } else {
            straight_controls(q, &c.start, controls);
        }
        return;
    }
    if delta < 2.0 * c.delta_min {
        let mut controls_elem: Vec<Control> = Vec::new();
        if cc_elementary_controls(c, q, delta, order, &mut controls_elem) {
            let mut controls_default: Vec<Control> = Vec::new();
            cc_default_controls(c, q, delta, order, &mut controls_default);
            let len_elem: f64 = controls_elem.iter().map(|x| x.delta_s.abs()).sum();
            let len_default: f64 = controls_default.iter().map(|x| x.delta_s.abs()).sum();
            if len_elem < len_default {
                controls.extend(controls_elem);
            } else {
                controls.extend(controls_default);
            }
        } else {
            cc_default_controls(c, q, delta, order, controls);
        }
    } else {
        cc_default_controls(c, q, delta, order, controls);
    }
}
