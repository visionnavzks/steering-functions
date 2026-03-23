use crate::utilities::{end_of_circular_arc, end_of_straight_line, sgn, twopify, TWO_PI};
use crate::{Control, State, StateSpace};
use std::f64;

const DUBINS_EPS: f64 = 1e-6;
const DUBINS_ZERO: f64 = -1e-9;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DubinsPathSegmentType {
    Left = 0,
    Straight = 1,
    Right = 2,
}

const DUBINS_PATH_TYPE: [[DubinsPathSegmentType; 3]; 6] = [
    [
        DubinsPathSegmentType::Left,
        DubinsPathSegmentType::Straight,
        DubinsPathSegmentType::Left,
    ],
    [
        DubinsPathSegmentType::Right,
        DubinsPathSegmentType::Straight,
        DubinsPathSegmentType::Right,
    ],
    [
        DubinsPathSegmentType::Right,
        DubinsPathSegmentType::Straight,
        DubinsPathSegmentType::Left,
    ],
    [
        DubinsPathSegmentType::Left,
        DubinsPathSegmentType::Straight,
        DubinsPathSegmentType::Right,
    ],
    [
        DubinsPathSegmentType::Right,
        DubinsPathSegmentType::Left,
        DubinsPathSegmentType::Right,
    ],
    [
        DubinsPathSegmentType::Left,
        DubinsPathSegmentType::Right,
        DubinsPathSegmentType::Left,
    ],
];

#[derive(Debug, Clone, Copy)]
pub struct DubinsPath {
    pub types: [DubinsPathSegmentType; 3],
    pub length: [f64; 3],
}

impl DubinsPath {
    pub fn new(types: [DubinsPathSegmentType; 3], t: f64, p: f64, q: f64) -> Self {
        Self {
            types,
            length: [t, p, q],
        }
    }

    pub fn length(&self) -> f64 {
        self.length[0] + self.length[1] + self.length[2]
    }

    pub fn empty() -> Self {
        Self::new(DUBINS_PATH_TYPE[0], f64::MAX, f64::MAX, f64::MAX)
    }
}

pub struct DubinsStateSpace {
    pub kappa: f64,
    pub kappa_inv: f64,
    pub discretization: f64,
    pub forwards: bool,
}

impl DubinsStateSpace {
    pub fn new(kappa: f64, discretization: f64, forwards: bool) -> Self {
        assert!(kappa > 0.0 && discretization > 0.0);
        Self {
            kappa,
            kappa_inv: 1.0 / kappa,
            discretization,
            forwards,
        }
    }

    fn dubins_lsl(d: f64, alpha: f64, beta: f64) -> Option<DubinsPath> {
        let (ca, sa) = (alpha.cos(), alpha.sin());
        let (cb, sb) = (beta.cos(), beta.sin());
        let tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sa - sb));
        if tmp >= DUBINS_ZERO {
            let theta = (cb - ca).atan2(d + sa - sb);
            let t = twopify(-alpha + theta);
            let p = tmp.max(0.0).sqrt();
            let q = twopify(beta - theta);
            return Some(DubinsPath::new(DUBINS_PATH_TYPE[0], t, p, q));
        }
        None
    }

    fn dubins_rsr(d: f64, alpha: f64, beta: f64) -> Option<DubinsPath> {
        let (ca, sa) = (alpha.cos(), alpha.sin());
        let (cb, sb) = (beta.cos(), beta.sin());
        let tmp = 2.0 + d * d - 2.0 * (ca * cb + sa * sb - d * (sb - sa));
        if tmp >= DUBINS_ZERO {
            let theta = (ca - cb).atan2(d - sa + sb);
            let t = twopify(alpha - theta);
            let p = tmp.max(0.0).sqrt();
            let q = twopify(-beta + theta);
            return Some(DubinsPath::new(DUBINS_PATH_TYPE[1], t, p, q));
        }
        None
    }

    fn dubins_rsl(d: f64, alpha: f64, beta: f64) -> Option<DubinsPath> {
        let (ca, sa) = (alpha.cos(), alpha.sin());
        let (cb, sb) = (beta.cos(), beta.sin());
        let tmp = d * d - 2.0 + 2.0 * (ca * cb + sa * sb - d * (sa + sb));
        if tmp >= DUBINS_ZERO {
            let p = tmp.max(0.0).sqrt();
            let theta = (ca + cb).atan2(d - sa - sb) - (2.0_f64).atan2(p);
            let t = twopify(alpha - theta);
            let q = twopify(beta - theta);
            return Some(DubinsPath::new(DUBINS_PATH_TYPE[2], t, p, q));
        }
        None
    }

    fn dubins_lsr(d: f64, alpha: f64, beta: f64) -> Option<DubinsPath> {
        let (ca, sa) = (alpha.cos(), alpha.sin());
        let (cb, sb) = (beta.cos(), beta.sin());
        let tmp = -2.0 + d * d + 2.0 * (ca * cb + sa * sb + d * (sa + sb));
        if tmp >= DUBINS_ZERO {
            let p = tmp.max(0.0).sqrt();
            let theta = (-ca - cb).atan2(d + sa + sb) - (-2.0_f64).atan2(p);
            let t = twopify(-alpha + theta);
            let q = twopify(-beta + theta);
            return Some(DubinsPath::new(DUBINS_PATH_TYPE[3], t, p, q));
        }
        None
    }

    fn dubins_rlr(d: f64, alpha: f64, beta: f64) -> Option<DubinsPath> {
        let (ca, sa) = (alpha.cos(), alpha.sin());
        let (cb, sb) = (beta.cos(), beta.sin());
        let tmp = 0.125 * (6.0 - d * d + 2.0 * (ca * cb + sa * sb + d * (sa - sb)));
        if tmp.abs() < 1.0 {
            let p = TWO_PI - tmp.acos();
            let theta = (ca - cb).atan2(d - sa + sb);
            let t = twopify(alpha - theta + 0.5 * p);
            let q = twopify(alpha - beta - t + p);
            return Some(DubinsPath::new(DUBINS_PATH_TYPE[4], t, p, q));
        }
        None
    }

    fn dubins_lrl(d: f64, alpha: f64, beta: f64) -> Option<DubinsPath> {
        let (ca, sa) = (alpha.cos(), alpha.sin());
        let (cb, sb) = (beta.cos(), beta.sin());
        let tmp = 0.125 * (6.0 - d * d + 2.0 * (ca * cb + sa * sb - d * (sa - sb)));
        if tmp.abs() < 1.0 {
            let p = TWO_PI - tmp.acos();
            let theta = (-ca + cb).atan2(d + sa - sb);
            let t = twopify(-alpha + theta + 0.5 * p);
            let q = twopify(beta - alpha - t + p);
            return Some(DubinsPath::new(DUBINS_PATH_TYPE[5], t, p, q));
        }
        None
    }

    fn dubins_core(d: f64, alpha: f64, beta: f64) -> DubinsPath {
        if d < DUBINS_EPS && (alpha - beta).abs() < DUBINS_EPS {
            return DubinsPath::new(DUBINS_PATH_TYPE[0], 0.0, d, 0.0);
        }

        let mut best_path = DubinsPath::empty();
        let mut min_len = f64::MAX;

        let funcs = [
            Self::dubins_lsl,
            Self::dubins_rsr,
            Self::dubins_rsl,
            Self::dubins_lsr,
            Self::dubins_rlr,
            Self::dubins_lrl,
        ];

        for &func in &funcs {
            if let Some(path) = func(d, alpha, beta) {
                let len = path.length();
                if len < min_len {
                    min_len = len;
                    best_path = path;
                }
            }
        }

        best_path
    }

    pub fn dubins(&self, state1: &State, state2: &State) -> DubinsPath {
        let dx = state2.x - state1.x;
        let dy = state2.y - state1.y;
        let th = dy.atan2(dx);
        let d = (dx * dx + dy * dy).sqrt() * self.kappa;
        let alpha = twopify(state1.theta - th);
        let beta = twopify(state2.theta - th);
        Self::dubins_core(d, alpha, beta)
    }

    pub fn get_distance(&self, state1: &State, state2: &State) -> f64 {
        if self.forwards {
            self.kappa_inv * self.dubins(state1, state2).length()
        } else {
            self.kappa_inv * self.dubins(state2, state1).length()
        }
    }

    fn get_controls_core(
        &self,
        state1: &State,
        state2: &State,
        reverse_direction: bool,
    ) -> Vec<Control> {
        let path = if !reverse_direction {
            self.dubins(state1, state2)
        } else {
            self.dubins(state2, state1)
        };

        let mut controls = Vec::with_capacity(3);

        for i in 0..3 {
            let mut control = Control::default();
            match path.types[i] {
                DubinsPathSegmentType::Left => {
                    control.delta_s = self.kappa_inv * path.length[i];
                    control.kappa = self.kappa;
                    control.sigma = 0.0;
                }
                DubinsPathSegmentType::Straight => {
                    control.delta_s = self.kappa_inv * path.length[i];
                    control.kappa = 0.0;
                    control.sigma = 0.0;
                }
                DubinsPathSegmentType::Right => {
                    control.delta_s = self.kappa_inv * path.length[i];
                    control.kappa = -self.kappa;
                    control.sigma = 0.0;
                }
            }
            if control.delta_s.abs() > DUBINS_EPS {
                controls.push(control);
            }
        }

        if reverse_direction {
            controls.reverse();
            for control in &mut controls {
                control.delta_s = -control.delta_s;
            }
        }

        controls
    }

    pub fn get_controls_reverse(&self, state1: &State, state2: &State) -> Vec<Control> {
        self.get_controls_core(state1, state2, true)
    }

    #[inline]
    pub fn integrate_ode(&self, state: &State, control: &Control, integration_step: f64) -> State {
        let mut state_next = *state;
        let kappa = control.kappa;
        let d = sgn(control.delta_s);

        if kappa.abs() > crate::EPSILON {
            let (x_f, y_f, theta_f) =
                end_of_circular_arc(state.x, state.y, state.theta, kappa, d, integration_step);
            state_next.x = x_f;
            state_next.y = y_f;
            state_next.theta = theta_f;
        } else {
            let (x_f, y_f) =
                end_of_straight_line(state.x, state.y, state.theta, d, integration_step);
            state_next.x = x_f;
            state_next.y = y_f;
        }

        state_next.kappa = kappa;
        state_next.d = d;
        state_next.s = state.s + integration_step;
        state_next
    }
}

impl StateSpace for DubinsStateSpace {
    fn get_controls(&self, state1: &State, state2: &State) -> Vec<Control> {
        self.get_controls_core(state1, state2, false)
    }

    fn get_all_controls(&self, state1: &State, state2: &State) -> Vec<Vec<Control>> {
        vec![
            self.get_controls(state1, state2),
            self.get_controls_reverse(state1, state2),
        ]
    }

    fn get_path(&self, state1: &State, state2: &State) -> Vec<State> {
        let forward_len = self.get_distance(state1, state2);
        let reverse_len = self.kappa_inv * self.dubins(state2, state1).length();

        let controls = if forward_len <= reverse_len {
            self.get_controls_core(state1, state2, false)
        } else {
            self.get_controls_core(state1, state2, true)
        };

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
        path.reserve(n_states + 3);

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
