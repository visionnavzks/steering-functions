use steering_functions::base_state_space::StateSpace;
use steering_functions::hc_reeds_shepp::{
    Hc0pmRsStateSpace, HcRsStateSpace, Hcpm0RsStateSpace, HcpmpmRsStateSpace,
};
use steering_functions::paths::reverse_control;
use steering_functions::state::{Control, State};

const KAPPA: f64 = 1.0;
const SIGMA: f64 = 1.0;
const DISC: f64 = 0.05;
const EPS: f64 = 1e-6;

fn approx_eq(left: f64, right: f64, eps: f64) {
    assert!((left - right).abs() <= eps, "left={left}, right={right}, eps={eps}");
}

fn assert_state_close(actual: &State, expected: &State, eps_xy: f64, eps_theta: f64) {
    approx_eq(actual.x, expected.x, eps_xy);
    approx_eq(actual.y, expected.y, eps_xy);
    approx_eq(actual.theta, expected.theta, eps_theta);
    approx_eq(actual.kappa, expected.kappa, EPS);
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

fn assert_controls_match(actual: &[Control], expected: &[Control], eps: f64) {
    assert_eq!(actual.len(), expected.len());
    for (actual_control, expected_control) in actual.iter().zip(expected.iter()) {
        approx_eq(actual_control.delta_s, expected_control.delta_s, eps);
        approx_eq(actual_control.kappa, expected_control.kappa, eps);
        approx_eq(actual_control.sigma, expected_control.sigma, eps);
    }
}

#[test]
fn hcpm0_matches_reversed_hc0pm() {
    let hc0pm = Hc0pmRsStateSpace::new(KAPPA, SIGMA, DISC);
    let hcpm0 = Hcpm0RsStateSpace::new(KAPPA, SIGMA, DISC);

    let start = State {
        x: -1.4,
        y: 0.8,
        theta: 0.35,
        kappa: 1.0,
        ..State::default()
    };
    let goal = State {
        x: 1.9,
        y: -0.6,
        theta: -0.9,
        kappa: 0.0,
        ..State::default()
    };

    let controls = hcpm0.get_controls(&start, &goal);
    let mirrored = reverse_controls(&hc0pm.get_controls(&goal, &start));
    assert_controls_match(&controls, &mirrored, EPS);

    let path = hcpm0.get_path(&start, &goal);
    assert_state_close(path.last().unwrap(), &goal, 1e-2, 1e-2);
}

#[test]
fn hcpmpm_regression_sample_matches_expected_length_and_endpoint() {
    let hcpmpm = HcpmpmRsStateSpace::new(KAPPA, SIGMA, DISC);

    let start = State {
        x: -1.2,
        y: 0.7,
        theta: 0.4,
        kappa: 1.0,
        ..State::default()
    };
    let goal = State {
        x: 2.1,
        y: -0.4,
        theta: -0.8,
        kappa: -1.0,
        ..State::default()
    };

    let controls = hcpmpm.get_controls(&start, &goal);
    let path = hcpmpm.get_path(&start, &goal);

    assert!(controls.len() >= 3);
    approx_eq(controls_length(&controls), 4.675812123995, 1e-6);
    assert_state_close(path.last().unwrap(), &goal, 1e-2, 1e-2);
}

#[test]
fn hc_rs_delegates_to_hcpmpm_when_both_states_have_max_curvature() {
    let hc_rs = HcRsStateSpace::new(KAPPA, SIGMA, DISC);
    let hcpmpm = HcpmpmRsStateSpace::new(KAPPA, SIGMA, DISC);

    let start = State {
        x: -1.2,
        y: 0.7,
        theta: 0.4,
        kappa: 1.0,
        ..State::default()
    };
    let goal = State {
        x: 2.1,
        y: -0.4,
        theta: -0.8,
        kappa: -1.0,
        ..State::default()
    };

    let hc_rs_controls = hc_rs.get_controls(&start, &goal);
    let hcpmpm_controls = hcpmpm.get_controls(&start, &goal);
    assert_controls_match(&hc_rs_controls, &hcpmpm_controls, EPS);
}

#[test]
fn hc_rs_reaches_goal_for_mixed_curvature_case() {
    let hc_rs = HcRsStateSpace::new(KAPPA, SIGMA, DISC);

    let start = State {
        x: -2.0,
        y: -0.5,
        theta: 0.2,
        kappa: 0.35,
        ..State::default()
    };
    let goal = State {
        x: 1.7,
        y: 1.4,
        theta: 1.1,
        kappa: -0.4,
        ..State::default()
    };

    let controls = hc_rs.get_controls(&start, &goal);
    let path = hc_rs.get_path(&start, &goal);

    assert!(!controls.is_empty());
    assert!(controls_length(&controls) > 0.0);
    assert_state_close(path.last().unwrap(), &goal, 1e-2, 1e-2);
}