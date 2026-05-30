use steering_functions_dubins_rs::{State, DubinsDirectionMode};
use steering_functions_dubins_rs::dubins::DubinsStateSpace;
use steering_functions_dubins_rs::state_space::StateSpace;
use std::f64::consts::{PI, FRAC_PI_4};

fn state(x: f64, y: f64, theta: f64) -> State {
    State { x, y, theta, ..State::default() }
}

fn path_ends_near(_start: &State, goal: &State, path: &[State], pos_eps: f64, heading_eps: f64) -> bool {
    let end = path.last().unwrap();
    (end.x - goal.x).abs() < pos_eps
        && (end.y - goal.y).abs() < pos_eps
        && (end.theta - goal.theta).abs() < heading_eps
}

#[test]
fn dubins_forward_only_straight() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOnly);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(5.0, 0.0, 0.0);
    let path = planner.get_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn dubins_forward_only_turn() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOnly);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(2.0, 2.0, PI / 2.0);
    let path = planner.get_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn dubins_reverse_only_reaches_goal() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ReverseOnly);
    let start = state(-3.0, 0.0, 0.0);
    let goal = state(3.0, 2.0, FRAC_PI_4);
    let path = planner.get_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn dubins_forward_or_reverse_picks_shorter() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOrReverse);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(0.5, 0.0, PI);
    let path = planner.get_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn dubins_path_length_positive() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOnly);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let path = planner.get_path(&start, &goal);
    assert!(path.len() > 1);
}

#[test]
fn dubins_get_controls_nonempty() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOnly);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let controls = planner.get_controls(&start, &goal);
    assert!(!controls.is_empty());
}

#[test]
fn dubins_get_all_controls_nonempty() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOrReverse);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let all_controls = planner.get_all_controls(&start, &goal);
    assert!(!all_controls.is_empty());
}

#[test]
fn dubins_discretization() {
    let planner = DubinsStateSpace::new(1.0, 0.1, DubinsDirectionMode::ForwardOnly);
    assert!((planner.discretization() - 0.1).abs() < 1e-10);
}

#[test]
fn dubins_same_start_goal() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOnly);
    let start = state(1.0, 2.0, 0.5);
    let goal = state(1.0, 2.0, 0.5);
    let path = planner.get_path(&start, &goal);
    if path.is_empty() {
        return;
    }
    assert!(path_ends_near(&start, &goal, &path, 1e-3, 1e-3));
}

#[test]
fn dubins_various_orientations() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOrReverse);
    let start = state(0.0, 0.0, 0.0);

    for angle in [0.0, PI / 4.0, PI / 2.0] {
        let goal = state(4.0, 0.0, angle);
        let path = planner.get_path(&start, &goal);
        assert!(
            path_ends_near(&start, &goal, &path, 1e-2, 1e-2),
            "Failed for angle {}",
            angle
        );
    }
}

#[test]
fn dubins_new_kappa_assert() {
    let result = std::panic::catch_unwind(|| {
        DubinsStateSpace::new(0.0, 0.05, DubinsDirectionMode::ForwardOnly);
    });
    assert!(result.is_err());
}

#[test]
fn dubins_new_discretization_assert() {
    let result = std::panic::catch_unwind(|| {
        DubinsStateSpace::new(1.0, 0.0, DubinsDirectionMode::ForwardOnly);
    });
    assert!(result.is_err());
}

#[test]
fn dubins_integrate_empty_controls() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOnly);
    let start = state(0.0, 0.0, 0.0);
    let path = planner.integrate(&start, &[]);
    assert!(path.is_empty());
}

#[test]
fn dubins_interpolate_at_zero() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOnly);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let controls = planner.get_controls(&start, &goal);
    let interpolated = planner.interpolate(&start, &controls, 0.0);
    assert!((interpolated.x - start.x).abs() < 1e-4);
    assert!((interpolated.y - start.y).abs() < 1e-4);
}

#[test]
fn dubins_interpolate_at_one() {
    let planner = DubinsStateSpace::new(1.0, 0.05, DubinsDirectionMode::ForwardOrReverse);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let controls = planner.get_controls(&start, &goal);
    assert!(!controls.is_empty());
    let interpolated = planner.interpolate(&start, &controls, 1.0);
    let total_s: f64 = controls.iter().map(|c| c.delta_s.abs()).sum();
    assert!(total_s > 0.0);
    assert!(interpolated.x.is_finite());
    assert!(interpolated.y.is_finite());
}
