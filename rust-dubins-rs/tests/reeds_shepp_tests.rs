use steering_functions_dubins_rs::{State, reeds_shepp::ReedsSheppStateSpace};
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
fn rs_straight_forward() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(5.0, 0.0, 0.0);
    let path = planner.get_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn rs_straight_reverse() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(-5.0, 0.0, 0.0);
    let path = planner.get_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn rs_turn_left() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(0.0, 2.0, PI / 2.0);
    let path = planner.get_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn rs_turn_right() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(0.0, -2.0, -PI / 2.0);
    let path = planner.get_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn rs_reverse_with_turn() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(-2.0, 2.0, PI);
    let path = planner.get_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn rs_complex_maneuver() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 2.0, FRAC_PI_4);
    let path = planner.get_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn rs_get_controls_nonempty() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let controls = planner.get_controls(&start, &goal);
    assert!(!controls.is_empty());
}

#[test]
fn rs_get_all_controls_nonempty() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let all_controls = planner.get_all_controls(&start, &goal);
    assert!(!all_controls.is_empty());
}

#[test]
fn rs_get_distance_positive() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let dist = planner.get_distance(&start, &goal);
    assert!(dist > 0.0);
}

#[test]
fn rs_discretization() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.1);
    assert!((planner.discretization() - 0.1).abs() < 1e-10);
}

#[test]
fn rs_new_kappa_assert() {
    let result = std::panic::catch_unwind(|| {
        ReedsSheppStateSpace::new(0.0, 0.05);
    });
    assert!(result.is_err());
}

#[test]
fn rs_new_discretization_assert() {
    let result = std::panic::catch_unwind(|| {
        ReedsSheppStateSpace::new(1.0, 0.0);
    });
    assert!(result.is_err());
}

#[test]
fn rs_integrate_empty_controls() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let path = planner.integrate(&start, &[]);
    assert!(path.is_empty());
}

#[test]
fn rs_interpolate_at_zero() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let controls = planner.get_controls(&start, &goal);
    let interpolated = planner.interpolate(&start, &controls, 0.0);
    assert!((interpolated.x - start.x).abs() < 1e-4);
    assert!((interpolated.y - start.y).abs() < 1e-4);
}

#[test]
fn rs_interpolate_at_one() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let controls = planner.get_controls(&start, &goal);
    assert!(!controls.is_empty());
    let interpolated = planner.interpolate(&start, &controls, 1.0);
    assert!(interpolated.x.is_finite());
    assert!(interpolated.y.is_finite());
}

#[test]
fn rs_same_start_goal() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(1.0, 2.0, 0.5);
    let goal = state(1.0, 2.0, 0.5);
    let path = planner.get_path(&start, &goal);
    if path.is_empty() {
        return;
    }
    assert!(path_ends_near(&start, &goal, &path, 1e-3, 1e-3));
}

#[test]
fn rs_various_orientations() {
    let planner = ReedsSheppStateSpace::new(1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);

    for angle in [0.0, PI / 4.0, PI / 2.0] {
        let goal = state(4.0, 0.0, angle);
        let path = planner.get_path(&start, &goal);
        if path.is_empty() {
            continue;
        }
        assert!(
            path_ends_near(&start, &goal, &path, 1e-2, 1e-2),
            "Failed for angle {}",
            angle
        );
    }
}
