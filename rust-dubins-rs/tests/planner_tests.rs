use steering_functions_dubins_rs::{State, PathType, SteeringPath, DubinsDirectionMode};
use std::f64::consts::PI;

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
fn planner_new_dubins() {
    let p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    assert_eq!(p.path_type(), PathType::Dubins);
    assert!((p.kappa_max() - 1.0).abs() < 1e-10);
    assert!((p.discretization() - 0.05).abs() < 1e-10);
}

#[test]
fn planner_new_rs() {
    let p = SteeringPath::new(PathType::Rs, 1.0, 0.05);
    assert_eq!(p.path_type(), PathType::Rs);
}

#[test]
fn planner_try_new_valid() {
    let p = SteeringPath::try_new(PathType::Dubins, 1.0, 0.05);
    assert!(p.is_ok());
}

#[test]
fn planner_try_new_invalid_kappa() {
    let p = SteeringPath::try_new(PathType::Dubins, -1.0, 0.05);
    assert!(p.is_err());
}

#[test]
fn planner_try_new_invalid_disc() {
    let p = SteeringPath::try_new(PathType::Dubins, 1.0, 0.0);
    assert!(p.is_err());
}

#[test]
fn planner_set_kappa_max_valid() {
    let mut p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    assert!(p.set_kappa_max(2.0).is_ok());
    assert!((p.kappa_max() - 2.0).abs() < 1e-10);
}

#[test]
fn planner_set_kappa_max_invalid() {
    let mut p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    assert!(p.set_kappa_max(-1.0).is_err());
    assert!((p.kappa_max() - 1.0).abs() < 1e-10);
}

#[test]
fn planner_set_discretization_valid() {
    let mut p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    assert!(p.set_discretization(0.1).is_ok());
    assert!((p.discretization() - 0.1).abs() < 1e-10);
}

#[test]
fn planner_set_discretization_invalid() {
    let mut p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    assert!(p.set_discretization(-1.0).is_err());
    assert!((p.discretization() - 0.05).abs() < 1e-10);
}

#[test]
fn planner_set_path_type() {
    let mut p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    p.set_path_type(PathType::Rs);
    assert_eq!(p.path_type(), PathType::Rs);
}

#[test]
fn planner_supported_path_types() {
    let types = SteeringPath::supported_path_types();
    assert!(types.contains(&PathType::Dubins));
    assert!(types.contains(&PathType::Rs));
}

#[test]
fn planner_with_dubins_direction_mode() {
    let p = SteeringPath::new(PathType::Dubins, 1.0, 0.05)
        .with_dubins_direction_mode(DubinsDirectionMode::ReverseOnly);
    assert_eq!(p.dubins_direction_mode(), DubinsDirectionMode::ReverseOnly);
}

#[test]
fn planner_set_dubins_direction_mode() {
    let mut p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    p.set_dubins_direction_mode(DubinsDirectionMode::ForwardOrReverse);
    assert_eq!(p.dubins_direction_mode(), DubinsDirectionMode::ForwardOrReverse);
}

#[test]
fn planner_dubins_compute_path() {
    let p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let path = p.compute_shortest_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn planner_rs_compute_path() {
    let p = SteeringPath::new(PathType::Rs, 1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let path = p.compute_shortest_path(&start, &goal);
    assert!(path_ends_near(&start, &goal, &path, 1e-4, 1e-4));
}

#[test]
fn planner_dubins_compute_controls() {
    let p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let controls = p.compute_shortest_control_sequence(&start, &goal);
    assert!(!controls.is_empty());
}

#[test]
fn planner_rs_compute_controls() {
    let p = SteeringPath::new(PathType::Rs, 1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let controls = p.compute_shortest_control_sequence(&start, &goal);
    assert!(!controls.is_empty());
}

#[test]
fn planner_dubins_compute_all_paths() {
    let p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let all_paths = p.compute_all_paths(&start, &goal);
    assert!(!all_paths.is_empty());
}

#[test]
fn planner_rs_compute_all_paths() {
    let p = SteeringPath::new(PathType::Rs, 1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let all_paths = p.compute_all_paths(&start, &goal);
    assert!(!all_paths.is_empty());
}

#[test]
fn planner_dubins_compute_all_controls() {
    let p = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let all_controls = p.compute_all_control_sequences(&start, &goal);
    assert!(!all_controls.is_empty());
}

#[test]
fn planner_rs_compute_all_controls() {
    let p = SteeringPath::new(PathType::Rs, 1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);
    let goal = state(3.0, 3.0, PI / 4.0);
    let all_controls = p.compute_all_control_sequences(&start, &goal);
    assert!(!all_controls.is_empty());
}

#[test]
fn planner_new_kappa_assert() {
    let result = std::panic::catch_unwind(|| {
        SteeringPath::new(PathType::Dubins, 0.0, 0.05);
    });
    assert!(result.is_err());
}

#[test]
fn planner_new_disc_assert() {
    let result = std::panic::catch_unwind(|| {
        SteeringPath::new(PathType::Dubins, 1.0, 0.0);
    });
    assert!(result.is_err());
}

fn planner_various_orientations(path_type: PathType) {
    let p = SteeringPath::new(path_type, 1.0, 0.05);
    let start = state(0.0, 0.0, 0.0);

    for angle in [0.0, PI / 4.0, PI / 2.0] {
        let goal = state(4.0, 0.0, angle);
        let path = p.compute_shortest_path(&start, &goal);
        if path.is_empty() {
            continue;
        }
        assert!(
            path_ends_near(&start, &goal, &path, 1e-2, 1e-2),
            "{:?} failed for angle {}",
            path_type,
            angle
        );
    }
}

#[test]
fn planner_various_orientations_dubins() {
    planner_various_orientations(PathType::Dubins);
}

#[test]
fn planner_various_orientations_rs() {
    planner_various_orientations(PathType::Rs);
}
