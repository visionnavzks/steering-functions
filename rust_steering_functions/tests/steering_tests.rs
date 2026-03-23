use rust_steering_functions::*;
use std::f64::consts::PI;

#[test]
fn test_dubins_straight_line() {
    let state1 = State::new(0.0, 0.0, 0.0, 0.0);
    let state2 = State::new(10.0, 0.0, 0.0, 0.0);

    let path = SteeringFunctions::get_dubins_path(&state1, &state2, 1.0, 0.1, true);

    assert!(!path.is_empty());

    let end_state = path.last().unwrap();
    assert!((end_state.x - 10.0).abs() < 1e-4);
    assert!((end_state.y - 0.0).abs() < 1e-4);
    assert!((end_state.theta - 0.0).abs() < 1e-4);
}

#[test]
fn test_reeds_shepp_reverse() {
    let state1 = State::new(0.0, 0.0, 0.0, 0.0);
    let state2 = State::new(-10.0, 0.0, 0.0, 0.0);

    let path = SteeringFunctions::get_reeds_shepp_path(&state1, &state2, 1.0, 0.1);

    assert!(!path.is_empty());

    let end_state = path.last().unwrap();
    assert!((end_state.x - (-10.0)).abs() < 1e-4);
    assert!((end_state.y - 0.0).abs() < 1e-4);
    assert!((end_state.theta - 0.0).abs() < 1e-4);
}

#[test]
fn test_dubins_curve() {
    let state1 = State::new(0.0, 0.0, 0.0, 0.0);
    let state2 = State::new(0.0, 10.0, PI, 0.0);

    let path = SteeringFunctions::get_dubins_path(&state1, &state2, 0.2, 0.1, true);

    assert!(!path.is_empty());

    let end_state = path.last().unwrap();
    assert!((end_state.x - 0.0).abs() < 1e-4);
    assert!((end_state.y - 10.0).abs() < 1e-4);
    // Note: angles might wrap
    let theta_diff = (end_state.theta - PI).rem_euclid(2.0 * PI);
    assert!(theta_diff < 1e-4 || (theta_diff - 2.0 * PI).abs() < 1e-4);
}
