use steering_functions_dubins_rs::{State, Control};

#[test]
fn test_state_default() {
    let s = State::default();
    assert_eq!(s.x, 0.0);
    assert_eq!(s.y, 0.0);
    assert_eq!(s.theta, 0.0);
    assert_eq!(s.kappa, 0.0);
    assert_eq!(s.sigma, 0.0);
}

#[test]
fn test_state_nearly_equal_same() {
    let s1 = State { x: 1.0, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    let s2 = State { x: 1.0, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    assert!(s1.nearly_equal(&s2));
}

#[test]
fn test_state_nearly_equal_small_diff() {
    let s1 = State { x: 1.0, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    let s2 = State { x: 1.0 + 1e-7, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    assert!(s1.nearly_equal(&s2));
}

#[test]
fn test_state_nearly_equal_large_diff() {
    let s1 = State { x: 1.0, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    let s2 = State { x: 1.1, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    assert!(!s1.nearly_equal(&s2));
}

#[test]
fn test_state_eq() {
    let s1 = State { x: 1.0, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    let s2 = State { x: 1.0, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    assert_eq!(s1, s2);
}

#[test]
fn test_state_eq_with_tolerance() {
    let s1 = State { x: 1.0, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    let s2 = State { x: 1.0 + 1e-7, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    assert_eq!(s1, s2);
}

#[test]
fn test_state_clone() {
    let s1 = State { x: 1.0, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    let s2 = s1.clone();
    assert_eq!(s1, s2);
}

#[test]
fn test_state_debug() {
    let s = State { x: 1.0, y: 2.0, theta: 3.0, kappa: 4.0, sigma: 5.0 };
    let debug_str = format!("{:?}", s);
    assert!(debug_str.contains("1.0"));
    assert!(debug_str.contains("2.0"));
}

#[test]
fn test_control_default() {
    let c = Control::default();
    assert_eq!(c.delta_s, 0.0);
    assert_eq!(c.kappa, 0.0);
    assert_eq!(c.sigma, 0.0);
}

#[test]
fn test_control_new() {
    let c = Control { delta_s: 1.0, kappa: 2.0, sigma: 3.0 };
    assert_eq!(c.delta_s, 1.0);
    assert_eq!(c.kappa, 2.0);
    assert_eq!(c.sigma, 3.0);
}

#[test]
fn test_control_clone() {
    let c1 = Control { delta_s: 1.0, kappa: 2.0, sigma: 3.0 };
    let c2 = c1.clone();
    assert_eq!(c1.delta_s, c2.delta_s);
    assert_eq!(c1.kappa, c2.kappa);
    assert_eq!(c1.sigma, c2.sigma);
}

#[test]
fn test_control_debug() {
    let c = Control { delta_s: 1.0, kappa: 2.0, sigma: 3.0 };
    let debug_str = format!("{:?}", c);
    assert!(debug_str.contains("1.0"));
    assert!(debug_str.contains("2.0"));
    assert!(debug_str.contains("3.0"));
}
