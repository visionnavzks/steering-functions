use steering_functions_dubins_rs::math::*;
use std::f64::consts::PI;

#[test]
fn test_sgn_positive() {
    assert_eq!(sgn(5.0), 1.0);
}

#[test]
fn test_sgn_negative() {
    assert_eq!(sgn(-3.0), -1.0);
}

#[test]
fn test_sgn_zero() {
    assert_eq!(sgn(0.0), 1.0);
}

#[test]
fn test_polar_origin() {
    let (r, theta) = polar(0.0, 0.0);
    assert!((r - 0.0).abs() < 1e-10);
    assert!((theta - 0.0).abs() < 1e-10);
}

#[test]
fn test_polar_positive_x() {
    let (r, theta) = polar(3.0, 0.0);
    assert!((r - 3.0).abs() < 1e-10);
    assert!((theta - 0.0).abs() < 1e-10);
}

#[test]
fn test_polar_positive_y() {
    let (r, theta) = polar(0.0, 4.0);
    assert!((r - 4.0).abs() < 1e-10);
    assert!((theta - PI / 2.0).abs() < 1e-10);
}

#[test]
fn test_polar_quadrant1() {
    let (r, theta) = polar(3.0, 4.0);
    assert!((r - 5.0).abs() < 1e-10);
    assert!((theta - (4.0_f64).atan2(3.0)).abs() < 1e-10);
}

#[test]
fn test_twopify_zero() {
    assert!((twopify(0.0) - 0.0).abs() < 1e-10);
}

#[test]
fn test_twopify_positive() {
    let result = twopify(PI);
    assert!((result - PI).abs() < 1e-10);
}

#[test]
fn test_twopify_large() {
    let result = twopify(5.0 * PI);
    assert!(result >= 0.0 && result < TWO_PI);
}

#[test]
fn test_twopify_negative() {
    let result = twopify(-PI);
    assert!(result >= 0.0 && result < TWO_PI);
}

#[test]
fn test_pify_zero() {
    assert!((pify(0.0) - 0.0).abs() < 1e-10);
}

#[test]
fn test_pify_positive_pi() {
    let result = pify(PI);
    assert!((result - PI).abs() < 1e-10);
}

#[test]
fn test_pify_negative_pi() {
    let result = pify(-PI);
    assert!((result - PI).abs() < 1e-10 || (result + PI).abs() < 1e-10);
}

#[test]
fn test_pify_beyond_pi() {
    let result = pify(3.0 * PI);
    assert!(result > -PI && result <= PI);
}

#[test]
fn test_pify_negative_beyond() {
    let result = pify(-3.0 * PI);
    assert!((result - PI).abs() < 1e-10 || (result + PI).abs() < 1e-10);
}

#[test]
fn test_fresnel_zero() {
    let (s, c) = fresnel(0.0);
    assert!(s.abs() < 1e-10);
    assert!(c.abs() < 1e-10);
}

#[test]
fn test_fresnel_positive() {
    let (s, c) = fresnel(1.0);
    assert!(s > 0.0 && s < 1.0);
    assert!(c > 0.0 && c < 1.0);
}

#[test]
fn test_fresnel_negative() {
    let (s_pos, c_pos) = fresnel(1.0);
    let (s_neg, c_neg) = fresnel(-1.0);
    assert!((s_pos + s_neg).abs() < 1e-10);
    assert!((c_pos + c_neg).abs() < 1e-10);
}

#[test]
fn test_fresnel_large() {
    let (s, c) = fresnel(10.0);
    assert!((s - 0.5).abs() < 0.1);
    assert!((c - 0.5).abs() < 0.1);
}

#[test]
fn test_end_of_straight_line_forward() {
    let (x_f, y_f) = end_of_straight_line(0.0, 0.0, 0.0, 1.0, 5.0);
    assert!((x_f - 5.0).abs() < 1e-10);
    assert!((y_f - 0.0).abs() < 1e-10);
}

#[test]
fn test_end_of_straight_line_backward() {
    let (x_f, y_f) = end_of_straight_line(5.0, 0.0, 0.0, -1.0, 5.0);
    assert!((x_f - 0.0).abs() < 1e-10);
    assert!((y_f - 0.0).abs() < 1e-10);
}

#[test]
fn test_end_of_straight_line_angle() {
    let (x_f, y_f) = end_of_straight_line(0.0, 0.0, PI / 2.0, 1.0, 3.0);
    assert!(x_f.abs() < 1e-10);
    assert!((y_f - 3.0).abs() < 1e-10);
}

#[test]
fn test_end_of_circular_arc_zero_length() {
    let (x_f, y_f, _theta_f) = end_of_circular_arc(0.0, 0.0, 0.0, 1.0, 1.0, 0.0);
    assert!((x_f - 0.0).abs() < 1e-6);
    assert!((y_f - 0.0).abs() < 1e-6);
}

#[test]
fn test_end_of_circular_arc_half_circle() {
    let kappa = 1.0;
    let length = PI;
    let (x_f, y_f, _theta_f) = end_of_circular_arc(0.0, 0.0, 0.0, kappa, 1.0, length);
    assert!((x_f - 0.0).abs() < 1e-6);
    assert!((y_f - 2.0).abs() < 1e-6);
}

#[test]
fn test_end_of_clothoid_zero_length() {
    let (x_f, y_f, theta_f, kappa_f) =
        end_of_clothoid(0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0);
    assert!(x_f.abs() < 1e-10);
    assert!(y_f.abs() < 1e-10);
    assert!(theta_f.abs() < 1e-10);
    assert!(kappa_f.abs() < 1e-10);
}

#[test]
fn test_end_of_clothoid_kappa_change() {
    let sigma = 0.5;
    let length = 2.0;
    let (_, _, _, kappa_f) = end_of_clothoid(0.0, 0.0, 0.0, 0.0, sigma, 1.0, length);
    assert!((kappa_f - sigma * length).abs() < 1e-10);
}
