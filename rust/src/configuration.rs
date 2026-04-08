use crate::utilities::{twopify, point_distance, get_epsilon};

/// Kinematic configuration: position, orientation, curvature.
#[derive(Clone, Copy, Debug, Default)]
pub struct Configuration {
    pub x: f64,
    pub y: f64,
    pub theta: f64,
    pub kappa: f64,
}

impl Configuration {
    pub fn new(x: f64, y: f64, theta: f64, kappa: f64) -> Self {
        Self { x, y, theta: twopify(theta), kappa }
    }
}

pub fn configuration_distance(q1: &Configuration, q2: &Configuration) -> f64 {
    point_distance(q1.x, q1.y, q2.x, q2.y)
}

pub fn configuration_aligned(q1: &Configuration, q2: &Configuration) -> bool {
    if (q2.theta - q1.theta).abs() > get_epsilon() {
        return false;
    }
    let angle = twopify((q2.y - q1.y).atan2(q2.x - q1.x));
    (angle - q1.theta).abs() <= get_epsilon()
}

pub fn configuration_equal(q1: &Configuration, q2: &Configuration) -> bool {
    if (q2.theta - q1.theta).abs() > get_epsilon() {
        return false;
    }
    configuration_distance(q1, q2) <= get_epsilon()
}
