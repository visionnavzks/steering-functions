use core::f64::consts::PI as CORE_PI;

pub const PI: f64 = CORE_PI;
pub const TWO_PI: f64 = 2.0 * PI;

pub fn point_distance(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    (x2 - x1).hypot(y2 - y1)
}

pub fn twopify(alpha: f64) -> f64 {
    alpha - TWO_PI * (alpha / TWO_PI).floor()
}

pub fn pify(alpha: f64) -> f64 {
    let mut v = alpha % TWO_PI;
    if v < -PI {
        v += TWO_PI;
    } else if v > PI {
        v -= TWO_PI;
    }
    v
}
