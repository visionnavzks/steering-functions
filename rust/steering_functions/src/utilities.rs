use core::f64::consts::PI as CORE_PI;

pub const PI: f64 = CORE_PI;
pub const HALF_PI: f64 = PI / 2.0;
pub const TWO_PI: f64 = 2.0 * PI;
pub const SQRT_PI: f64 = 1.772_453_850_905_516_f64;
pub const SQRT_PI_INV: f64 = 0.564_189_583_547_756_3_f64;
pub const SQRT_TWO_PI_INV: f64 = 0.398_942_280_401_432_7_f64;
pub const EPSILON: f64 = 1e-4;

pub fn get_epsilon() -> f64 {
    EPSILON
}

pub fn sgn(x: f64) -> f64 {
    if x < 0.0 { -1.0 } else { 1.0 }
}

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

pub fn polar(x: f64, y: f64) -> (f64, f64) {
    ((x * x + y * y).sqrt(), y.atan2(x))
}

pub fn fresnel(s: f64) -> (f64, f64) {
    let n = 10_000usize;
    let sign = if s < 0.0 { -1.0 } else { 1.0 };
    let s_abs = s.abs();
    if s_abs == 0.0 {
        return (0.0, 0.0);
    }
    let h = s_abs / n as f64;
    let mut sf = 0.0;
    let mut cf = 0.0;
    for i in 0..n {
        let u0 = i as f64 * h;
        let u1 = (i as f64 + 1.0) * h;
        let um = 0.5 * (u0 + u1);
        let f_s0 = (HALF_PI * u0 * u0).sin();
        let f_s1 = (HALF_PI * u1 * u1).sin();
        let f_sm = (HALF_PI * um * um).sin();
        sf += (f_s0 + 4.0 * f_sm + f_s1) * h / 6.0;
        let f_c0 = (HALF_PI * u0 * u0).cos();
        let f_c1 = (HALF_PI * u1 * u1).cos();
        let f_cm = (HALF_PI * um * um).cos();
        cf += (f_c0 + 4.0 * f_cm + f_c1) * h / 6.0;
    }
    (sign * sf, sign * cf)
}

#[allow(clippy::too_many_arguments)]
pub fn end_of_clothoid(
    x_i: f64,
    y_i: f64,
    theta_i: f64,
    kappa_i: f64,
    sigma: f64,
    direction: f64,
    length: f64,
) -> (f64, f64, f64, f64) {
    let sgn_sigma = sgn(sigma);
    let abs_sigma = sigma.abs();
    let sqrt_sigma_inv = 1.0 / abs_sigma.sqrt();
    let k1 = theta_i - 0.5 * direction * kappa_i * kappa_i / sigma;
    let k2 = SQRT_PI_INV * sqrt_sigma_inv * (abs_sigma * length + sgn_sigma * kappa_i);
    let k3 = SQRT_PI_INV * sqrt_sigma_inv * sgn_sigma * kappa_i;
    let cos_k1 = k1.cos();
    let sin_k1 = k1.sin();
    let (fresnel_s_k2, fresnel_c_k2) = fresnel(k2);
    let (fresnel_s_k3, fresnel_c_k3) = fresnel(k3);
    let x_f = x_i
        + SQRT_PI
            * sqrt_sigma_inv
            * (direction * cos_k1 * (fresnel_c_k2 - fresnel_c_k3)
                - sgn_sigma * sin_k1 * (fresnel_s_k2 - fresnel_s_k3));
    let y_f = y_i
        + SQRT_PI
            * sqrt_sigma_inv
            * (direction * sin_k1 * (fresnel_c_k2 - fresnel_c_k3)
                + sgn_sigma * cos_k1 * (fresnel_s_k2 - fresnel_s_k3));
    let theta_f =
        pify(theta_i + kappa_i * direction * length + 0.5 * sigma * direction * length * length);
    let kappa_f = kappa_i + sigma * length;
    (x_f, y_f, theta_f, kappa_f)
}

pub fn end_of_circular_arc(
    x_i: f64,
    y_i: f64,
    theta_i: f64,
    kappa: f64,
    direction: f64,
    length: f64,
) -> (f64, f64, f64) {
    let x_f = x_i + (1.0 / kappa) * (-theta_i.sin() + (theta_i + direction * length * kappa).sin());
    let y_f = y_i + (1.0 / kappa) * (theta_i.cos() - (theta_i + direction * length * kappa).cos());
    let theta_f = pify(theta_i + kappa * direction * length);
    (x_f, y_f, theta_f)
}

pub fn end_of_straight_line(
    x_i: f64,
    y_i: f64,
    theta: f64,
    direction: f64,
    length: f64,
) -> (f64, f64) {
    (
        x_i + direction * length * theta.cos(),
        y_i + direction * length * theta.sin(),
    )
}

pub fn global_frame_change(x: f64, y: f64, theta: f64, local_x: f64, local_y: f64) -> (f64, f64) {
    let sin_th = theta.sin();
    let cos_th = theta.cos();
    (
        local_x * cos_th - local_y * sin_th + x,
        local_x * sin_th + local_y * cos_th + y,
    )
}

pub fn local_frame_change(x: f64, y: f64, theta: f64, global_x: f64, global_y: f64) -> (f64, f64) {
    let sin_th = theta.sin();
    let cos_th = theta.cos();
    (
        (global_x - x) * cos_th + (global_y - y) * sin_th,
        -(global_x - x) * sin_th + (global_y - y) * cos_th,
    )
}

pub fn array_index_min(array: &[f64]) -> usize {
    let mut min = array[0];
    let mut index_min = 0usize;
    for (i, &v) in array.iter().enumerate().skip(1) {
        if v < min {
            index_min = i;
            min = v;
        }
    }
    index_min
}

pub fn double_array_init(size: usize, value: f64) -> Vec<f64> {
    vec![value; size]
}

pub fn pointer_array_init<T>(size: usize) -> Vec<Option<T>> {
    (0..size).map(|_| None).collect()
}
