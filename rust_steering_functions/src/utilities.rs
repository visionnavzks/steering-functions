pub const EPSILON: f64 = 1e-4;
pub const PI: f64 = std::f64::consts::PI;
pub const HALF_PI: f64 = PI / 2.0;
pub const TWO_PI: f64 = PI * 2.0;
pub const SQRT_PI: f64 = 1.772_453_850_905_516;
pub const SQRT_PI_INV: f64 = 0.564_189_583_547_756_3;
pub const SQRT_TWO_PI_INV: f64 = 0.398_942_280_401_432_7;

#[inline]
pub fn sgn(x: f64) -> f64 {
    if x < 0.0 {
        -1.0
    } else {
        1.0
    }
}

#[inline]
pub fn point_distance(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt()
}

#[inline]
pub fn polar(x: f64, y: f64) -> (f64, f64) {
    let r = (x * x + y * y).sqrt();
    let theta = y.atan2(x);
    (r, theta)
}

#[inline]
pub fn twopify(alpha: f64) -> f64 {
    alpha - TWO_PI * (alpha / TWO_PI).floor()
}

#[inline]
pub fn pify(alpha: f64) -> f64 {
    let mut v = alpha.rem_euclid(TWO_PI);
    if v < -PI {
        v += TWO_PI;
    } else if v > PI {
        v -= TWO_PI;
    }
    v
}

const CHEBEV_A: [f64; 18] = [
    0.764_351_386_641_86,
    -0.431_355_475_476_601_8,
    0.432_881_999_797_266_55,
    -0.269_733_103_383_871_13,
    0.084_160_453_208_769_36,
    -0.015_465_244_844_613_82,
    0.001_878_554_234_398_220_2,
    -0.00016264977618887547,
    0.000_010_573_976_563_832_6,
    -0.00000053609339889243,
    0.00000002181658454933,
    -0.00000000072901621186,
    0.00000000002037332548,
    -0.00000000000048344033,
    0.00000000000000986533,
    -0.00000000000000017502,
    0.00000000000000000272,
    -0.00000000000000000004,
];

const CHEBEV_B: [f64; 17] = [
    0.630_414_043_145_705_4,
    -0.423_445_114_057_053_3,
    0.376_171_726_433_436_55,
    -0.162_494_891_545_095_67,
    0.038_222_557_786_330_09,
    -0.005_645_634_771_321_909,
    0.000_574_549_519_768_973_7,
    -0.00004287071532102004,
    0.00000245120749923299,
    -0.00000011098841840868,
    0.00000000408249731696,
    -0.00000000012449830219,
    0.00000000000320048425,
    -0.00000000000007032416,
    0.00000000000000133638,
    -0.00000000000000002219,
    0.00000000000000000032,
];

const CHEBEV_E: [f64; 41] = [
    0.974_627_790_932_968_3,
    -0.024_247_018_739_693_215,
    0.001_034_009_068_429_773_1,
    -0.00008052450246908016,
    0.00000905962481966582,
    -0.00000131016996757743,
    0.00000022770820391497,
    -0.00000004558623552026,
    0.00000001021567537083,
    -0.00000000251114508133,
    0.00000000066704761275,
    -0.00000000018931512852,
    0.00000000005689898935,
    -0.00000000001798219359,
    0.00000000000594162963,
    -0.00000000000204285065,
    0.00000000000072797580,
    -0.00000000000026797428,
    0.00000000000010160694,
    -0.00000000000003958559,
    0.00000000000001581262,
    -0.00000000000000646411,
    0.00000000000000269981,
    -0.00000000000000115038,
    0.00000000000000049942,
    -0.00000000000000022064,
    0.00000000000000009910,
    -0.00000000000000004520,
    0.00000000000000002092,
    -0.00000000000000000982,
    0.00000000000000000467,
    -0.00000000000000000225,
    0.00000000000000000110,
    -0.00000000000000000054,
    0.00000000000000000027,
    -0.00000000000000000014,
    0.00000000000000000007,
    -0.00000000000000000004,
    0.00000000000000000002,
    -0.00000000000000000001,
    0.00000000000000000001,
];

const CHEBEV_F: [f64; 35] = [
    0.994_615_451_794_079_3,
    -0.005_242_767_660_842_972,
    0.000_133_258_642_298_839_1,
    -0.00000770856452642713,
    0.00000070848077032045,
    -0.00000008812517411602,
    0.00000001359784717148,
    -0.00000000246858295747,
    0.00000000050925789921,
    -0.00000000011653400634,
    0.00000000002906578309,
    -0.00000000000779847361,
    0.00000000000222802542,
    -0.00000000000067239338,
    0.00000000000021296411,
    -0.00000000000007041482,
    0.00000000000002419805,
    -0.00000000000000861080,
    0.00000000000000316287,
    -0.00000000000000119596,
    0.00000000000000046444,
    -0.00000000000000018485,
    0.00000000000000007527,
    -0.00000000000000003131,
    0.00000000000000001328,
    -0.00000000000000000574,
    0.00000000000000000252,
    -0.00000000000000000113,
    0.00000000000000000051,
    -0.00000000000000000024,
    0.00000000000000000011,
    -0.00000000000000000005,
    0.00000000000000000002,
    -0.00000000000000000001,
    0.00000000000000000001,
];

fn fresnel_0_8(x: f64) -> (f64, f64) {
    let quarter_x = 0.25 * x;
    let arg = 0.03125 * x * x - 1.0;
    let t0 = 1.0;
    let t1 = 0.125 * x;
    let t2 = arg;
    let t3 = quarter_x * t2 - t1;
    let mut a = CHEBEV_A[0] * t0 + CHEBEV_A[1] * t2;
    let mut b = CHEBEV_B[0] * t1 + CHEBEV_B[1] * t3;
    let mut t2n_m4 = t0;
    let mut t2n_m2 = t2;
    let mut t2n_m1 = t3;

    for n in 2..17 {
        let t2n = 2.0 * arg * t2n_m2 - t2n_m4;
        let t2n_p1 = quarter_x * t2n - t2n_m1;
        a += CHEBEV_A[n] * t2n;
        b += CHEBEV_B[n] * t2n_p1;
        t2n_m4 = t2n_m2;
        t2n_m2 = t2n;
        t2n_m1 = t2n_p1;
    }
    let t34 = 2.0 * arg * t2n_m2 - t2n_m4;
    a += CHEBEV_A[17] * t34;

    let sqrt_x = x.sqrt();
    let c_f = SQRT_TWO_PI_INV * sqrt_x * a;
    let s_f = SQRT_TWO_PI_INV * sqrt_x * b;
    (s_f, c_f)
}

fn fresnel_8_inf(x: f64) -> (f64, f64) {
    let arg = 128.0 / (x * x) - 1.0;
    let t0 = 1.0;
    let t2 = arg;
    let mut e = CHEBEV_E[0] * t0 + CHEBEV_E[1] * t2;
    let mut f = CHEBEV_F[0] * t0 + CHEBEV_F[1] * t2;
    let mut t2n_m4 = t0;
    let mut t2n_m2 = t2;

    for n in 2..35 {
        let t2n = 2.0 * arg * t2n_m2 - t2n_m4;
        e += CHEBEV_E[n] * t2n;
        f += CHEBEV_F[n] * t2n;
        t2n_m4 = t2n_m2;
        t2n_m2 = t2n;
    }
    for n in 35..41 {
        let t2n = 2.0 * arg * t2n_m2 - t2n_m4;
        e += CHEBEV_E[n] * t2n;
        t2n_m4 = t2n_m2;
        t2n_m2 = t2n;
    }

    let sin_x = x.sin();
    let cos_x = x.cos();
    let sqrt_x = x.sqrt();
    let c_f = 0.5 - SQRT_TWO_PI_INV * (e * cos_x / (2.0 * x) - f * sin_x) / sqrt_x;
    let s_f = 0.5 - SQRT_TWO_PI_INV * (e * sin_x / (2.0 * x) + f * cos_x) / sqrt_x;
    (s_f, c_f)
}

pub fn fresnel(s: f64) -> (f64, f64) {
    let x = HALF_PI * s * s;
    let (mut s_f, mut c_f) = if x <= 8.0 {
        fresnel_0_8(x)
    } else {
        fresnel_8_inf(x)
    };
    if s < 0.0 {
        s_f = -s_f;
        c_f = -c_f;
    }
    (s_f, c_f)
}

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
    let x_f = x_i + direction * length * theta.cos();
    let y_f = y_i + direction * length * theta.sin();
    (x_f, y_f)
}

pub fn global_frame_change(x: f64, y: f64, theta: f64, local_x: f64, local_y: f64) -> (f64, f64) {
    let sin_th = theta.sin();
    let cos_th = theta.cos();
    let global_x = local_x * cos_th - local_y * sin_th + x;
    let global_y = local_x * sin_th + local_y * cos_th + y;
    (global_x, global_y)
}

pub fn local_frame_change(x: f64, y: f64, theta: f64, global_x: f64, global_y: f64) -> (f64, f64) {
    let sin_th = theta.sin();
    let cos_th = theta.cos();
    let local_x = (global_x - x) * cos_th + (global_y - y) * sin_th;
    let local_y = -(global_x - x) * sin_th + (global_y - y) * cos_th;
    (local_x, local_y)
}
