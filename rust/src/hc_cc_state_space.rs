use crate::hc_cc_circle::HcCcCircleParam;
use crate::state::{State, Control};
use crate::utilities::{get_epsilon, end_of_clothoid, point_distance};

/// Shared parameters for HC/CC state spaces.
#[derive(Clone, Debug)]
pub struct HcCcStateSpaceParams {
    pub kappa_: f64,
    pub sigma_: f64,
    pub hc_cc_circle_param_: HcCcCircleParam,
}

impl HcCcStateSpaceParams {
    pub fn new(kappa: f64, sigma: f64) -> Self {
        assert!(kappa > 0.0 && sigma > 0.0);
        let mut param = HcCcCircleParam::default();

        let length_min = kappa / sigma;
        let (x_i, y_i, theta_i) = if length_min > get_epsilon() {
            let (xi, yi, ti, _) = end_of_clothoid(0.0, 0.0, 0.0, 0.0, sigma, 1.0, length_min);
            (xi, yi, ti)
        } else {
            (0.0, 0.0, 0.0)
        };

        let xc = x_i - theta_i.sin() / kappa;
        let yc = y_i + theta_i.cos() / kappa;
        let radius = point_distance(xc, yc, 0.0, 0.0);

        let mu = xc.abs().atan2(yc.abs());
        let sin_mu = mu.sin();
        let cos_mu = mu.cos();
        let delta_min = 0.5 * kappa * kappa / sigma;

        param.set_param(kappa, sigma, radius, mu, sin_mu, cos_mu, delta_min);

        Self {
            kappa_: kappa,
            sigma_: sigma,
            hc_cc_circle_param_: param,
        }
    }
}

/// Filter negligible control segments and integrate the path.
pub fn hc_cc_get_path(
    state1: &State,
    controls_raw: Vec<Control>,
    disc: f64,
    integrate: impl Fn(&State, &[Control], f64) -> Vec<State>,
) -> Vec<State> {
    let controls: Vec<Control> = controls_raw
        .into_iter()
        .filter(|c| c.delta_s.abs() > 1e-6)
        .collect();
    integrate(state1, &controls, disc)
}
