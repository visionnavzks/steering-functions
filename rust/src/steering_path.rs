use crate::state::{State, Control};

/// Supported path planning algorithm types.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathType {
    None,
    CcDubins,
    Cc00Dubins,
    Cc0pmDubins,
    Ccpm0Dubins,
    CcpmpmDubins,
    Dubins,
    Cc00Rs,
    HcRs,
    Hc00Rs,
    Hc0pmRs,
    Hcpm0Rs,
    HcpmpmRs,
    Rs,
}

/// High-level wrapper for steering path computation.
pub struct SteeringPath {
    pub path_type: PathType,
    pub kappa_max: f64,
    pub sigma_max: f64,
    pub discretization: f64,
}

impl SteeringPath {
    pub fn new(path_type: PathType, kappa_max: f64, sigma_max: f64, discretization: f64) -> Self {
        assert!(kappa_max > 0.0 && sigma_max >= 0.0 && discretization > 0.0);
        Self { path_type, kappa_max, sigma_max, discretization }
    }
}
