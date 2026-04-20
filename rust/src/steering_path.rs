
use crate::base_state_space::StateSpace;
use crate::cc_dubins::{
    CC00DubinsStateSpace, CC0pmDubinsStateSpace, CCDubinsStateSpace,
    CCpm0DubinsStateSpace, CCpmpmDubinsStateSpace,
};
use crate::cc_reeds_shepp::CcReedsSheppStateSpace;
use crate::dubins::DubinsStateSpace;
use crate::hc_reeds_shepp::{
    Hc00RsStateSpace, Hc0pmRsStateSpace, HcRsStateSpace, Hcpm0RsStateSpace,
    HcpmpmRsStateSpace,
};
use crate::reeds_shepp::ReedsSheppStateSpace;
use crate::state::{Control, State};

/// Supported path planning algorithm types.
#[repr(u8)]
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

enum Planner {
    CcDubins(CCDubinsStateSpace),
    Cc00Dubins(CC00DubinsStateSpace),
    Cc0pmDubins(CC0pmDubinsStateSpace),
    Ccpm0Dubins(CCpm0DubinsStateSpace),
    CcpmpmDubins(CCpmpmDubinsStateSpace),
    Dubins(DubinsStateSpace),
    Cc00Rs(CcReedsSheppStateSpace),
    HcRs(HcRsStateSpace),
    Hc00Rs(Hc00RsStateSpace),
    Hc0pmRs(Hc0pmRsStateSpace),
    Hcpm0Rs(Hcpm0RsStateSpace),
    HcpmpmRs(HcpmpmRsStateSpace),
    Rs(ReedsSheppStateSpace),
}

impl Planner {
    fn get_controls(&self, start: &State, goal: &State) -> Vec<Control> {
        match self {
            Self::CcDubins(planner) => planner.get_controls(start, goal),
            Self::Cc00Dubins(planner) => planner.get_controls(start, goal),
            Self::Cc0pmDubins(planner) => planner.get_controls(start, goal),
            Self::Ccpm0Dubins(planner) => planner.get_controls(start, goal),
            Self::CcpmpmDubins(planner) => planner.get_controls(start, goal),
            Self::Dubins(planner) => planner.get_controls(start, goal),
            Self::Cc00Rs(planner) => planner.get_controls(start, goal),
            Self::HcRs(planner) => planner.get_controls(start, goal),
            Self::Hc00Rs(planner) => planner.get_controls(start, goal),
            Self::Hc0pmRs(planner) => planner.get_controls(start, goal),
            Self::Hcpm0Rs(planner) => planner.get_controls(start, goal),
            Self::HcpmpmRs(planner) => planner.get_controls(start, goal),
            Self::Rs(planner) => planner.get_controls(start, goal),
        }
    }

    fn get_path(&self, start: &State, goal: &State) -> Vec<State> {
        match self {
            Self::CcDubins(planner) => planner.get_path(start, goal),
            Self::Cc00Dubins(planner) => planner.get_path(start, goal),
            Self::Cc0pmDubins(planner) => planner.get_path(start, goal),
            Self::Ccpm0Dubins(planner) => planner.get_path(start, goal),
            Self::CcpmpmDubins(planner) => planner.get_path(start, goal),
            Self::Dubins(planner) => planner.get_path(start, goal),
            Self::Cc00Rs(planner) => planner.get_path(start, goal),
            Self::HcRs(planner) => planner.get_path(start, goal),
            Self::Hc00Rs(planner) => planner.get_path(start, goal),
            Self::Hc0pmRs(planner) => planner.get_path(start, goal),
            Self::Hcpm0Rs(planner) => planner.get_path(start, goal),
            Self::HcpmpmRs(planner) => planner.get_path(start, goal),
            Self::Rs(planner) => planner.get_path(start, goal),
        }
    }

    fn get_all_controls(&self, start: &State, goal: &State) -> Vec<Vec<Control>> {
        match self {
            Self::CcDubins(planner) => planner.get_all_controls(start, goal),
            Self::Cc00Dubins(planner) => planner.get_all_controls(start, goal),
            Self::Cc0pmDubins(planner) => planner.get_all_controls(start, goal),
            Self::Ccpm0Dubins(planner) => planner.get_all_controls(start, goal),
            Self::CcpmpmDubins(planner) => planner.get_all_controls(start, goal),
            Self::Dubins(planner) => planner.get_all_controls(start, goal),
            Self::Cc00Rs(planner) => planner.get_all_controls(start, goal),
            Self::HcRs(planner) => planner.get_all_controls(start, goal),
            Self::Hc00Rs(planner) => planner.get_all_controls(start, goal),
            Self::Hc0pmRs(planner) => planner.get_all_controls(start, goal),
            Self::Hcpm0Rs(planner) => planner.get_all_controls(start, goal),
            Self::HcpmpmRs(planner) => planner.get_all_controls(start, goal),
            Self::Rs(planner) => planner.get_all_controls(start, goal),
        }
    }

    fn get_all_paths(&self, start: &State, goal: &State) -> Vec<Vec<State>> {
        match self {
            Self::CcDubins(planner) => planner.get_all_paths(start, goal),
            Self::Cc00Dubins(planner) => planner.get_all_paths(start, goal),
            Self::Cc0pmDubins(planner) => planner.get_all_paths(start, goal),
            Self::Ccpm0Dubins(planner) => planner.get_all_paths(start, goal),
            Self::CcpmpmDubins(planner) => planner.get_all_paths(start, goal),
            Self::Dubins(planner) => planner.get_all_paths(start, goal),
            Self::Cc00Rs(planner) => planner.get_all_paths(start, goal),
            Self::HcRs(planner) => planner.get_all_paths(start, goal),
            Self::Hc00Rs(planner) => planner.get_all_paths(start, goal),
            Self::Hc0pmRs(planner) => planner.get_all_paths(start, goal),
            Self::Hcpm0Rs(planner) => planner.get_all_paths(start, goal),
            Self::HcpmpmRs(planner) => planner.get_all_paths(start, goal),
            Self::Rs(planner) => planner.get_all_paths(start, goal),
        }
    }
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

    pub fn try_new(
        path_type: PathType,
        kappa_max: f64,
        sigma_max: f64,
        discretization: f64,
    ) -> Result<Self, String> {
        if kappa_max <= 0.0 || sigma_max < 0.0 || discretization <= 0.0 {
            return Err("invalid SteeringPath parameters".to_string());
        }
        Ok(Self { path_type, kappa_max, sigma_max, discretization })
    }

    pub fn supported_path_types() -> Vec<PathType> {
        vec![
            PathType::CcDubins,
            PathType::Cc00Dubins,
            PathType::Cc0pmDubins,
            PathType::Ccpm0Dubins,
            PathType::CcpmpmDubins,
            PathType::Dubins,
            PathType::Cc00Rs,
            PathType::HcRs,
            PathType::Hc00Rs,
            PathType::Hc0pmRs,
            PathType::Hcpm0Rs,
            PathType::HcpmpmRs,
            PathType::Rs,
        ]
    }

    pub fn is_supported(path_type: PathType) -> bool {
        Self::supported_path_types().contains(&path_type)
    }

    fn planner(&self) -> Result<Planner, String> {
        let kappa = self.kappa_max;
        let sigma = self.sigma_max;
        let disc = self.discretization;

        match self.path_type {
            PathType::CcDubins => Ok(Planner::CcDubins(CCDubinsStateSpace::new(kappa, sigma, disc, true))),
            PathType::Cc00Dubins => Ok(Planner::Cc00Dubins(CC00DubinsStateSpace::new(kappa, sigma, disc, true))),
            PathType::Cc0pmDubins => Ok(Planner::Cc0pmDubins(CC0pmDubinsStateSpace::new(kappa, sigma, disc, true))),
            PathType::Ccpm0Dubins => Ok(Planner::Ccpm0Dubins(CCpm0DubinsStateSpace::new(kappa, sigma, disc, true))),
            PathType::CcpmpmDubins => Ok(Planner::CcpmpmDubins(CCpmpmDubinsStateSpace::new(kappa, sigma, disc, true))),
            PathType::Dubins => Ok(Planner::Dubins(DubinsStateSpace::new(kappa, disc, true))),
            PathType::Cc00Rs => Ok(Planner::Cc00Rs(CcReedsSheppStateSpace::new(kappa, sigma, disc))),
            PathType::HcRs => Ok(Planner::HcRs(HcRsStateSpace::new(kappa, sigma, disc))),
            PathType::Hc00Rs => Ok(Planner::Hc00Rs(Hc00RsStateSpace::new(kappa, sigma, disc))),
            PathType::Hc0pmRs => Ok(Planner::Hc0pmRs(Hc0pmRsStateSpace::new(kappa, sigma, disc))),
            PathType::Hcpm0Rs => Ok(Planner::Hcpm0Rs(Hcpm0RsStateSpace::new(kappa, sigma, disc))),
            PathType::HcpmpmRs => Ok(Planner::HcpmpmRs(HcpmpmRsStateSpace::new(kappa, sigma, disc))),
            PathType::Rs => Ok(Planner::Rs(ReedsSheppStateSpace::new(kappa, disc))),
            PathType::None => Err("PathType::None is not a runnable planner".to_string()),
        }
    }

    pub fn compute_shortest_control_sequence(
        &self,
        start: &State,
        goal: &State,
    ) -> Result<Vec<Control>, String> {
        Ok(self.planner()?.get_controls(start, goal))
    }

    pub fn compute_shortest_path(
        &self,
        start: &State,
        goal: &State,
    ) -> Result<Vec<State>, String> {
        Ok(self.planner()?.get_path(start, goal))
    }

    pub fn compute_all_control_sequences(
        &self,
        start: &State,
        goal: &State,
    ) -> Result<Vec<Vec<Control>>, String> {
        Ok(self.planner()?.get_all_controls(start, goal))
    }

    pub fn compute_all_paths(
        &self,
        start: &State,
        goal: &State,
    ) -> Result<Vec<Vec<State>>, String> {
        Ok(self.planner()?.get_all_paths(start, goal))
    }
}
