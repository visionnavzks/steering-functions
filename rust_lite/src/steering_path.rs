use crate::base_state_space::StateSpace;
use crate::dubins::{DubinsDirectionMode, DubinsStateSpace};
use crate::reeds_shepp::ReedsSheppStateSpace;
use crate::state::{Control, State};

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathType {
    Dubins,
    Rs,
}

enum Planner {
    Dubins(DubinsStateSpace),
    Rs(ReedsSheppStateSpace),
}

impl Planner {
    fn get_controls(&self, start: &State, goal: &State) -> Vec<Control> {
        match self {
            Self::Dubins(planner) => planner.get_controls(start, goal),
            Self::Rs(planner) => planner.get_controls(start, goal),
        }
    }

    fn get_path(&self, start: &State, goal: &State) -> Vec<State> {
        match self {
            Self::Dubins(planner) => planner.get_path(start, goal),
            Self::Rs(planner) => planner.get_path(start, goal),
        }
    }

    fn get_all_controls(&self, start: &State, goal: &State) -> Vec<Vec<Control>> {
        match self {
            Self::Dubins(planner) => planner.get_all_controls(start, goal),
            Self::Rs(planner) => planner.get_all_controls(start, goal),
        }
    }

    fn get_all_paths(&self, start: &State, goal: &State) -> Vec<Vec<State>> {
        match self {
            Self::Dubins(planner) => planner.get_all_paths(start, goal),
            Self::Rs(planner) => planner.get_all_paths(start, goal),
        }
    }
}

pub struct SteeringPath {
    pub path_type: PathType,
    pub kappa_max: f64,
    pub discretization: f64,
    pub dubins_direction_mode: DubinsDirectionMode,
}

impl SteeringPath {
    pub fn new(path_type: PathType, kappa_max: f64, discretization: f64) -> Self {
        assert!(kappa_max > 0.0 && discretization > 0.0);
        Self {
            path_type,
            kappa_max,
            discretization,
            dubins_direction_mode: DubinsDirectionMode::ForwardOnly,
        }
    }

    pub fn try_new(path_type: PathType, kappa_max: f64, discretization: f64) -> Result<Self, String> {
        if kappa_max <= 0.0 || discretization <= 0.0 {
            return Err("invalid SteeringPath parameters".to_string());
        }

        Ok(Self {
            path_type,
            kappa_max,
            discretization,
            dubins_direction_mode: DubinsDirectionMode::ForwardOnly,
        })
    }

    pub fn with_dubins_direction_mode(mut self, direction_mode: DubinsDirectionMode) -> Self {
        self.dubins_direction_mode = direction_mode;
        self
    }

    pub fn supported_path_types() -> Vec<PathType> {
        vec![PathType::Dubins, PathType::Rs]
    }

    fn planner(&self) -> Planner {
        match self.path_type {
            PathType::Dubins => Planner::Dubins(DubinsStateSpace::new(
                self.kappa_max,
                self.discretization,
                self.dubins_direction_mode,
            )),
            PathType::Rs => Planner::Rs(ReedsSheppStateSpace::new(self.kappa_max, self.discretization)),
        }
    }

    pub fn compute_shortest_control_sequence(&self, start: &State, goal: &State) -> Vec<Control> {
        self.planner().get_controls(start, goal)
    }

    pub fn compute_shortest_path(&self, start: &State, goal: &State) -> Vec<State> {
        self.planner().get_path(start, goal)
    }

    pub fn compute_all_control_sequences(&self, start: &State, goal: &State) -> Vec<Vec<Control>> {
        self.planner().get_all_controls(start, goal)
    }

    pub fn compute_all_paths(&self, start: &State, goal: &State) -> Vec<Vec<State>> {
        self.planner().get_all_paths(start, goal)
    }
}
