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

pub struct SteeringPath {
    path_type: PathType,
    kappa_max: f64,
    discretization: f64,
    dubins_direction_mode: DubinsDirectionMode,
    planner: Box<dyn StateSpace>,
}

impl SteeringPath {
    pub fn new(path_type: PathType, kappa_max: f64, discretization: f64) -> Self {
        assert!(kappa_max > 0.0 && discretization > 0.0);
        let dubins_direction_mode = DubinsDirectionMode::ForwardOnly;
        Self {
            path_type,
            kappa_max,
            discretization,
            dubins_direction_mode,
            planner: Self::build_planner(path_type, kappa_max, discretization, dubins_direction_mode),
        }
    }

    pub fn try_new(path_type: PathType, kappa_max: f64, discretization: f64) -> Result<Self, String> {
        if kappa_max <= 0.0 || discretization <= 0.0 {
            return Err("invalid SteeringPath parameters".to_string());
        }

        Ok(Self::new(path_type, kappa_max, discretization))
    }

    pub fn with_dubins_direction_mode(mut self, direction_mode: DubinsDirectionMode) -> Self {
        self.set_dubins_direction_mode(direction_mode);
        self
    }

    pub fn path_type(&self) -> PathType {
        self.path_type
    }

    pub fn kappa_max(&self) -> f64 {
        self.kappa_max
    }

    pub fn discretization(&self) -> f64 {
        self.discretization
    }

    pub fn dubins_direction_mode(&self) -> DubinsDirectionMode {
        self.dubins_direction_mode
    }

    pub fn set_path_type(&mut self, path_type: PathType) {
        self.path_type = path_type;
        self.rebuild_planner();
    }

    pub fn set_kappa_max(&mut self, kappa_max: f64) -> Result<(), String> {
        if kappa_max <= 0.0 {
            return Err("invalid kappa_max".to_string());
        }
        self.kappa_max = kappa_max;
        self.rebuild_planner();
        Ok(())
    }

    pub fn set_discretization(&mut self, discretization: f64) -> Result<(), String> {
        if discretization <= 0.0 {
            return Err("invalid discretization".to_string());
        }
        self.discretization = discretization;
        self.rebuild_planner();
        Ok(())
    }

    pub fn set_dubins_direction_mode(&mut self, direction_mode: DubinsDirectionMode) {
        self.dubins_direction_mode = direction_mode;
        self.rebuild_planner();
    }

    pub fn supported_path_types() -> Vec<PathType> {
        vec![PathType::Dubins, PathType::Rs]
    }

    fn build_planner(
        path_type: PathType,
        kappa_max: f64,
        discretization: f64,
        dubins_direction_mode: DubinsDirectionMode,
    ) -> Box<dyn StateSpace> {
        match path_type {
            PathType::Dubins => Box::new(DubinsStateSpace::new(
                kappa_max,
                discretization,
                dubins_direction_mode,
            )),
            PathType::Rs => Box::new(ReedsSheppStateSpace::new(kappa_max, discretization)),
        }
    }

    fn rebuild_planner(&mut self) {
        self.planner = Self::build_planner(
            self.path_type,
            self.kappa_max,
            self.discretization,
            self.dubins_direction_mode,
        );
    }

    fn planner(&self) -> &dyn StateSpace {
        self.planner.as_ref()
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
