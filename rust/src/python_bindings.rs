use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;

use crate::state::{Control, State};
use crate::steering_path::{PathType, SteeringPath};

#[pyclass(name = "State", module = "steering_functions_rust")]
#[derive(Clone)]
pub struct PyState {
    #[pyo3(get, set)]
    pub x: f64,
    #[pyo3(get, set)]
    pub y: f64,
    #[pyo3(get, set)]
    pub theta: f64,
    #[pyo3(get, set)]
    pub kappa: f64,
    #[pyo3(get, set)]
    pub sigma: f64,
    #[pyo3(get, set)]
    pub d: f64,
    #[pyo3(get, set)]
    pub s: f64,
    #[pyo3(get, set)]
    pub vel: f64,
    #[pyo3(get, set)]
    pub acc: f64,
    #[pyo3(get, set)]
    pub time: f64,
    #[pyo3(get, set)]
    pub fork_y: f64,
}

#[pymethods]
impl PyState {
    #[new]
    #[pyo3(signature = (
        x=0.0,
        y=0.0,
        theta=0.0,
        kappa=0.0,
        sigma=0.0,
        d=0.0,
        s=0.0,
        vel=0.0,
        acc=0.0,
        time=0.0,
        fork_y=0.0,
    ))]
    fn new(
        x: f64,
        y: f64,
        theta: f64,
        kappa: f64,
        sigma: f64,
        d: f64,
        s: f64,
        vel: f64,
        acc: f64,
        time: f64,
        fork_y: f64,
    ) -> Self {
        Self { x, y, theta, kappa, sigma, d, s, vel, acc, time, fork_y }
    }

    fn __repr__(&self) -> String {
        format!(
            "State(x={:.6}, y={:.6}, theta={:.6}, kappa={:.6})",
            self.x, self.y, self.theta, self.kappa
        )
    }
}

impl From<PyState> for State {
    fn from(value: PyState) -> Self {
        Self {
            x: value.x,
            y: value.y,
            theta: value.theta,
            kappa: value.kappa,
            sigma: value.sigma,
            d: value.d,
            s: value.s,
            vel: value.vel,
            acc: value.acc,
            time: value.time,
            fork_y: value.fork_y,
        }
    }
}

impl From<State> for PyState {
    fn from(value: State) -> Self {
        Self {
            x: value.x,
            y: value.y,
            theta: value.theta,
            kappa: value.kappa,
            sigma: value.sigma,
            d: value.d,
            s: value.s,
            vel: value.vel,
            acc: value.acc,
            time: value.time,
            fork_y: value.fork_y,
        }
    }
}

#[pyclass(name = "Control", module = "steering_functions_rust")]
#[derive(Clone)]
pub struct PyControl {
    #[pyo3(get, set)]
    pub delta_s: f64,
    #[pyo3(get, set)]
    pub kappa: f64,
    #[pyo3(get, set)]
    pub sigma: f64,
}

#[pymethods]
impl PyControl {
    #[new]
    #[pyo3(signature = (delta_s=0.0, kappa=0.0, sigma=0.0))]
    fn new(delta_s: f64, kappa: f64, sigma: f64) -> Self {
        Self { delta_s, kappa, sigma }
    }

    fn __repr__(&self) -> String {
        format!(
            "Control(delta_s={:.6}, kappa={:.6}, sigma={:.6})",
            self.delta_s, self.kappa, self.sigma
        )
    }
}

impl From<Control> for PyControl {
    fn from(value: Control) -> Self {
        Self {
            delta_s: value.delta_s,
            kappa: value.kappa,
            sigma: value.sigma,
        }
    }
}

#[pyclass(name = "PathType", eq, eq_int, module = "steering_functions_rust")]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum PyPathType {
    NONE = 0,
    CC_DUBINS = 1,
    CC00_DUBINS = 2,
    CC0PM_DUBINS = 3,
    CCPM0_DUBINS = 4,
    CCPMPM_DUBINS = 5,
    DUBINS = 6,
    CC00_RS = 7,
    HC_RS = 8,
    HC00_RS = 9,
    HC0PM_RS = 10,
    HCPM0_RS = 11,
    HCPMPM_RS = 12,
    RS = 13,
}

impl From<PyPathType> for PathType {
    fn from(value: PyPathType) -> Self {
        match value {
            PyPathType::NONE => PathType::None,
            PyPathType::CC_DUBINS => PathType::CcDubins,
            PyPathType::CC00_DUBINS => PathType::Cc00Dubins,
            PyPathType::CC0PM_DUBINS => PathType::Cc0pmDubins,
            PyPathType::CCPM0_DUBINS => PathType::Ccpm0Dubins,
            PyPathType::CCPMPM_DUBINS => PathType::CcpmpmDubins,
            PyPathType::DUBINS => PathType::Dubins,
            PyPathType::CC00_RS => PathType::Cc00Rs,
            PyPathType::HC_RS => PathType::HcRs,
            PyPathType::HC00_RS => PathType::Hc00Rs,
            PyPathType::HC0PM_RS => PathType::Hc0pmRs,
            PyPathType::HCPM0_RS => PathType::Hcpm0Rs,
            PyPathType::HCPMPM_RS => PathType::HcpmpmRs,
            PyPathType::RS => PathType::Rs,
        }
    }
}

impl From<PathType> for PyPathType {
    fn from(value: PathType) -> Self {
        match value {
            PathType::None => PyPathType::NONE,
            PathType::CcDubins => PyPathType::CC_DUBINS,
            PathType::Cc00Dubins => PyPathType::CC00_DUBINS,
            PathType::Cc0pmDubins => PyPathType::CC0PM_DUBINS,
            PathType::Ccpm0Dubins => PyPathType::CCPM0_DUBINS,
            PathType::CcpmpmDubins => PyPathType::CCPMPM_DUBINS,
            PathType::Dubins => PyPathType::DUBINS,
            PathType::Cc00Rs => PyPathType::CC00_RS,
            PathType::HcRs => PyPathType::HC_RS,
            PathType::Hc00Rs => PyPathType::HC00_RS,
            PathType::Hc0pmRs => PyPathType::HC0PM_RS,
            PathType::Hcpm0Rs => PyPathType::HCPM0_RS,
            PathType::HcpmpmRs => PyPathType::HCPMPM_RS,
            PathType::Rs => PyPathType::RS,
        }
    }
}

#[pyclass(name = "SteeringPath", module = "steering_functions_rust")]
pub struct PySteeringPath {
    inner: SteeringPath,
}

fn map_steering_error(err: String) -> PyErr {
    if err.contains("not implemented") {
        PyNotImplementedError::new_err(err)
    } else {
        PyValueError::new_err(err)
    }
}

fn controls_to_py(controls: Vec<Control>) -> Vec<PyControl> {
    controls.into_iter().map(Into::into).collect()
}

fn states_to_py(states: Vec<State>) -> Vec<PyState> {
    states.into_iter().map(Into::into).collect()
}

#[pymethods]
impl PySteeringPath {
    #[new]
    #[pyo3(signature = (path_type, kappa_max, sigma_max, discretization))]
    fn new(
        path_type: PyPathType,
        kappa_max: f64,
        sigma_max: f64,
        discretization: f64,
    ) -> PyResult<Self> {
        let inner = SteeringPath::try_new(path_type.into(), kappa_max, sigma_max, discretization)
            .map_err(map_steering_error)?;
        if !SteeringPath::is_supported(inner.path_type) {
            return Err(PyNotImplementedError::new_err(format!(
                "{:?} is not implemented in the Rust port yet",
                inner.path_type
            )));
        }
        Ok(Self { inner })
    }

    #[getter]
    fn path_type(&self) -> PyPathType {
        self.inner.path_type.into()
    }

    #[getter]
    fn kappa_max(&self) -> f64 {
        self.inner.kappa_max
    }

    #[getter]
    fn sigma_max(&self) -> f64 {
        self.inner.sigma_max
    }

    #[getter]
    fn discretization(&self) -> f64 {
        self.inner.discretization
    }

    fn compute_shortest_control_sequence(
        &self,
        start: PyState,
        goal: PyState,
    ) -> PyResult<Vec<PyControl>> {
        let controls = self
            .inner
            .compute_shortest_control_sequence(&start.into(), &goal.into())
            .map_err(map_steering_error)?;
        Ok(controls_to_py(controls))
    }

    fn compute_shortest_path(&self, start: PyState, goal: PyState) -> PyResult<Vec<PyState>> {
        let path = self
            .inner
            .compute_shortest_path(&start.into(), &goal.into())
            .map_err(map_steering_error)?;
        Ok(states_to_py(path))
    }

    fn compute_all_control_sequences(
        &self,
        start: PyState,
        goal: PyState,
    ) -> PyResult<Vec<Vec<PyControl>>> {
        let sequences = self
            .inner
            .compute_all_control_sequences(&start.into(), &goal.into())
            .map_err(map_steering_error)?;
        Ok(sequences.into_iter().map(controls_to_py).collect())
    }

    fn compute_all_paths(&self, start: PyState, goal: PyState) -> PyResult<Vec<Vec<PyState>>> {
        let paths = self
            .inner
            .compute_all_paths(&start.into(), &goal.into())
            .map_err(map_steering_error)?;
        Ok(paths.into_iter().map(states_to_py).collect())
    }

    #[staticmethod]
    fn supported_path_types() -> Vec<PyPathType> {
        SteeringPath::supported_path_types()
            .into_iter()
            .map(Into::into)
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "SteeringPath(path_type={:?}, kappa_max={}, sigma_max={}, discretization={})",
            self.inner.path_type,
            self.inner.kappa_max,
            self.inner.sigma_max,
            self.inner.discretization,
        )
    }
}

#[pyfunction]
fn supported_path_types() -> Vec<PyPathType> {
    PySteeringPath::supported_path_types()
}

#[pymodule]
fn _core(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyState>()?;
    m.add_class::<PyControl>()?;
    m.add_class::<PyPathType>()?;
    m.add_class::<PySteeringPath>()?;
    m.add_function(wrap_pyfunction!(supported_path_types, m)?)?;
    Ok(())
}