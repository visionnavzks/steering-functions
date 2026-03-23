use crate::steering_functions::{Controller, MeasurementNoise, MotionNoise, State};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix2x2(pub [[f64; 2]; 2]);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix2x3(pub [[f64; 3]; 2]);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix3x2(pub [[f64; 2]; 3]);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix3x3(pub [[f64; 3]; 3]);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix3x4(pub [[f64; 4]; 3]);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix4x3(pub [[f64; 3]; 4]);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix4x4(pub [[f64; 4]; 4]);

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Ekf {
    pub motion_noise: MotionNoise,
    pub measurement_noise: MeasurementNoise,
    pub controller: Controller,
}

impl Ekf {
    pub fn set_parameters(
        &mut self,
        motion_noise: MotionNoise,
        measurement_noise: MeasurementNoise,
        controller: Controller,
    ) {
        self.motion_noise = motion_noise;
        self.measurement_noise = measurement_noise;
        self.controller = controller;
    }

    pub fn get_motion_jacobi(
        &self,
        _state: &State,
        _control_kappa: f64,
        _control_sigma: f64,
        _control_delta_s: f64,
    ) -> (Matrix4x4, Matrix4x3) {
        (
            Matrix4x4([
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]),
            Matrix4x3([[0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]]),
        )
    }

    pub fn get_observation_jacobi(&self, _state: &State) -> Matrix3x4 {
        Matrix3x4([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ])
    }

    pub fn get_motion_covariance(
        &self,
        _control_kappa: f64,
        _control_sigma: f64,
        _control_delta_s: f64,
    ) -> Matrix3x3 {
        Matrix3x3([[0.0; 3], [0.0; 3], [0.0; 3]])
    }

    pub fn get_observation_covariance(&self) -> Matrix3x3 {
        Matrix3x3([
            [
                self.measurement_noise.std_x * self.measurement_noise.std_x,
                0.0,
                0.0,
            ],
            [
                0.0,
                self.measurement_noise.std_y * self.measurement_noise.std_y,
                0.0,
            ],
            [
                0.0,
                0.0,
                self.measurement_noise.std_theta * self.measurement_noise.std_theta,
            ],
        ])
    }

    pub fn get_controller_gain(&self) -> Matrix2x3 {
        Matrix2x3([
            [self.controller.k1, 0.0, 0.0],
            [0.0, self.controller.k2, self.controller.k3],
        ])
    }
}
