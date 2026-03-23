use approx::assert_relative_eq;
use steering_functions::dubins_state_space::DubinsStateSpace;
use steering_functions::ekf::Ekf;
use steering_functions::hc_cc_state_space::{
    Cc00DubinsStateSpace, Cc00ReedsSheppStateSpace, Cc0pmDubinsStateSpace, CcDubinsStateSpace,
    Ccpm0DubinsStateSpace, CcpmpmDubinsStateSpace, Hc00ReedsSheppStateSpace,
    Hc0pmReedsSheppStateSpace, HcCcStateSpace, HcReedsSheppStateSpace, Hcpm0ReedsSheppStateSpace,
    HcpmpmReedsSheppStateSpace,
};
use steering_functions::reeds_shepp_state_space::ReedsSheppStateSpace;
use steering_functions::steering_functions::{
    BaseStateSpace, Control, Controller, MeasurementNoise, MotionNoise, State,
};
use steering_functions::utilities::{
    EPSILON, HALF_PI, PI, TWO_PI, array_index_min, end_of_circular_arc, end_of_straight_line,
    get_epsilon, global_frame_change, local_frame_change, point_distance, pointer_array_init,
    polar, twopify,
};

#[test]
fn test_all_core_modules_are_usable() {
    let start = State {
        x: 0.0,
        y: 0.0,
        theta: 0.0,
        kappa: 0.0,
        d: 1,
    };
    let goal = State {
        x: 2.0,
        y: 0.0,
        theta: 0.0,
        kappa: 0.0,
        d: 1,
    };

    let dubins = DubinsStateSpace::new(1.0, 0.1);
    let mut controls = Vec::new();
    let path = dubins.get_path(&start, &goal, &mut controls);
    assert_eq!(path.len(), 2);
    assert_eq!(controls.len(), 1);

    let rs = ReedsSheppStateSpace::new(1.0, 0.1);
    let _ = rs.get_all_controls(&start, &goal);

    let hc_spaces: Vec<Box<dyn HcCcStateSpace>> = vec![
        Box::new(Cc00DubinsStateSpace),
        Box::new(Cc0pmDubinsStateSpace),
        Box::new(Ccpm0DubinsStateSpace),
        Box::new(CcpmpmDubinsStateSpace),
        Box::new(CcDubinsStateSpace),
        Box::new(HcReedsSheppStateSpace),
        Box::new(Hc00ReedsSheppStateSpace),
        Box::new(Hc0pmReedsSheppStateSpace),
        Box::new(Hcpm0ReedsSheppStateSpace),
        Box::new(HcpmpmReedsSheppStateSpace),
        Box::new(Cc00ReedsSheppStateSpace),
    ];
    for space in hc_spaces {
        let c = space.get_controls(&start, &goal);
        let s = space.integrate(&start, &c, 0.1);
        assert!(!s.is_empty());
    }
}

#[test]
fn test_utilities_api_surface() {
    assert_relative_eq!(point_distance(0.0, 0.0, 3.0, 4.0), 5.0, epsilon = 1e-12);
    assert_relative_eq!(twopify(-0.1), TWO_PI - 0.1, epsilon = 1e-12);
    assert!(get_epsilon() >= EPSILON);
    let (r, th) = polar(0.0, 1.0);
    assert_relative_eq!(r, 1.0, epsilon = 1e-12);
    assert_relative_eq!(th, HALF_PI, epsilon = 1e-12);

    let (x2, y2) = end_of_straight_line(1.0, 2.0, 0.0, 1.0, 3.0);
    assert_relative_eq!(x2, 4.0, epsilon = 1e-12);
    assert_relative_eq!(y2, 2.0, epsilon = 1e-12);

    let (_x, _y, t) = end_of_circular_arc(0.0, 0.0, 0.0, 1.0, 1.0, PI / 2.0);
    assert_relative_eq!(t, PI / 2.0, epsilon = 1e-12);

    let (gx, gy) = global_frame_change(1.0, 2.0, 0.0, 3.0, 4.0);
    assert_relative_eq!(gx, 4.0, epsilon = 1e-12);
    assert_relative_eq!(gy, 6.0, epsilon = 1e-12);
    let (lx, ly) = local_frame_change(1.0, 2.0, 0.0, gx, gy);
    assert_relative_eq!(lx, 3.0, epsilon = 1e-12);
    assert_relative_eq!(ly, 4.0, epsilon = 1e-12);

    let idx = array_index_min(&[5.0, 1.0, 2.0]);
    assert_eq!(idx, 1);
    let pointers: Vec<Option<Control>> = pointer_array_init(3);
    assert_eq!(pointers, vec![None, None, None]);
}

#[test]
fn test_ekf_api_surface() {
    let mut ekf = Ekf::default();
    ekf.set_parameters(
        MotionNoise {
            alpha1: 1.0,
            alpha2: 2.0,
            alpha3: 3.0,
            alpha4: 4.0,
        },
        MeasurementNoise {
            std_x: 0.1,
            std_y: 0.2,
            std_theta: 0.3,
        },
        Controller {
            k1: 10.0,
            k2: 20.0,
            k3: 30.0,
        },
    );
    let s = State::default();
    let (_fx, _fu) = ekf.get_motion_jacobi(&s, 0.0, 0.0, 0.0);
    let _hx = ekf.get_observation_jacobi(&s);
    let _q = ekf.get_motion_covariance(0.0, 0.0, 0.0);
    let r = ekf.get_observation_covariance();
    assert_relative_eq!(r.0[0][0], 0.01, epsilon = 1e-12);
    assert_relative_eq!(r.0[1][1], 0.04, epsilon = 1e-12);
    assert_relative_eq!(r.0[2][2], 0.09, epsilon = 1e-12);
    let k = ekf.get_controller_gain();
    assert_relative_eq!(k.0[0][0], 10.0, epsilon = 1e-12);
    assert_relative_eq!(k.0[1][1], 20.0, epsilon = 1e-12);
    assert_relative_eq!(k.0[1][2], 30.0, epsilon = 1e-12);
}
