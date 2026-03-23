use rust_steering_functions::*;

fn main() {
    let state1 = State::new(0.0, 0.0, 0.0, 0.0);
    let state2 = State::new(10.0, 10.0, std::f64::consts::PI / 2.0, 0.0);

    let kappa = 1.0;
    let discretization = 0.5;

    let dubins_path =
        SteeringFunctions::get_dubins_path(&state1, &state2, kappa, discretization, true);

    println!("Dubins Path Generated with {} states", dubins_path.len());
    if let Some(last) = dubins_path.last() {
        println!(
            "  Target: x={}, y={}, theta={}",
            state2.x, state2.y, state2.theta
        );
        println!(
            "  Actual: x={:.4}, y={:.4}, theta={:.4}",
            last.x, last.y, last.theta
        );
    }

    println!("--------------------------------------------------");

    let rs_path = SteeringFunctions::get_reeds_shepp_path(&state1, &state2, kappa, discretization);

    println!("Reeds-Shepp Path Generated with {} states", rs_path.len());
    if let Some(last) = rs_path.last() {
        println!(
            "  Target: x={}, y={}, theta={}",
            state2.x, state2.y, state2.theta
        );
        println!(
            "  Actual: x={:.4}, y={:.4}, theta={:.4}",
            last.x, last.y, last.theta
        );
    }
}
