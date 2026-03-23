pub mod utilities;

#[cfg(test)]
mod tests {
    use crate::utilities::{PI, TWO_PI, pify, point_distance, twopify};
    use approx::assert_relative_eq;

    #[test]
    fn test_point_distance() {
        assert_relative_eq!(point_distance(0.0, 0.0, 3.0, 4.0), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn test_twopify_wraps_to_0_2pi() {
        assert_relative_eq!(twopify(-0.25), TWO_PI - 0.25, epsilon = 1e-12);
        assert_relative_eq!(twopify(TWO_PI + 0.5), 0.5, epsilon = 1e-12);
    }

    #[test]
    fn test_pify_wraps_to_minus_pi_pi() {
        assert_relative_eq!(pify(PI + 0.1), -PI + 0.1, epsilon = 1e-12);
        assert_relative_eq!(pify(-PI - 0.1), PI - 0.1, epsilon = 1e-12);
    }
}
