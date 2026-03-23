#ifndef STEERING_REFACTOR_ANGLE_UTILS_HPP
#define STEERING_REFACTOR_ANGLE_UTILS_HPP

namespace steering_refactor
{

constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 6.28318530717958647692;

// Normalize any angle to [0, 2*pi)
double normalize_zero_to_two_pi(double angle);

// Normalize any angle to [-pi, pi]
double normalize_minus_pi_to_pi(double angle);

// Cartesian distance between two 2D points
double point_distance(double x1, double y1, double x2, double y2);

} // namespace steering_refactor

#endif
