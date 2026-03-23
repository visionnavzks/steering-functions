#include "steering_refactor/angle_utils.hpp"

#include <cmath>

namespace steering_refactor
{

double normalize_zero_to_two_pi(double angle)
{
    return angle - kTwoPi * std::floor(angle / kTwoPi);
}

double normalize_minus_pi_to_pi(double angle)
{
    double normalized = std::fmod(angle, kTwoPi);
    if (normalized < -kPi)
    {
        normalized += kTwoPi;
    }
    else if (normalized > kPi)
    {
        normalized -= kTwoPi;
    }
    return normalized;
}

double point_distance(double x1, double y1, double x2, double y2)
{
    const double dx = x2 - x1;
    const double dy = y2 - y1;
    return std::sqrt(dx * dx + dy * dy);
}

} // namespace steering_refactor
