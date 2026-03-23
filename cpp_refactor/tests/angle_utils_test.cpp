#include "steering_refactor/angle_utils.hpp"

#include <cmath>
#include <iostream>

namespace
{

bool almost_equal(double a, double b, double eps = 1e-12)
{
    return std::fabs(a - b) <= eps;
}

} // namespace

int main()
{
    using namespace steering_refactor;

    if (!almost_equal(normalize_zero_to_two_pi(0.0), 0.0))
    {
        std::cerr << "normalize_zero_to_two_pi(0.0) failed" << std::endl;
        return 1;
    }

    if (!almost_equal(normalize_zero_to_two_pi(-kPi), kPi))
    {
        std::cerr << "normalize_zero_to_two_pi(-pi) failed" << std::endl;
        return 1;
    }

    if (!almost_equal(normalize_minus_pi_to_pi(3.0 * kPi), kPi))
    {
        std::cerr << "normalize_minus_pi_to_pi(3*pi) failed" << std::endl;
        return 1;
    }

    if (!almost_equal(normalize_minus_pi_to_pi(-3.0 * kPi), -kPi))
    {
        std::cerr << "normalize_minus_pi_to_pi(-3*pi) failed" << std::endl;
        return 1;
    }

    if (!almost_equal(point_distance(0.0, 0.0, 3.0, 4.0), 5.0))
    {
        std::cerr << "point_distance(3-4-5) failed" << std::endl;
        return 1;
    }

    if (!almost_equal(point_distance(-1.0, -1.0, 2.0, 3.0), point_distance(2.0, 3.0, -1.0, -1.0)))
    {
        std::cerr << "point_distance symmetry failed" << std::endl;
        return 1;
    }

    if (!almost_equal(point_distance(1.25, -8.0, 1.25, -8.0), 0.0))
    {
        std::cerr << "point_distance zero-length failed" << std::endl;
        return 1;
    }

    std::cout << "All angle_utils tests passed." << std::endl;
    return 0;
}
