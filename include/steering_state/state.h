#pragma once

#include <cmath>
#include <ostream>
#include <sstream>
#include <string>

namespace steering
{

    /** \brief Description of a kinematic car's state */
    struct State
    {
        /** \brief Position in x of the robot */
        double x{0};

        /** \brief Position in y of the robot */
        double y{0};

        /** \brief Orientation of the robot */
        double theta{0};

        /** \brief Curvature at position (x,y) */
        double kappa{0};
        double sigma{0};
        /** \brief Driving direction {-1,0,1} */
        double d{0};
        double s{0};

        double vel{0};
        double acc{0};
        double time{0};

        double fork_y{0};

        State() = default;
        State(double x_value, double y_value, double theta_value, double kappa_value = 0)
            : x(x_value), y(y_value), theta(theta_value), kappa(kappa_value)
        {
        }

        bool operator==(const State& other) const
        {
            constexpr double epsilon = 1e-6;
            return (std::fabs(x - other.x) < epsilon) && (std::fabs(y - other.y) < epsilon) &&
                   (std::fabs(theta - other.theta) < epsilon) &&
                   (std::fabs(kappa - other.kappa) < epsilon) &&
                   (std::fabs(sigma - other.sigma) < epsilon) &&
                   (std::fabs(d - other.d) < epsilon) && (std::fabs(s - other.s) < epsilon) &&
                   (std::fabs(vel - other.vel) < epsilon) &&
                   (std::fabs(acc - other.acc) < epsilon) &&
                   (std::fabs(time - other.time) < epsilon) &&
                   (std::fabs(fork_y - other.fork_y) < epsilon);
        }

        friend std::ostream& operator<<(std::ostream& os, const State& state)
        {
            os << "State("
               << "x: " << state.x << ", "
               << "y: " << state.y << ", "
               << "theta: " << state.theta << ", "
               << "kappa: " << state.kappa << ", "
               << "sigma: " << state.sigma << ", "
               << "d: " << state.d << ", "
               << "s: " << state.s << ", "
               << "vel: " << state.vel << ", "
               << "acc: " << state.acc << ", "
               << "time: " << state.time << ", "
               << "fork_y: " << state.fork_y << ")";
            return os;
        }

        /**
         * @brief Convert state to string representation
         * @return String representation of the state in a human-readable format
         */
        [[nodiscard]] std::string toString() const
        {
            std::stringstream ss;
            ss << *this;
            return ss.str();
        }
    };

} // namespace steering