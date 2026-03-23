#pragma once

#include <cmath>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <sstream>

#include "common/json_if_exist.h"
#include "geometry/geometry_types.h"

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
        /** \Driving direction {-1,0,1} */
        double d{0};
        double s{0};

        double vel{0};
        double acc{0};
        double time{0};

        double fork_y{0};

        State() = default;
        State(double _x, double _y, double _theta, double _kappa = 0)
            : x(_x), y(_y), theta(_theta), kappa(_kappa)
        {
        }

        /**
         * @brief Convert State to Pose2D
         * @return The corresponding Pose2D
         */
        Pose2D toPose2D() const { return Pose2D(x, y, theta); }

        /**
         * @brief Convert from Pose2D to State
         * @param pose The pose to convert from
         * @return The corresponding State
         */
        static State fromPose2D(const Pose2D& pose)
        {
            State state;
            state.x     = pose.x;
            state.y     = pose.y;
            state.theta = pose.theta;
            return state;
        }

        // Add equality operator for State comparison with tolerance
        bool operator==(const State& other) const
        {
            const double epsilon = 1e-6; // Define a tolerance level
            return (fabs(x - other.x) < epsilon) && (fabs(y - other.y) < epsilon) &&
                   (fabs(theta - other.theta) < epsilon) && (fabs(kappa - other.kappa) < epsilon) &&
                   (fabs(sigma - other.sigma) < epsilon) && (fabs(d - other.d) < epsilon) &&
                   (fabs(s - other.s) < epsilon) && (fabs(vel - other.vel) < epsilon) &&
                   (fabs(acc - other.acc) < epsilon) && (fabs(time - other.time) < epsilon) &&
                   (fabs(fork_y - other.fork_y) < epsilon);
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
        std::string toString() const
        {
            std::stringstream ss;
            ss << *this;
            return ss.str();
        }

        NLOHMANN_DEFINE_TYPE_INTRUSIVE_IF_EXISTS(
            State, x, y, theta, kappa, d, sigma, s, vel, acc, time, fork_y);
    };

} // namespace steering