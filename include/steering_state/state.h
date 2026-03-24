#pragma once

#include <cmath>
#include <ostream>
#include <sstream>
#include <string>

#if __has_include(<nlohmann/json.hpp>)
#include <nlohmann/json.hpp>
#define STEERING_STATE_HAS_NLOHMANN_JSON 1
#endif

#if __has_include("common/json_if_exist.h")
#include "common/json_if_exist.h"
#define STEERING_STATE_HAS_JSON_IF_EXIST 1
#endif

#if __has_include("geometry/geometry_types.h")
#include "geometry/geometry_types.h"
#define STEERING_STATE_HAS_POSE2D 1
#endif

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

#if defined(STEERING_STATE_HAS_POSE2D)
        /**
         * @brief Convert State to Pose2D
         * @return The corresponding Pose2D
         */
        [[nodiscard]] Pose2D toPose2D() const { return Pose2D(x, y, theta); }

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
#endif

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

#if defined(STEERING_STATE_HAS_NLOHMANN_JSON) && defined(STEERING_STATE_HAS_JSON_IF_EXIST)
        NLOHMANN_DEFINE_TYPE_INTRUSIVE_IF_EXISTS(
            State, x, y, theta, kappa, d, sigma, s, vel, acc, time, fork_y);
#elif defined(STEERING_STATE_HAS_NLOHMANN_JSON)
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(
            State, x, y, theta, kappa, d, sigma, s, vel, acc, time, fork_y);
#endif
    };

} // namespace steering