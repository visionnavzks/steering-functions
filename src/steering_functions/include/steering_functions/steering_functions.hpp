/*********************************************************************
 *  Copyright (c) 2017 - for information on the respective copyright
 *  owner see the NOTICE file and/or the repository
 *
 *      https://github.com/hbanzhaf/steering_functions.git
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 *  implied. See the License for the specific language governing
 *  permissions and limitations under the License.
 ***********************************************************************/

#ifndef STEERING_FUNCTIONS_HPP
#define STEERING_FUNCTIONS_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "steering_state/state.h"

namespace steering
{

    /** \brief Description of a kinematic car's state with covariance */
    struct State_With_Covariance
    {
        /** \brief Expected state of the robot */
        State state;

        /** \brief Covariance of the state estimation due to motion and measurement noise */
        double Sigma[16] = {0.0};

        /** \brief Covariance of the state estimate due to the absence of measurements */
        double Lambda[16] = {0.0};

        /** \brief Covariance of the state given by Sigma + Lambda: (x_x      x_y      x_theta
           x_kappa y_x      y_y      y_theta      y_kappa theta_x  theta_y  theta_theta  theta_kappa
                                                                     kappa_x  kappa_y  kappa_theta
           kappa_kappa) */
        double covariance[16] = {0.0};
    };

    /** \brief Description of a path segment with its corresponding control inputs */
    struct Control
    {
        /** \brief Signed arc length of a segment */
        double delta_s;

        /** \brief Curvature at the beginning of a segment */
        double kappa;

        /** \brief Sharpness (derivative of curvature with respect to arc length) of a segment */
        double sigma;

        std::string to_str() const noexcept
        {
            return "Control Segment [delta_s: " + std::to_string(delta_s) +
                   ", curvature: " + std::to_string(kappa) +
                   ", sharpness: " + std::to_string(sigma) + "]";
        }
    };

    /** \brief Parameters of the motion noise model according to the book:
        Probabilistic Robotics, S. Thrun and others, MIT Press, 2006, p. 127-128 and p.204-206. */
    struct Motion_Noise
    {
        /** \brief Variance in longitudinal direction: alpha1*delta_s*delta_s + alpha2*kappa*kappa
         */
        double alpha1;
        double alpha2;

        /** \brief Variance in lateral direction: alpha3*delta_s*delta_s + alpha4*kappa*kappa */
        double alpha3;
        double alpha4;
    };

    /** \brief Parameters of the measurement noise */
    struct Measurement_Noise
    {
        /** \brief Standard deviation of localization in x */
        double std_x;

        /** \brief Standard deviation of localization in y */
        double std_y;

        /** \brief Standard deviation of localization in theta */
        double std_theta;
    };

    /** \brief Parameters of the feedback controller */
    struct Controller
    {
        /** \brief Weight on longitudinal error */
        double k1;

        /** \brief Weight on lateral error */
        double k2;

        /** \brief Weight on heading error */
        double k3;
    };

    /** \brief Base class for all steering state spaces.
     *
     *  Provides common path integration, interpolation, and ODE integration.
     *  Subclasses must implement get_controls() and get_all_controls().
     */
    struct BaseStateSpace
    {
    public:
        virtual ~BaseStateSpace() = default;

        /**
         * @brief Single-step ODE integration for state propagation.
         *
         * Handles clothoid (variable curvature), circular arc (constant curvature),
         * and straight line segments.
         *
         * @param state Current state
         * @param control Control input for this segment
         * @param integration_step Step size for integration
         * @return Next state after integration
         */
        static State
        integrate_ODE(const State& state, const Control& control, double integration_step);

        /**
         * @brief Get controls for the shortest path from state1 to state2
         * @param state1 Starting state
         * @param state2 Goal state
         * @return Vector of controls representing the shortest path
         */
        virtual std::vector<Control> get_controls(const State& state1,
                                                  const State& state2) const = 0;

        /**
         * @brief Get path between two states (convenience overload).
         *        Delegates to the 3-parameter version.
         * @param state1 Starting state
         * @param state2 Goal state
         * @return Vector of states representing the path
         */
        std::vector<State> get_path(const State& state1, const State& state2);

        /**
         * @brief Get path between two states with controls output.
         *        Default: calls get_controls() then integrate().
         * @param state1 Starting state
         * @param state2 Goal state
         * @param controls Output vector to store control sequence
         * @return Vector of states representing the path
         */
        virtual std::vector<State>
        get_path(const State& state1, const State& state2, std::vector<Control>& controls);

        /**
         * @brief Interpolate state along a path at given normalized parameter.
         *        Default implementation uses integrate_ODE.
         * @param state Starting state
         * @param controls Control sequence defining the path
         * @param t Interpolation parameter in [0,1]
         * @return Interpolated state
         */
        virtual State
        interpolate(const State& state, const std::vector<Control>& controls, double t) const;

        /**
         * @brief Integrate controls to generate a discretized path.
         *        Default implementation uses integrate_ODE.
         * @param state Starting state
         * @param controls Control sequence to integrate
         * @param discretization Step size for integration
         * @return Vector of states representing the discretized path
         */
        virtual std::vector<State> integrate(const State&                state,
                                             const std::vector<Control>& controls,
                                             double                      discretization) const;

        /**
         * @brief Integrate controls using the stored discretization step.
         * @param state Starting state
         * @param controls Control sequence to integrate
         * @return Vector of states representing the discretized path
         */
        std::vector<State> integrate(const State&                state,
                                     const std::vector<Control>& controls) const;

        /**
         * @brief Get all possible control sequences between two states
         * @param state1 Starting state
         * @param state2 Goal state
         * @return Vector of all possible control sequences
         */
        virtual std::vector<std::vector<Control>> get_all_controls(const State& state1,
                                                                   const State& state2) const = 0;

    protected:
        explicit BaseStateSpace(double discretization = 0.1) : discretization_(discretization) {}
        BaseStateSpace(const BaseStateSpace&)            = default;
        BaseStateSpace& operator=(const BaseStateSpace&) = default;

        /** \brief Discretization step size for path integration */
        double discretization_;
    };

} // namespace steering

#endif
