/*********************************************************************
 *  Copyright (c) 2017 Robert Bosch GmbH.
 *  All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 ***********************************************************************/

#include "steering_functions/hc_cc_state_space/hc_cc_state_space.hpp"

#include <cmath>

#include "steering_functions/utilities/utilities.hpp"

using namespace std;

namespace steering
{

    HC_CC_State_Space::HC_CC_State_Space(double kappa, double sigma, double discretization)
        : BaseStateSpace(discretization), kappa_(kappa), sigma_(sigma)
    {
        // assert positive inputs
        assert(kappa > 0.0 && sigma > 0.0 && discretization > 0.0);
        // intermediate configuration after first clothoid
        double length_min = kappa / sigma;
        double x_i, y_i, theta_i;
        if (length_min > get_epsilon())
        {
            double kappa_i;
            end_of_clothoid(0, 0, 0, 0, sigma, 1, length_min, &x_i, &y_i, &theta_i, &kappa_i);
        }
        else
        {
            x_i     = 0;
            y_i     = 0;
            theta_i = 0;
        }
        // radius
        double xc, yc;
        xc            = x_i - sin(theta_i) / kappa;
        yc            = y_i + cos(theta_i) / kappa;
        double radius = point_distance(xc, yc, 0.0, 0.0);
        // mu
        double mu     = atan(fabs(xc / yc));
        double sin_mu = sin(mu);
        double cos_mu = cos(mu);
        // delta_min
        double delta_min = 0.5 * pow(kappa, 2) / sigma;
        // assign
        hc_cc_circle_param_.set_param(kappa, sigma, radius, mu, sin_mu, cos_mu, delta_min);
    }

    void HC_CC_State_Space::set_filter_parameters(const Motion_Noise&      motion_noise,
                                                  const Measurement_Noise& measurement_noise,
                                                  const Controller&        controller)
    {
        ekf_.set_parameters(motion_noise, measurement_noise, controller);
    }

    vector<State> HC_CC_State_Space::get_path(const State&          state1,
                                              const State&          state2,
                                              std::vector<Control>& controls)
    {
        auto un_filtered_controls = get_controls(state1, state2);

        for (auto& control : un_filtered_controls)
        {
            if (fabs(control.delta_s) > 1e-6)
            {
                controls.push_back(control);
            }
        }

        return integrate(state1, controls);
    }

    vector<State_With_Covariance>
    HC_CC_State_Space::get_path_with_covariance(const State_With_Covariance& state1,
                                                const State&                 state2) const
    {
        vector<Control> controls = get_controls(state1.state, state2);
        return integrate_with_covariance(state1, controls);
    }

    vector<State_With_Covariance>
    HC_CC_State_Space::integrate_with_covariance(const State_With_Covariance& state,
                                                 const vector<Control>&       controls) const
    {
        vector<State_With_Covariance> path_with_covariance;
        State_With_Covariance         state_curr, state_pred, state_next;
        // reserve capacity of path
        int n_states(0);
        for (const auto& control : controls)
        {
            double abs_delta_s(fabs(control.delta_s));
            n_states += ceil(abs_delta_s / discretization_);
        }
        path_with_covariance.reserve(n_states + 3);
        // push back first state
        state_curr.state.x     = state.state.x;
        state_curr.state.y     = state.state.y;
        state_curr.state.theta = state.state.theta;
        state_curr.state.kappa = controls.front().kappa;
        state_curr.state.d     = sgn(controls.front().delta_s);
        for (int i = 0; i < 16; i++)
        {
            state_curr.Sigma[i]      = state.Sigma[i];
            state_curr.Lambda[i]     = state.Lambda[i];
            state_curr.covariance[i] = state.covariance[i];
        }
        path_with_covariance.push_back(state_curr);

        for (const auto& control : controls)
        {
            double delta_s(control.delta_s);
            double abs_delta_s(fabs(delta_s));
            double kappa(control.kappa);
            double s_seg(0.0);
            double integration_step(0.0);
            // push_back current state if curvature discontinuity
            if (fabs(kappa - state_curr.state.kappa) > get_epsilon())
            {
                state_curr.state.kappa = kappa;
                state_curr.state.d     = sgn(delta_s);
                path_with_covariance.push_back(state_curr);
            }

            for (int i = 0, n = ceil(abs_delta_s / discretization_); i < n; ++i)
            {
                // get integration step
                s_seg += discretization_;
                if (s_seg > abs_delta_s)
                {
                    integration_step = discretization_ - (s_seg - abs_delta_s);
                    s_seg            = abs_delta_s;
                }
                else
                {
                    integration_step = discretization_;
                }
                // predict
                state_pred.state = BaseStateSpace::integrate_ODE(state_curr.state, control, integration_step);
                ekf_.predict(state_curr, control, integration_step, state_pred);
                // update
                state_next.state = state_pred.state;
                ekf_.update(state_pred, state_next);

                path_with_covariance.push_back(state_next);
                state_curr.state = state_next.state;
                for (int i = 0; i < 16; i++)
                {
                    state_curr.Sigma[i]      = state_next.Sigma[i];
                    state_curr.Lambda[i]     = state_next.Lambda[i];
                    state_curr.covariance[i] = state_next.covariance[i];
                }
            }
        }
        return path_with_covariance;
    }

    std::vector<std::vector<Control>> HC_CC_State_Space::get_all_controls(const State& state1,
                                                                          const State& state2) const
    {
        std::vector<std::vector<Control>> all_controls;
        auto                              control = get_controls(state1, state2);
        all_controls.push_back(control);
        return all_controls;
    }

} // namespace steering
