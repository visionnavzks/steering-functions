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

    std::vector<std::vector<Control>> HC_CC_State_Space::get_all_controls(const State& state1,
                                                                          const State& state2) const
    {
        std::vector<std::vector<Control>> all_controls;
        auto                              control = get_controls(state1, state2);
        all_controls.push_back(control);
        return all_controls;
    }

} // namespace steering
