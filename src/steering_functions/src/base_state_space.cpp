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

#include "steering_functions/steering_functions.hpp"

#include <cassert>
#include <cmath>

#include "steering_functions/utilities/utilities.hpp"

using namespace std;

namespace steering
{

    State BaseStateSpace::integrate_ODE(const State&   state,
                                        const Control& control,
                                        double         integration_step)
    {
        State  state_next = {};
        double kappa(control.kappa);
        double sigma(control.sigma);
        double d(sgn(control.delta_s));
        if (fabs(sigma) > get_epsilon())
        {
            end_of_clothoid(state.x,
                            state.y,
                            state.theta,
                            state.kappa,
                            sigma,
                            d,
                            integration_step,
                            &state_next.x,
                            &state_next.y,
                            &state_next.theta,
                            &state_next.kappa);
            state_next.sigma = sigma;
        }
        else
        {
            if (fabs(kappa) > get_epsilon())
            {
                end_of_circular_arc(state.x,
                                    state.y,
                                    state.theta,
                                    kappa,
                                    d,
                                    integration_step,
                                    &state_next.x,
                                    &state_next.y,
                                    &state_next.theta);
                state_next.kappa = kappa;
            }
            else
            {
                end_of_straight_line(state.x,
                                     state.y,
                                     state.theta,
                                     d,
                                     integration_step,
                                     &state_next.x,
                                     &state_next.y);
                state_next.theta = state.theta;
            }
        }
        state_next.d = d;
        state_next.s = state.s + integration_step;
        return state_next;
    }

    vector<State> BaseStateSpace::get_path(const State& state1, const State& state2)
    {
        vector<Control> controls;
        return get_path(state1, state2, controls);
    }

    vector<State> BaseStateSpace::get_path(const State&          state1,
                                           const State&          state2,
                                           std::vector<Control>& controls)
    {
        controls = get_controls(state1, state2);
        return integrate(state1, controls);
    }

    State BaseStateSpace::interpolate(const State&           state,
                                      const vector<Control>& controls,
                                      double                 t) const
    {
        State state_curr{state}, state_inter{state};

        if (controls.empty())
        {
            return state_inter;
        }

        // Calculate total path length
        double s_path = 0.0;
        for (const auto& control : controls)
        {
            s_path += fabs(control.delta_s);
        }

        t                    = min(max(t, 0.0), 1.0);
        double s_inter       = t * s_path;
        double s_accumulated = 0.0;

        for (const auto& control : controls)
        {
            double abs_delta_s = fabs(control.delta_s);

            if (s_inter - s_accumulated > abs_delta_s)
            {
                state_curr = integrate_ODE(state_curr, control, abs_delta_s);
                s_accumulated += abs_delta_s;
            }
            else
            {
                state_inter = integrate_ODE(state_curr, control, s_inter - s_accumulated);
                break;
            }
        }
        return state_inter;
    }

    vector<State> BaseStateSpace::integrate(const State&           state,
                                            const vector<Control>& controls,
                                            double                 discretization) const
    {
        vector<State> path;
        if (controls.empty())
        {
            return path;
        }

        State state_curr{}, state_next;
        // reserve capacity of path
        int n_states = 0;
        for (const auto& control : controls)
        {
            double abs_delta_s = fabs(control.delta_s);
            n_states += static_cast<int>(ceil(abs_delta_s / discretization));
        }
        path.reserve(n_states + 13);

        // initialize first state
        state_curr.x     = state.x;
        state_curr.y     = state.y;
        state_curr.theta = state.theta;
        state_curr.kappa = controls.front().kappa;
        state_curr.sigma = controls.front().sigma;
        state_curr.d     = sgn(controls.front().delta_s);
        path.push_back(state_curr);

        for (const auto& control : controls)
        {
            double abs_delta_s     = fabs(control.delta_s);
            double s_seg           = 0.0;
            double integration_step = 0.0;

            int    n             = static_cast<int>(max(2.0, ceil(abs_delta_s / discretization)));
            double step_size     = abs_delta_s / n;
            assert(step_size > 0);

            for (int i = 0; i < n; ++i)
            {
                s_seg += step_size;
                if (s_seg > abs_delta_s)
                {
                    integration_step = step_size - (s_seg - abs_delta_s);
                    s_seg            = abs_delta_s;
                }
                else
                {
                    integration_step = step_size;
                }
                state_next = integrate_ODE(state_curr, control, integration_step);
                path.push_back(state_next);
                state_curr = state_next;
            }
        }
        return path;
    }

    vector<State> BaseStateSpace::integrate(const State&           state,
                                            const vector<Control>& controls) const
    {
        return integrate(state, controls, discretization_);
    }

    vector<vector<State>> BaseStateSpace::get_all_paths(const State& state1, const State& state2)
    {
        auto                  all_controls = get_all_controls(state1, state2);
        vector<vector<State>> all_paths;
        all_paths.reserve(all_controls.size());

        for (const auto& controls : all_controls)
        {
            all_paths.push_back(integrate(state1, controls));
        }
        return all_paths;
    }

} // namespace steering
