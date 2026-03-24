#pragma once

#include "steering_functions/hc_cc_core/hc_cc_state_space.hpp"/*********************************************************************
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

#ifndef HC_CC_STATE_SPACE_HPP
#define HC_CC_STATE_SPACE_HPP

#include <vector>

#include "steering_functions/hc_cc_state_space/hc_cc_circle.hpp"
#include "steering_functions/steering_functions.hpp"

namespace steering
{

    class HC_CC_State_Space : public BaseStateSpace
    {
    public:
        using BaseStateSpace::get_path;

        /** \brief Constructor */
        HC_CC_State_Space(double kappa, double sigma, double discretization);

        /** \brief Virtual function that returns controls of the shortest path from state1 to state2
         */
        std::vector<Control> get_controls(const State& state1,
                                          const State& state2) const override = 0;

        /** \brief Returns path from state1 to state2.
         *  Overrides base to filter out negligibly small control segments. */
        std::vector<State>
        get_path(const State& state1, const State& state2, std::vector<Control>& controls) override;

        std::vector<std::vector<Control>> get_all_controls(const State& state1,
                                                           const State& state2) const override;

    protected:
        /** \brief Curvature, sharpness of clothoid */
        double kappa_, sigma_;

        /** \brief Parameters of a hc-/cc-circle */
        HC_CC_Circle_Param hc_cc_circle_param_;
    };

} // namespace steering

#endif
