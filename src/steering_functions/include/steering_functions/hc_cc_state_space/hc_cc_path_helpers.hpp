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
**********************************************************************/

#ifndef HC_CC_PATH_HELPERS_HPP
#define HC_CC_PATH_HELPERS_HPP

#include <cmath>

#include "steering_functions/hc_cc_state_space/configuration.hpp"
#include "steering_functions/hc_cc_state_space/hc_cc_circle.hpp"
#include "steering_functions/utilities/utilities.hpp"

namespace steering
{
namespace path_helpers
{

    // ===== Existence checks =====

    inline bool TT_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                          double radius)
    {
        if (c1.left == c2.left)
        {
            return false;
        }
        if (c1.forward == c2.forward)
        {
            return false;
        }
        return fabs(distance - 2 * radius) < get_epsilon();
    }

    inline bool TiST_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                             double radius)
    {
        if (c1.left == c2.left)
        {
            return false;
        }
        if (c1.forward == c2.forward)
        {
            return false;
        }
        return (distance >= 2 * radius);
    }

    inline bool TeST_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                             double radius_sin_mu)
    {
        if (c1.left != c2.left)
        {
            return false;
        }
        if (c1.forward == c2.forward)
        {
            return false;
        }
        return (distance >= 2 * radius_sin_mu);
    }

    inline bool TST_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                            double radius, double radius_sin_mu)
    {
        return TiST_exists(distance, c1, c2, radius) ||
               TeST_exists(distance, c1, c2, radius_sin_mu);
    }

    inline bool TTT_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                            double radius)
    {
        if (c1.left != c2.left)
        {
            return false;
        }
        if (c1.forward == c2.forward)
        {
            return false;
        }
        return distance <= 4 * radius;
    }

    // ===== Tangent circle computations =====

    inline void TT_tangent_circles(const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                                    double mu, Configuration** q)
    {
        double x     = (c1.xc + c2.xc) / 2;
        double y     = (c1.yc + c2.yc) / 2;
        double angle = atan2(c2.yc - c1.yc, c2.xc - c1.xc);
        double theta;
        if (c1.left)
        {
            if (c1.forward)
            {
                theta = angle + HALF_PI - mu;
            }
            else
            {
                theta = angle + HALF_PI + mu;
            }
        }
        else
        {
            if (c1.forward)
            {
                theta = angle - HALF_PI + mu;
            }
            else
            {
                theta = angle - HALF_PI - mu;
            }
        }
        *q = new Configuration(x, y, theta, 0);
    }

    inline void TiST_tangent_circles(const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                                      double radius, double sin_mu, double cos_mu,
                                      Configuration** q1, Configuration** q2)
    {
        double distance = center_distance(c1, c2);
        double angle    = atan2(c2.yc - c1.yc, c2.xc - c1.xc);
        double alpha    = asin(2 * radius * cos_mu / distance);
        double delta_x  = radius * sin_mu;
        double delta_y  = radius * cos_mu;
        double x, y, theta;
        if (c1.left && c1.forward)
        {
            theta = angle + alpha;
            global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
            *q1 = new Configuration(x, y, theta, 0);
            global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
            *q2 = new Configuration(x, y, theta, 0);
        }
        if (c1.left && !c1.forward)
        {
            theta = angle - alpha;
            global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
            *q1 = new Configuration(x, y, theta + PI, 0);
            global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
            *q2 = new Configuration(x, y, theta + PI, 0);
        }
        if (!c1.left && c1.forward)
        {
            theta = angle - alpha;
            global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
            *q1 = new Configuration(x, y, theta, 0);
            global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
            *q2 = new Configuration(x, y, theta, 0);
        }
        if (!c1.left && !c1.forward)
        {
            theta = angle + alpha;
            global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
            *q1 = new Configuration(x, y, theta + PI, 0);
            global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
            *q2 = new Configuration(x, y, theta + PI, 0);
        }
    }

    inline void TeST_tangent_circles(const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                                      double radius, double sin_mu, double cos_mu,
                                      Configuration** q1, Configuration** q2)
    {
        double delta_x = radius * sin_mu;
        double delta_y = radius * cos_mu;
        double theta   = atan2(c2.yc - c1.yc, c2.xc - c1.xc);
        double x, y;
        if (c1.left && c1.forward)
        {
            global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
            *q1 = new Configuration(x, y, theta, 0);
            global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
            *q2 = new Configuration(x, y, theta, 0);
        }
        if (c1.left && !c1.forward)
        {
            global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
            *q1 = new Configuration(x, y, theta + PI, 0);
            global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
            *q2 = new Configuration(x, y, theta + PI, 0);
        }
        if (!c1.left && c1.forward)
        {
            global_frame_change(c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
            *q1 = new Configuration(x, y, theta, 0);
            global_frame_change(c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
            *q2 = new Configuration(x, y, theta, 0);
        }
        if (!c1.left && !c1.forward)
        {
            global_frame_change(c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
            *q1 = new Configuration(x, y, theta + PI, 0);
            global_frame_change(c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
            *q2 = new Configuration(x, y, theta + PI, 0);
        }
    }

    // ===== RS-specific existence checks =====
    // These check direction flags and compare distance against a threshold.
    // The threshold formula differs between CC and HC variants, so the
    // caller computes the specific threshold and passes it in.

    inline bool TcT_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                            double two_radius_threshold)
    {
        if (c1.left == c2.left)
        {
            return false;
        }
        if (c1.forward != c2.forward)
        {
            return false;
        }
        return fabs(distance - two_radius_threshold) < get_epsilon();
    }

    inline bool TcTcT_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                              double upper_bound)
    {
        if (c1.left == c2.left)
        {
            return false;
        }
        if (c1.forward == c2.forward)
        {
            return false;
        }
        return distance <= upper_bound;
    }

    inline bool TcTT_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                             double upper_bound, double lower_bound)
    {
        if (c1.left == c2.left)
        {
            return false;
        }
        if (c1.forward != c2.forward)
        {
            return false;
        }
        return (distance <= upper_bound) && (distance >= lower_bound);
    }

    inline bool TTcT_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                             double upper_bound, double lower_bound)
    {
        if (c1.left == c2.left)
        {
            return false;
        }
        if (c1.forward != c2.forward)
        {
            return false;
        }
        return (distance <= upper_bound) && (distance >= lower_bound);
    }

    inline bool TTcTT_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                              double upper_bound)
    {
        if (c1.left == c2.left)
        {
            return false;
        }
        if (c1.forward != c2.forward)
        {
            return false;
        }
        return distance <= upper_bound;
    }

    inline bool TcTTcT_exists(double distance, const HC_CC_Circle& c1, const HC_CC_Circle& c2,
                               double upper_bound, double lower_bound)
    {
        if (c1.left == c2.left)
        {
            return false;
        }
        if (c1.forward == c2.forward)
        {
            return false;
        }
        return (distance <= upper_bound) && (distance >= lower_bound);
    }

} // namespace path_helpers
} // namespace steering

#endif // HC_CC_PATH_HELPERS_HPP
