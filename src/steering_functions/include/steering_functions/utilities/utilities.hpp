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

*  This source code is derived from Continuous Curvature (CC) Steer.
*  Copyright (c) 2016, Thierry Fraichard and Institut national de
*  recherche en informatique et en automatique (Inria), licensed under
*  the BSD license, cf. 3rd-party-licenses.txt file in the root
*  directory of this source tree.
**********************************************************************/

#pragma once

#include <algorithm>
#include <cmath>

#include "fresnel.data"

namespace steering
{

    constexpr double PI              = 3.1415926535897932384;
    constexpr double HALF_PI         = 1.5707963267948966192;
    constexpr double TWO_PI          = 6.2831853071795864770;
    constexpr double SQRT_PI         = 1.7724538509055160273;
    constexpr double SQRT_PI_INV     = 0.56418958354775628695;
    constexpr double SQRT_TWO_PI_INV = 0.39894228040143267794;
    constexpr double epsilon         = 1e-4;

    /** \brief Return value of epsilon */
    inline double get_epsilon() { return epsilon; }

    /** \brief Return sign of a number */
    inline double sgn(double x) { return (x < 0) ? -1.0 : 1.0; }

    /** \brief Cartesian distance between two points */
    inline double point_distance(double x1, double y1, double x2, double y2)
    {
        double dx = x2 - x1, dy = y2 - y1;
        return std::sqrt(dx * dx + dy * dy);
    }

    /** \brief Computation of a point's polar coordinates */
    void polar(double x, double y, double& r, double& theta);

    /** \brief Conversion of arbitrary angle given in [rad] to [0, 2*pi[ */
    double twopify(double alpha);

    /** \brief Conversion of arbitrary angle given in [rad] to [-pi, pi[ */
    double pify(double alpha);

    /** \brief Fresnel integrals: S_f = int_0_s(sin(pi/2 u*u)du), C_f = int_0_s(cos(pi/2 u*u)du)
       approximated with Chebyshev polynomials
        */
    void fresnel(double s, double& S_f, double& C_f);

    /** \brief Computation of the end point on a clothoid */
    void end_of_clothoid(double  x_i,
                         double  y_i,
                         double  theta_i,
                         double  kappa_i,
                         double  sigma,
                         double  direction,
                         double  length,
                         double* x_f,
                         double* y_f,
                         double* theta_f,
                         double* kappa_f);

    /** \brief Computation of the end point on a circular arc */
    void end_of_circular_arc(double  x_i,
                             double  y_i,
                             double  theta_i,
                             double  kappa,
                             double  direction,
                             double  length,
                             double* x_f,
                             double* y_f,
                             double* theta_f);

    /** \brief Computation of the end point on a straight line */
    void end_of_straight_line(double  x_i,
                              double  y_i,
                              double  theta,
                              double  direction,
                              double  length,
                              double* x_f,
                              double* y_f);

    /** \brief Transformation of (local_x, local_y) from local coordinate system to global one */
    void global_frame_change(double  x,
                             double  y,
                             double  theta,
                             double  local_x,
                             double  local_y,
                             double* global_x,
                             double* global_y);

    /** \brief Transformation of (global_x, global_y) from global coordinate system to local one */
    void local_frame_change(double  x,
                            double  y,
                            double  theta,
                            double  global_x,
                            double  global_y,
                            double* local_x,
                            double* local_y);

    /** \brief Find index with minimal value in double array */
    inline int array_index_min(double array[], int size)
    {
        return static_cast<int>(std::min_element(array, array + size) - array);
    }

    /** \brief Initialize an array with a given value */
    inline void double_array_init(double array[], int size, double value)
    {
        std::fill(array, array + size, value);
    }

    /** \brief Initialize an array with nullptr */
    inline void pointer_array_init(void* array[], int size)
    {
        std::fill(array, array + size, nullptr);
    }

} // namespace steering
