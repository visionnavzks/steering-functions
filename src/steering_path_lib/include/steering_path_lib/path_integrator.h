#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "steering_functions/steering_functions.hpp"

namespace steering
{

    /**
     * @brief Class for handling path integration operations
     *
     * This class provides static methods for integrating paths from controls,
     * interpolating between states, and performing ODE integration.
     */
    class PathIntegrator
    {
    public:
        PathIntegrator()  = default;
        ~PathIntegrator() = default;

        /**
         * @brief Performs single-step integration of ordinary differential equations
         * @param state Current state
         * @param control Control input
         * @param integration_step Integration step size
         * @return Next state after integration
         */
        static State
        computeNextState(const State& state, const Control& control, double integration_step);

        /**
         * @brief Computes an intermediate state along a path at a given normalized distance
         * @param state Starting state
         * @param controls Vector of controls defining the path
         * @param normalizedDistance Interpolation parameter in range [0,1]
         * @return Interpolated state at the specified distance
         */
        static State computeIntermediateState(const State&                state,
                                              const std::vector<Control>& controls,
                                              double                      normalizedDistance);

        /**
         * @brief Generates a discretized path by integrating a sequence of controls
         * @param state Starting state
         * @param controls Vector of controls defining the path
         * @param stepSize Discretization step size for path sampling
         * @return Vector of states representing the discretized path
         */
        static std::vector<State> generateDiscretizedPath(const State&                state,
                                                          const std::vector<Control>& controls,
                                                          double                      stepSize);

    private:
        static constexpr double kEpsilon          = 1e-10;
        static constexpr double kMinSegmentLength = 1e-6;
    };

} // namespace steering
