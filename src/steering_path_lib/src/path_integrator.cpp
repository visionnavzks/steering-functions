#include "steering_path_lib/path_integrator.h"

#include <cassert>
#include <cmath>

#include "steering_functions/utilities/utilities.hpp"

namespace steering
{

    State PathIntegrator::computeNextState(const State&   currentState,
                                           const Control& control,
                                           double         stepSize)
    {
        State        nextState         = {};
        const double curvature         = control.kappa;
        const double curvatureRate     = control.sigma;
        const double movementDirection = sgn(control.delta_s);

        if (std::fabs(curvatureRate) > kEpsilon)
        {
            // Clothoid curve: curvature changes continuously
            end_of_clothoid(currentState.x,
                            currentState.y,
                            currentState.theta,
                            currentState.kappa,
                            curvatureRate,
                            movementDirection,
                            stepSize,
                            &nextState.x,
                            &nextState.y,
                            &nextState.theta,
                            &nextState.kappa);
            nextState.sigma = curvatureRate;
        }
        else
        {
            if (std::fabs(curvature) > kEpsilon)
            {
                // Circular arc: constant curvature
                end_of_circular_arc(currentState.x,
                                    currentState.y,
                                    currentState.theta,
                                    curvature,
                                    movementDirection,
                                    stepSize,
                                    &nextState.x,
                                    &nextState.y,
                                    &nextState.theta);
                nextState.kappa = curvature;
            }
            else
            {
                // Straight line: zero curvature
                end_of_straight_line(currentState.x,
                                     currentState.y,
                                     currentState.theta,
                                     movementDirection,
                                     stepSize,
                                     &nextState.x,
                                     &nextState.y);
                nextState.theta = currentState.theta;
            }
        }

        nextState.d = movementDirection;
        nextState.s = currentState.s + stepSize;
        return nextState;
    }

    State PathIntegrator::computeIntermediateState(const State&                startState,
                                                   const std::vector<Control>& controls,
                                                   double                      normalizedDistance)
    {
        State currentState{startState}, interpolatedState{startState};

        if (controls.empty())
        {
            return interpolatedState;
        }

        // Calculate total path length and target interpolation distance
        double totalPathLength = 0.0;
        for (const auto& control : controls)
        {
            totalPathLength += std::fabs(control.delta_s);
        }

        normalizedDistance          = std::clamp(normalizedDistance, 0.0, 1.0);
        const double targetDistance = normalizedDistance * totalPathLength;

        // Track accumulated distance
        double accumulatedDistance = 0.0;

        for (const auto& control : controls)
        {
            const double segmentLength    = control.delta_s;
            const double absSegmentLength = std::fabs(segmentLength);

            if (targetDistance - accumulatedDistance > absSegmentLength)
            {
                // Move to next segment
                currentState = computeNextState(currentState, control, absSegmentLength);
                accumulatedDistance += absSegmentLength;
            }
            else
            {
                // Interpolate within current segment
                interpolatedState =
                    computeNextState(currentState, control, targetDistance - accumulatedDistance);
                break;
            }
        }
        return interpolatedState;
    }

    std::vector<State> PathIntegrator::generateDiscretizedPath(const State& startState,
                                                               const std::vector<Control>& controls,
                                                               double                      stepSize)
    {
        std::vector<State> discretizedPath;
        if (controls.empty())
        {
            return discretizedPath;
        }

        // Initialize states
        State currentState{}, nextState;

        // Pre-calculate path size for efficient memory allocation
        int estimatedPathSize = 0;
        for (const auto& control : controls)
        {
            const double segmentLength = std::fabs(control.delta_s);
            estimatedPathSize += static_cast<int>(std::ceil(segmentLength / stepSize));
        }
        discretizedPath.reserve(estimatedPathSize + 13); // Add buffer for safety

        // Initialize first state
        currentState.x     = startState.x;
        currentState.y     = startState.y;
        currentState.theta = startState.theta;
        currentState.kappa = controls.front().kappa;
        currentState.sigma = controls.front().sigma;
        currentState.d     = sgn(controls.front().delta_s);
        discretizedPath.push_back(currentState);

        // Process each control segment
        for (const auto& control : controls)
        {
            const double segmentLength    = control.delta_s;
            const double absSegmentLength = std::fabs(segmentLength);

            // Skip negligibly small segments
            if (absSegmentLength < kMinSegmentLength)
            {
                continue;
            }

            // Calculate number of points for this segment
            const int numPoints =
                static_cast<int>(std::max(2.0, std::ceil(absSegmentLength / stepSize)));
            const double adjustedStepSize = absSegmentLength / numPoints;
            assert(adjustedStepSize > 0);

            double accumulatedDistance = 0.0;
            for (int i = 0; i < numPoints; ++i)
            {
                // Calculate step size for this iteration
                accumulatedDistance += adjustedStepSize;
                const double currentStepSize =
                    (accumulatedDistance > absSegmentLength)
                        ? adjustedStepSize - (accumulatedDistance - absSegmentLength)
                        : adjustedStepSize;

                if (accumulatedDistance > absSegmentLength)
                {
                    accumulatedDistance = absSegmentLength;
                }

                nextState = computeNextState(currentState, control, currentStepSize);
                discretizedPath.push_back(nextState);
                currentState = nextState;
            }
        }
        return discretizedPath;
    }

} // namespace steering
