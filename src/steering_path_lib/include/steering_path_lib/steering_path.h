#pragma once

#include <memory>
#include <stdexcept>

#include "steering_functions/steering_functions.hpp"

namespace steering
{

    /**
     * @brief Represents different types of path planning algorithms
     */
    enum class PathType
    {
        NONE,
        CC_DUBINS,
        CC00_DUBINS,
        CC0PM_DUBINS,
        CCPM0_DUBINS,
        CCPMPM_DUBINS,
        DUBINS,
        CC00_RS,
        HC_RS,
        HC00_RS,
        HC0PM_RS,
        HCPM0_RS,
        HCPMPM_RS,
        RS
    };

    /**
     * @brief High-level wrapper for steering path computation.
     *
     * Provides a unified interface to compute optimal (shortest) or all possible
     * paths between two states using the specified steering function algorithm.
     */
    class SteeringPath
    {
    public:
        /**
         * @brief Construct a new Steering Path object
         * @param path_type Type of path planning algorithm to use
         * @param kappa_max Maximum curvature (must be > 0)
         * @param sigma_max Maximum rate of change of curvature (must be >= 0)
         * @param discretization Path discretization step (must be > 0)
         * @throws std::invalid_argument if parameters are invalid
         * @throws std::runtime_error if planner creation fails
         */
        SteeringPath(PathType path_type, double kappa_max, double sigma_max, float discretization);

        SteeringPath(const SteeringPath&)            = delete;
        SteeringPath& operator=(const SteeringPath&) = delete;
        SteeringPath(SteeringPath&&) noexcept            = default;
        SteeringPath& operator=(SteeringPath&&) noexcept = default;
        ~SteeringPath() = default;

        // --- Optimal (shortest) path interface ---

        /**
         * @brief Get the control sequence for the shortest path
         * @param start Starting state
         * @param goal Goal state
         * @return Control sequence of the shortest path
         */
        [[nodiscard]] std::vector<Control>
        computeShortestControlSequence(const State& start, const State& goal) const;

        /**
         * @brief Get the discretized shortest path
         * @param start Starting state
         * @param goal Goal state
         * @return Vector of states representing the shortest path
         */
        [[nodiscard]] std::vector<State>
        computeShortestPath(const State& start, const State& goal) const;

        // --- All paths interface ---

        /**
         * @brief Get all possible control sequences between two states
         * @param start Starting state
         * @param goal Goal state
         * @return Vector of all possible control sequences
         */
        [[nodiscard]] std::vector<std::vector<Control>>
        computeAllControlSequences(const State& start, const State& goal) const;

        /**
         * @brief Get all possible discretized paths between two states
         * @param start Starting state
         * @param goal Goal state
         * @return Vector of all discretized paths
         */
        [[nodiscard]] std::vector<std::vector<State>>
        computeAllPaths(const State& start, const State& goal) const;

        // --- Getters ---

        [[nodiscard]] PathType getPathType() const { return path_type_; }
        [[nodiscard]] double   getDiscretization() const { return discretization_; }
        [[nodiscard]] double   getKappaMax() const { return kappa_max_; }
        [[nodiscard]] double   getSigmaMax() const { return sigma_max_; }

    private:
        /** @brief Create appropriate state space planner based on path type */
        [[nodiscard]] std::unique_ptr<BaseStateSpace> createPlanner() const;

        PathType                        path_type_;
        double                          discretization_;
        double                          kappa_max_;
        double                          sigma_max_;
        std::unique_ptr<BaseStateSpace> base_state_space_;
    };

} // namespace steering
