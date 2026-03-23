#pragma once

#include <memory>
#include <stdexcept>

// Core dependencies
#include "steering_functions/steering_functions.hpp"

namespace steering
{

    /**
     * @brief Represents different types of path planning algorithms
     */
    enum class PathType
    {
        NONE, // No specific path type
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
     * @brief Class representing a path planning solution with various steering functions
     */
    class SteeringPath
    {
    public:
        /**
         * @brief Construct a new Steering Path object
         * @param path_type Type of path planning algorithm to use
         * @param kappa_max Maximum curvature
         * @param sigma_max Maximum rate of change of curvature
         * @param discretization Path discretization step
         * @throws std::invalid_argument if parameters are invalid
         */
        SteeringPath(PathType path_type, double kappa_max, double sigma_max, float discretization);

        // Disable copy operations due to unique_ptr
        SteeringPath(const SteeringPath&)            = delete;
        SteeringPath& operator=(const SteeringPath&) = delete;

        // Enable move operations
        SteeringPath(SteeringPath&&) noexcept            = default;
        SteeringPath& operator=(SteeringPath&&) noexcept = default;

        ~SteeringPath() = default;

        /**
         * @brief Get all possible controls between start and goal states
         * @param state_start Starting state
         * @param state_goal Goal state
         * @param use_shortest_path Whether to only compute shortest path
         * @return Vector of control sequences
         */
        [[nodiscard]] std::vector<std::vector<Control>>
        computeControlSequences(const State& state_start,
                                const State& state_goal,
                                bool         use_shortest_path = false) const;

        /**
         * @brief Compute all possible discretized paths between start and goal states
         * @param state_start Starting state
         * @param state_goal Goal state
         * @param use_shortest_path Whether to only compute shortest path
         * @return Vector of discretized paths, where each path is a vector of states
         */
        [[nodiscard]] std::vector<std::vector<State>>
        computeAllPaths(const State& state_start,
                        const State& state_goal,
                        bool         use_shortest_path = false) const;

        // Getters
        [[nodiscard]] PathType     getPathType() const { return path_type_; }
        [[nodiscard]] double       getDiscretization() const { return discretization_; }
        [[nodiscard]] const State& getStartState() const { return state_start_; }
        [[nodiscard]] const State& getGoalState() const { return state_goal_; }
        [[nodiscard]] double       getKappaMax() const { return kappa_max_; }
        [[nodiscard]] double       getSigmaMax() const { return sigma_max_; }

    private:
        /**
         * @brief Create appropriate state space planner based on path type
         * @return Unique pointer to base state space
         */
        [[nodiscard]] std::unique_ptr<BaseStateSpace> createPlanner() const;

        PathType                        path_type_;
        double                          discretization_;
        State                           state_start_{};
        State                           state_goal_{};
        double                          kappa_max_;
        double                          sigma_max_;
        std::unique_ptr<BaseStateSpace> base_state_space_;
    };

} // namespace steering
