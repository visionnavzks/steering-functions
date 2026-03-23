#include "steering_path_lib/steering_path.h"

#include "steering_path_lib/path_integrator.h"
#include "steering_functions/all_in_one.hpp"
namespace steering
{

    SteeringPath::SteeringPath(PathType path_type,
                               double   kappa_max,
                               double   sigma_max,
                               float    discretization)
        : path_type_(path_type)
        , discretization_(discretization)
        , kappa_max_(kappa_max)
        , sigma_max_(sigma_max)
    {
        if (kappa_max <= 0 || sigma_max < 0 || discretization <= 0)
        {
            throw std::invalid_argument("Invalid parameters in SteeringPath constructor");
        }

        base_state_space_ = createPlanner();
        if (!base_state_space_)
        {
            throw std::runtime_error("Failed to create path planner");
        }
    }

    std::unique_ptr<BaseStateSpace> SteeringPath::createPlanner() const
    {
        switch (path_type_)
        {
        case PathType::CC_DUBINS:
            return std::make_unique<CC_Dubins_State_Space>(
                kappa_max_, sigma_max_, discretization_, true);
        case PathType::CC00_DUBINS:
            return std::make_unique<CC00_Dubins_State_Space>(
                kappa_max_, sigma_max_, discretization_, true);
        case PathType::CC0PM_DUBINS:
            return std::make_unique<CC0pm_Dubins_State_Space>(
                kappa_max_, sigma_max_, discretization_, true);
        case PathType::CCPM0_DUBINS:
            return std::make_unique<CCpm0_Dubins_State_Space>(
                kappa_max_, sigma_max_, discretization_, true);
        case PathType::CCPMPM_DUBINS:
            return std::make_unique<CCpmpm_Dubins_State_Space>(
                kappa_max_, sigma_max_, discretization_, true);
        case PathType::DUBINS:
            return std::make_unique<Dubins_State_Space>(kappa_max_, discretization_, true);
        case PathType::CC00_RS:
            return std::make_unique<CC00_Reeds_Shepp_State_Space>(
                kappa_max_, sigma_max_, discretization_);
        case PathType::HC_RS:
            return std::make_unique<HC_Reeds_Shepp_State_Space>(
                kappa_max_, sigma_max_, discretization_);
        case PathType::HC00_RS:
            return std::make_unique<HC00_Reeds_Shepp_State_Space>(
                kappa_max_, sigma_max_, discretization_);
        case PathType::HC0PM_RS:
            return std::make_unique<HC0pm_Reeds_Shepp_State_Space>(
                kappa_max_, sigma_max_, discretization_);
        case PathType::HCPM0_RS:
            return std::make_unique<HCpm0_Reeds_Shepp_State_Space>(
                kappa_max_, sigma_max_, discretization_);
        case PathType::HCPMPM_RS:
            return std::make_unique<HCpmpm_Reeds_Shepp_State_Space>(
                kappa_max_, sigma_max_, discretization_);
        case PathType::RS:
            return std::make_unique<Reeds_Shepp_State_Space>(kappa_max_, discretization_);
        default:
            return std::make_unique<Reeds_Shepp_State_Space>(kappa_max_, discretization_);
        }
    }

    std::vector<std::vector<Control>> SteeringPath::computeControlSequences(
        const State& state_start, const State& state_goal, bool use_shortest_path) const
    {
        if (!base_state_space_)
        {
            throw std::runtime_error("State space not initialized");
        }

        if (use_shortest_path)
        {
            std::vector<Control> controls;
            // For shortest path, compute single optimal trajectory
            base_state_space_->get_path(state_start, state_goal, controls);
            return {controls};
        }

        // Compute all possible control sequences between states
        return base_state_space_->get_all_controls(state_start, state_goal);
    }

    std::vector<std::vector<State>> SteeringPath::computeAllPaths(const State& state_start,
                                                                  const State& state_goal,
                                                                  bool use_shortest_path) const
    {
        if (!base_state_space_)
        {
            throw std::runtime_error("State space not initialized");
        }

        // Get all possible control sequences
        auto control_sequences =
            computeControlSequences(state_start, state_goal, use_shortest_path);

        // Generate discretized paths for each control sequence
        std::vector<std::vector<State>> all_paths;
        all_paths.reserve(control_sequences.size());

        for (const auto& controls : control_sequences)
        {
            all_paths.push_back(
                PathIntegrator::generateDiscretizedPath(state_start, controls, discretization_));
        }

        return all_paths;
    }

} // namespace steering
