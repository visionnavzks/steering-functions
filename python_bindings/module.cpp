#include <sstream>

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "steering_functions/all_in_one.hpp"
#include "steering_path_lib/steering_path.h"

namespace py = pybind11;

namespace steering
{
    namespace
    {
        Control make_control(double delta_s, double kappa, double sigma)
        {
            Control control{};
            control.delta_s = delta_s;
            control.kappa   = kappa;
            control.sigma   = sigma;
            return control;
        }

        template <typename Planner>
        void bind_common_planner_api(py::class_<Planner>& planner_class)
        {
            planner_class
                .def("get_distance", &Planner::get_distance)
                .def("get_controls", &Planner::get_controls)
                .def("get_path",
                     [](Planner& planner, const State& start, const State& goal)
                     { return planner.get_path(start, goal); })
                .def("get_all_controls", &Planner::get_all_controls)
                .def("get_all_paths",
                     [](Planner& planner, const State& start, const State& goal)
                     { return planner.get_all_paths(start, goal); });
        }
    } // namespace
} // namespace steering

PYBIND11_MODULE(_core, module)
{
    using steering::Control;
    using steering::CC00_Dubins_State_Space;
    using steering::CC00_Reeds_Shepp_State_Space;
    using steering::CC0pm_Dubins_State_Space;
    using steering::CC_Dubins_State_Space;
    using steering::CCpm0_Dubins_State_Space;
    using steering::CCpmpm_Dubins_State_Space;
    using steering::PathType;
    using steering::Dubins_State_Space;
    using steering::HC00_Reeds_Shepp_State_Space;
    using steering::HC0pm_Reeds_Shepp_State_Space;
    using steering::HC_Reeds_Shepp_State_Space;
    using steering::HCpm0_Reeds_Shepp_State_Space;
    using steering::HCpmpm_Reeds_Shepp_State_Space;
    using steering::Reeds_Shepp_State_Space;
    using steering::State;
    using steering::SteeringPath;

    module.doc() = "Python bindings for steering functions";

    py::class_<State>(module, "State")
        .def(py::init<>())
        .def(py::init<double, double, double, double>(),
             py::arg("x"),
             py::arg("y"),
             py::arg("theta"),
             py::arg("kappa") = 0.0)
        .def_readwrite("x", &State::x)
        .def_readwrite("y", &State::y)
        .def_readwrite("theta", &State::theta)
        .def_readwrite("kappa", &State::kappa)
        .def_readwrite("sigma", &State::sigma)
        .def_readwrite("d", &State::d)
        .def_readwrite("s", &State::s)
        .def_readwrite("vel", &State::vel)
        .def_readwrite("acc", &State::acc)
        .def_readwrite("time", &State::time)
        .def_readwrite("fork_y", &State::fork_y)
        .def("to_string", &State::toString)
        .def("__repr__", [](const State& state) { return state.toString(); })
        .def("__str__", [](const State& state) { return state.toString(); })
        .def(py::self == py::self);

    py::class_<Control>(module, "Control")
        .def(py::init(&steering::make_control),
             py::arg("delta_s") = 0.0,
             py::arg("kappa") = 0.0,
             py::arg("sigma") = 0.0)
        .def_readwrite("delta_s", &Control::delta_s)
        .def_readwrite("kappa", &Control::kappa)
        .def_readwrite("sigma", &Control::sigma)
        .def("to_string", &Control::to_str)
        .def("__repr__", [](const Control& control) { return control.to_str(); })
        .def("__str__", [](const Control& control) { return control.to_str(); });

    py::enum_<PathType>(module, "PathType")
        .value("NONE", PathType::NONE)
        .value("CC_DUBINS", PathType::CC_DUBINS)
        .value("CC00_DUBINS", PathType::CC00_DUBINS)
        .value("CC0PM_DUBINS", PathType::CC0PM_DUBINS)
        .value("CCPM0_DUBINS", PathType::CCPM0_DUBINS)
        .value("CCPMPM_DUBINS", PathType::CCPMPM_DUBINS)
        .value("DUBINS", PathType::DUBINS)
        .value("CC00_RS", PathType::CC00_RS)
        .value("HC_RS", PathType::HC_RS)
        .value("HC00_RS", PathType::HC00_RS)
        .value("HC0PM_RS", PathType::HC0PM_RS)
        .value("HCPM0_RS", PathType::HCPM0_RS)
        .value("HCPMPM_RS", PathType::HCPMPM_RS)
        .value("RS", PathType::RS)
        .export_values();

    py::class_<SteeringPath>(module, "SteeringPath")
        .def(py::init<PathType, double, double, float>(),
             py::arg("path_type"),
             py::arg("kappa_max"),
             py::arg("sigma_max"),
             py::arg("discretization"))
        .def("computeShortestControlSequence", &SteeringPath::computeShortestControlSequence)
        .def("computeShortestPath", &SteeringPath::computeShortestPath)
        .def("computeAllControlSequences", &SteeringPath::computeAllControlSequences)
        .def("computeAllPaths", &SteeringPath::computeAllPaths)
        .def_property_readonly("path_type", &SteeringPath::getPathType)
        .def_property_readonly("discretization", &SteeringPath::getDiscretization)
        .def_property_readonly("kappa_max", &SteeringPath::getKappaMax)
        .def_property_readonly("sigma_max", &SteeringPath::getSigmaMax);

    auto dubins_state_space =
        py::class_<Dubins_State_Space>(module, "DubinsStateSpace")
            .def(py::init<double, double, bool>(),
                 py::arg("kappa"),
                 py::arg("discretization") = 0.1,
                 py::arg("forwards") = true)
            .def("get_controls_reverse", &Dubins_State_Space::get_controls_reverse);
    steering::bind_common_planner_api(dubins_state_space);

    auto reeds_shepp_state_space =
        py::class_<Reeds_Shepp_State_Space>(module, "ReedsSheppStateSpace")
            .def(py::init<double, double>(),
                 py::arg("kappa"),
                 py::arg("discretization") = 0.1);
    steering::bind_common_planner_api(reeds_shepp_state_space);

    auto cc_dubins_state_space =
        py::class_<CC_Dubins_State_Space>(module, "CCDubinsStateSpace")
            .def(py::init<double, double, double, bool>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1,
                 py::arg("forwards") = true);
    steering::bind_common_planner_api(cc_dubins_state_space);

    auto cc00_dubins_state_space =
        py::class_<CC00_Dubins_State_Space>(module, "CC00DubinsStateSpace")
            .def(py::init<double, double, double, bool>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1,
                 py::arg("forwards") = true);
    steering::bind_common_planner_api(cc00_dubins_state_space);

    auto cc0pm_dubins_state_space =
        py::class_<CC0pm_Dubins_State_Space>(module, "CC0pmDubinsStateSpace")
            .def(py::init<double, double, double, bool>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1,
                 py::arg("forwards") = true);
    steering::bind_common_planner_api(cc0pm_dubins_state_space);

    auto ccpm0_dubins_state_space =
        py::class_<CCpm0_Dubins_State_Space>(module, "CCpm0DubinsStateSpace")
            .def(py::init<double, double, double, bool>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1,
                 py::arg("forwards") = true);
    steering::bind_common_planner_api(ccpm0_dubins_state_space);

    auto ccpmpm_dubins_state_space =
        py::class_<CCpmpm_Dubins_State_Space>(module, "CCpmpmDubinsStateSpace")
            .def(py::init<double, double, double, bool>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1,
                 py::arg("forwards") = true);
    steering::bind_common_planner_api(ccpmpm_dubins_state_space);

    auto cc00_reeds_shepp_state_space =
        py::class_<CC00_Reeds_Shepp_State_Space>(module, "CC00ReedsSheppStateSpace")
            .def(py::init<double, double, double>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1);
    steering::bind_common_planner_api(cc00_reeds_shepp_state_space);

    auto hc_reeds_shepp_state_space =
        py::class_<HC_Reeds_Shepp_State_Space>(module, "HCReedsSheppStateSpace")
            .def(py::init<double, double, double>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1);
    steering::bind_common_planner_api(hc_reeds_shepp_state_space);

    auto hc00_reeds_shepp_state_space =
        py::class_<HC00_Reeds_Shepp_State_Space>(module, "HC00ReedsSheppStateSpace")
            .def(py::init<double, double, double>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1);
    steering::bind_common_planner_api(hc00_reeds_shepp_state_space);

    auto hc0pm_reeds_shepp_state_space =
        py::class_<HC0pm_Reeds_Shepp_State_Space>(module, "HC0pmReedsSheppStateSpace")
            .def(py::init<double, double, double>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1);
    steering::bind_common_planner_api(hc0pm_reeds_shepp_state_space);

    auto hcpm0_reeds_shepp_state_space =
        py::class_<HCpm0_Reeds_Shepp_State_Space>(module, "HCpm0ReedsSheppStateSpace")
            .def(py::init<double, double, double>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1);
    steering::bind_common_planner_api(hcpm0_reeds_shepp_state_space);

    auto hcpmpm_reeds_shepp_state_space =
        py::class_<HCpmpm_Reeds_Shepp_State_Space>(module, "HCpmpmReedsSheppStateSpace")
            .def(py::init<double, double, double>(),
                 py::arg("kappa"),
                 py::arg("sigma"),
                 py::arg("discretization") = 0.1);
    steering::bind_common_planner_api(hcpmpm_reeds_shepp_state_space);
}
