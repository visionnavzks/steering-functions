from __future__ import annotations

import math
import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from steering_functions import (  # noqa: E402
    DubinsStateSpace,
    PathType,
    ReedsSheppStateSpace,
    State,
    SteeringPath,
    pify,
)


def path_length(path: list[State]) -> float:
    total = 0.0
    previous = path[0]
    for current in path[1:]:
        total += math.hypot(current.x - previous.x, current.y - previous.y)
        previous = current
    return total


def control_length(controls) -> float:
    return sum(abs(control.delta_s) for control in controls)


class SteeringFunctionsPythonPortTest(unittest.TestCase):
    def assertStateNear(self, state: State, target: State, *, distance: float = 0.05, theta: float = 0.05) -> None:
        self.assertLessEqual(math.hypot(state.x - target.x, state.y - target.y), distance)
        self.assertLessEqual(abs(pify(state.theta - target.theta)), theta)

    def test_state_equality_uses_tolerance(self) -> None:
        self.assertEqual(State(x=1.0, y=2.0), State(x=1.0 + 1e-7, y=2.0 - 1e-7))

    def test_dubins_path_reaches_goal_and_matches_control_length(self) -> None:
        planner = DubinsStateSpace(1.0, 0.05, True)
        start = State(0.0, 0.0, 0.0)
        goal = State(4.0, 4.0, math.pi / 2.0)
        controls = planner.get_controls(start, goal)
        path = planner.get_path(start, goal)
        self.assertStateNear(path[-1], goal)
        self.assertAlmostEqual(path_length(path), control_length(controls), delta=0.1)
        self.assertTrue(all(control.sigma == 0.0 for control in controls))

    def test_dubins_get_all_controls_returns_forward_and_reverse_candidates(self) -> None:
        planner = DubinsStateSpace(1.0, 0.05, True)
        start = State(0.0, 0.0, 0.0)
        goal = State(2.0, 1.0, 0.25)
        all_controls = planner.get_all_controls(start, goal)
        self.assertEqual(len(all_controls), 2)
        self.assertGreaterEqual(control_length(all_controls[0]), 0.0)
        self.assertGreaterEqual(control_length(all_controls[1]), 0.0)

    def test_reeds_shepp_supports_reverse_segments(self) -> None:
        planner = ReedsSheppStateSpace(1.0, 0.05)
        start = State(0.0, 0.0, 0.0)
        goal = State(-2.0, 0.0, 0.0)
        controls = planner.get_controls(start, goal)
        path = planner.get_path(start, goal)
        self.assertStateNear(path[-1], goal)
        self.assertAlmostEqual(path_length(path), control_length(controls), delta=0.1)
        self.assertTrue(any(control.delta_s < 0.0 for control in controls))

    def test_steering_path_facade_supports_dubins_and_rs(self) -> None:
        start = State(0.0, 0.0, 0.0)
        goal = State(3.0, 2.0, 0.5)
        dubins = SteeringPath(PathType.DUBINS, 1.0, 0.0, 0.05)
        rs = SteeringPath(PathType.RS, 1.0, 0.0, 0.05)
        self.assertStateNear(dubins.compute_shortest_path(start, goal)[-1], goal)
        self.assertStateNear(rs.compute_shortest_path(start, goal)[-1], goal)

    def test_unsupported_python_planner_types_raise(self) -> None:
        with self.assertRaises(NotImplementedError):
            SteeringPath(PathType.CC_DUBINS, 1.0, 1.0, 0.05)


if __name__ == "__main__":
    unittest.main()
