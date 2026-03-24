import math
import unittest

from steering_functions import Control, PathType, State, SteeringPath
from steering_functions import DubinsStateSpace, ReedsSheppStateSpace


KAPPA_MAX = 1.0
SIGMA_MAX = 1.0
DISCRETIZATION = 0.1
POSITION_TOLERANCE = 0.15
THETA_TOLERANCE = 0.15


def angle_error(angle_a: float, angle_b: float) -> float:
    return abs(math.atan2(math.sin(angle_a - angle_b), math.cos(angle_a - angle_b)))


class PythonBindingsTest(unittest.TestCase):
    def test_state_equality_and_repr(self) -> None:
        state_a = State(1.0, 2.0, 0.5, 0.25)
        state_b = State(1.0, 2.0, 0.5, 0.25)

        self.assertEqual(state_a, state_b)
        self.assertIn("State(", str(state_a))

    def test_control_repr(self) -> None:
        control = Control(1.5, 0.5, 0.0)
        self.assertIn("Control Segment", str(control))

    def test_all_supported_path_types_reach_goal(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)

        supported_path_types = [
            PathType.CC_DUBINS,
            PathType.CC00_DUBINS,
            PathType.CC0PM_DUBINS,
            PathType.CCPM0_DUBINS,
            PathType.CCPMPM_DUBINS,
            PathType.DUBINS,
            PathType.CC00_RS,
            PathType.HC_RS,
            PathType.HC00_RS,
            PathType.HC0PM_RS,
            PathType.HCPM0_RS,
            PathType.HCPMPM_RS,
            PathType.RS,
        ]

        for path_type in supported_path_types:
            with self.subTest(path_type=path_type):
                planner = SteeringPath(path_type, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)

                shortest_controls = planner.computeShortestControlSequence(start, goal)
                shortest_path = planner.computeShortestPath(start, goal)
                all_controls = planner.computeAllControlSequences(start, goal)
                all_paths = planner.computeAllPaths(start, goal)

                self.assertGreater(len(shortest_controls), 0)
                self.assertGreater(len(shortest_path), 1)
                self.assertGreater(len(all_controls), 0)
                self.assertGreater(len(all_paths), 0)

                end_state = shortest_path[-1]
                self.assertLess(math.hypot(end_state.x - goal.x, end_state.y - goal.y), POSITION_TOLERANCE)
                self.assertLess(angle_error(end_state.theta, goal.theta), THETA_TOLERANCE)

    def test_low_level_state_spaces_are_available(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)

        planners = [
            DubinsStateSpace(KAPPA_MAX, DISCRETIZATION, True),
            ReedsSheppStateSpace(KAPPA_MAX, DISCRETIZATION),
        ]

        for planner in planners:
            with self.subTest(planner=type(planner).__name__):
                path = planner.get_path(start, goal)
                self.assertGreater(len(path), 1)
                self.assertLess(math.hypot(path[-1].x - goal.x, path[-1].y - goal.y), POSITION_TOLERANCE)


if __name__ == "__main__":
    unittest.main()
