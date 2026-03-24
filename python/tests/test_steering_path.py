import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from steering_functions import PathType, State, SteeringPath


class SteeringPathTest(unittest.TestCase):
    def test_constructor_validation(self):
        with self.assertRaises(ValueError):
            SteeringPath(PathType.DUBINS, 0.0, 0.1, 0.1)
        with self.assertRaises(ValueError):
            SteeringPath(PathType.DUBINS, 1.0, -0.1, 0.1)
        with self.assertRaises(ValueError):
            SteeringPath(PathType.DUBINS, 1.0, 0.1, 0.0)

    def test_compute_shortest_control_sequence(self):
        planner = SteeringPath(PathType.DUBINS, 1.0, 1.0, 0.1)
        start = State(0.0, 0.0, 0.0)
        goal = State(3.0, 4.0, 0.0)
        controls = planner.compute_shortest_control_sequence(start, goal)
        self.assertEqual(1, len(controls))
        self.assertAlmostEqual(5.0, controls[0].delta_s)
        self.assertAlmostEqual(0.0, controls[0].kappa)
        self.assertAlmostEqual(0.0, controls[0].sigma)

    def test_constructor_rejects_none_path_type(self):
        with self.assertRaises(ValueError):
            SteeringPath(PathType.NONE, 1.0, 0.1, 0.1)

    def test_compute_shortest_path_contains_start(self):
        planner = SteeringPath(PathType.RS, 1.0, 1.0, 0.5)
        start = State(0.0, 0.0, 0.0)
        goal = State(1.0, 0.0, 0.0)
        path = planner.compute_shortest_path(start, goal)
        self.assertGreaterEqual(len(path), 2)
        self.assertAlmostEqual(start.x, path[0].x)
        self.assertAlmostEqual(start.y, path[0].y)


if __name__ == "__main__":
    unittest.main()
