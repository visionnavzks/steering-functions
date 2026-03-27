"""Tests for SteeringPath high-level API."""

import math
import unittest

from steering_functions import PathType, State, SteeringPath


KAPPA_MAX = 1.0
SIGMA_MAX = 1.0
DISCRETIZATION = 0.1
POSITION_TOLERANCE = 0.15
THETA_TOLERANCE = 0.15


def angle_error(a: float, b: float) -> float:
    return abs(math.atan2(math.sin(a - b), math.cos(a - b)))


ALL_PATH_TYPES = [
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


class TestSteeringPathProperties(unittest.TestCase):
    """Test SteeringPath read-only properties."""

    def test_path_type_property(self) -> None:
        for pt in ALL_PATH_TYPES:
            planner = SteeringPath(pt, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)
            self.assertEqual(planner.path_type, pt)

    def test_kappa_max_property(self) -> None:
        planner = SteeringPath(PathType.RS, 2.0, 0.5, 0.1)
        self.assertAlmostEqual(planner.kappa_max, 2.0)

    def test_sigma_max_property(self) -> None:
        planner = SteeringPath(PathType.HC_RS, 1.0, 0.75, 0.1)
        self.assertAlmostEqual(planner.sigma_max, 0.75)

    def test_discretization_property(self) -> None:
        planner = SteeringPath(PathType.RS, 1.0, 1.0, 0.05)
        self.assertAlmostEqual(planner.discretization, 0.05)


class TestSteeringPathValidation(unittest.TestCase):
    """Test SteeringPath constructor parameter validation."""

    def test_zero_kappa_raises(self) -> None:
        with self.assertRaises(Exception):
            SteeringPath(PathType.RS, 0.0, 1.0, 0.1)

    def test_negative_kappa_raises(self) -> None:
        with self.assertRaises(Exception):
            SteeringPath(PathType.RS, -1.0, 1.0, 0.1)

    def test_negative_sigma_raises(self) -> None:
        with self.assertRaises(Exception):
            SteeringPath(PathType.HC_RS, 1.0, -1.0, 0.1)

    def test_zero_discretization_raises(self) -> None:
        with self.assertRaises(Exception):
            SteeringPath(PathType.RS, 1.0, 1.0, 0.0)

    def test_negative_discretization_raises(self) -> None:
        with self.assertRaises(Exception):
            SteeringPath(PathType.RS, 1.0, 1.0, -0.1)

    def test_valid_params_no_raise(self) -> None:
        planner = SteeringPath(PathType.RS, 1.0, 0.0, 0.1)
        self.assertIsNotNone(planner)

    def test_zero_sigma_valid_for_non_cc(self) -> None:
        """sigma_max=0 is valid for Dubins and RS (no curvature rate needed)."""
        for pt in [PathType.DUBINS, PathType.RS]:
            planner = SteeringPath(pt, 1.0, 0.0, 0.1)
            self.assertIsNotNone(planner)


class TestSteeringPathNoneType(unittest.TestCase):
    """Test behavior with PathType.NONE."""

    def test_none_type_creates_planner(self) -> None:
        """PathType.NONE falls through to default RS planner."""
        planner = SteeringPath(PathType.NONE, 1.0, 1.0, 0.1)
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(3.0, 2.0, 0.5, 0.0)
        path = planner.computeShortestPath(start, goal)
        self.assertGreater(len(path), 0)


class TestSteeringPathGoalReaching(unittest.TestCase):
    """Test that all path types reach their goals for various scenarios."""

    def _check_goal_reached(self, path_type: PathType, start: State, goal: State) -> None:
        planner = SteeringPath(path_type, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)
        path = planner.computeShortestPath(start, goal)

        self.assertGreater(len(path), 0, f"{path_type} returned empty path")
        end = path[-1]
        dist = math.hypot(end.x - goal.x, end.y - goal.y)
        self.assertLess(dist, POSITION_TOLERANCE,
                        f"{path_type}: position error {dist:.4f}")
        self.assertLess(angle_error(end.theta, goal.theta), THETA_TOLERANCE,
                        f"{path_type}: angle error too large")

    def test_straight_ahead(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(5.0, 0.0, 0.0, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                self._check_goal_reached(pt, start, goal)

    def test_right_turn(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(3.0, -3.0, -math.pi / 2, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                self._check_goal_reached(pt, start, goal)

    def test_left_turn(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(3.0, 3.0, math.pi / 2, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                self._check_goal_reached(pt, start, goal)

    def test_u_turn_reeds_shepp(self) -> None:
        """U-turn requires backward motion; only RS-like planners support this."""
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(-2.0, 0.0, math.pi, 0.0)
        rs_types = [PathType.RS, PathType.CC00_RS, PathType.HC_RS,
                    PathType.HC00_RS, PathType.HC0PM_RS,
                    PathType.HCPM0_RS, PathType.HCPMPM_RS]
        for pt in rs_types:
            with self.subTest(path_type=pt):
                self._check_goal_reached(pt, start, goal)

    def test_with_initial_curvature(self) -> None:
        """Start with non-zero curvature."""
        start = State(0.0, 0.0, 0.0, 0.5)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                self._check_goal_reached(pt, start, goal)

    def test_diagonal_goal(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 4.0, math.pi / 4, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                self._check_goal_reached(pt, start, goal)


class TestSteeringPathSameStartGoal(unittest.TestCase):
    """Test paths when start equals goal.

    Note: Most planners return empty paths when start == goal (distance is 0).
    This is expected behavior from the C++ library, not a bug.
    """

    def test_same_start_goal_distance_is_zero(self) -> None:
        state = State(1.5, -2.0, 0.75, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                planner = SteeringPath(pt, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)
                controls = planner.computeShortestControlSequence(state, state)
                total_length = sum(abs(c.delta_s) for c in controls)
                # Distance should be very small (zero or near-zero)
                self.assertLess(total_length, 0.1)

    def test_same_start_goal_path_is_valid_if_nonempty(self) -> None:
        """If a path is returned, its end state should be close to start."""
        state = State(1.5, -2.0, 0.75, 0.2)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                planner = SteeringPath(pt, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)
                path = planner.computeShortestPath(state, state)
                if len(path) > 0:
                    dist = math.hypot(path[-1].x - state.x, path[-1].y - state.y)
                    self.assertLess(dist, POSITION_TOLERANCE)


class TestSteeringPathAllPaths(unittest.TestCase):
    """Test computeAllPaths and computeAllControlSequences."""

    def test_all_paths_returns_multiple(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                planner = SteeringPath(pt, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)
                all_paths = planner.computeAllPaths(start, goal)
                self.assertGreater(len(all_paths), 0)
                for path in all_paths:
                    self.assertGreater(len(path), 0)

    def test_all_controls_returns_multiple(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                planner = SteeringPath(pt, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)
                all_controls = planner.computeAllControlSequences(start, goal)
                self.assertGreater(len(all_controls), 0)
                for controls in all_controls:
                    self.assertGreater(len(controls), 0)

    def test_shortest_path_is_not_longer_than_any_alternative(self) -> None:
        """The shortest path should be <= all alternative paths in length."""
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                planner = SteeringPath(pt, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)
                shortest_controls = planner.computeShortestControlSequence(start, goal)
                all_controls = planner.computeAllControlSequences(start, goal)
                shortest_len = sum(abs(c.delta_s) for c in shortest_controls)
                for controls in all_controls:
                    alt_len = sum(abs(c.delta_s) for c in controls)
                    self.assertLessEqual(shortest_len, alt_len + 0.01)


class TestSteeringPathControlsPathConsistency(unittest.TestCase):
    """Test that control sequences and paths are consistent."""

    def test_controls_nonempty_when_path_nonempty(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                planner = SteeringPath(pt, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)
                controls = planner.computeShortestControlSequence(start, goal)
                path = planner.computeShortestPath(start, goal)
                if len(path) > 1:
                    self.assertGreater(len(controls), 0)

    def test_path_length_positive(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for pt in ALL_PATH_TYPES:
            with self.subTest(path_type=pt):
                planner = SteeringPath(pt, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)
                controls = planner.computeShortestControlSequence(start, goal)
                total = sum(abs(c.delta_s) for c in controls)
                self.assertGreater(total, 0.0)


class TestDifferentDiscretizations(unittest.TestCase):
    """Test that different discretization values produce valid paths."""

    def test_fine_discretization(self) -> None:
        planner = SteeringPath(PathType.RS, KAPPA_MAX, SIGMA_MAX, 0.01)
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(3.0, 1.0, 0.5, 0.0)
        path = planner.computeShortestPath(start, goal)
        self.assertGreater(len(path), 0)
        end = path[-1]
        self.assertLess(math.hypot(end.x - goal.x, end.y - goal.y), POSITION_TOLERANCE)

    def test_coarse_discretization(self) -> None:
        planner = SteeringPath(PathType.RS, KAPPA_MAX, SIGMA_MAX, 0.5)
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(3.0, 1.0, 0.5, 0.0)
        path = planner.computeShortestPath(start, goal)
        self.assertGreater(len(path), 0)

    def test_finer_discretization_yields_more_states(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(5.0, 3.0, 1.0, 0.0)
        fine = SteeringPath(PathType.RS, KAPPA_MAX, SIGMA_MAX, 0.05)
        coarse = SteeringPath(PathType.RS, KAPPA_MAX, SIGMA_MAX, 0.5)
        fine_path = fine.computeShortestPath(start, goal)
        coarse_path = coarse.computeShortestPath(start, goal)
        self.assertGreater(len(fine_path), len(coarse_path))


if __name__ == "__main__":
    unittest.main()
