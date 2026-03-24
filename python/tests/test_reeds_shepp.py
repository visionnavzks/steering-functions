"""Tests for the Reeds-Shepp state space."""

import math
import random
import sys
import unittest

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent.parent))

from steering_functions import ReedsSheppStateSpace, State

KAPPA = 1.0
DISCRETIZATION = 0.1
OPERATING_REGION_X = 20.0
OPERATING_REGION_Y = 20.0
OPERATING_REGION_THETA = 2 * math.pi
EPS_DISTANCE = 0.05
EPS_YAW = 0.05
SAMPLES = 1000


def random_state(rng: random.Random) -> State:
    return State(
        x=rng.uniform(-OPERATING_REGION_X / 2, OPERATING_REGION_X / 2),
        y=rng.uniform(-OPERATING_REGION_Y / 2, OPERATING_REGION_Y / 2),
        theta=rng.uniform(-OPERATING_REGION_THETA / 2, OPERATING_REGION_THETA / 2),
    )


class TestReedsSheppStateSpace(unittest.TestCase):

    def setUp(self) -> None:
        self.ss = ReedsSheppStateSpace(KAPPA, DISCRETIZATION)
        self.rng = random.Random(42)

    def test_trivial_same_state(self) -> None:
        """Zero-length path for identical start and goal."""
        s = State(1.0, 2.0, math.pi / 4)
        controls = self.ss.get_controls(s, s)
        path = self.ss.integrate(s, controls)
        self.assertLessEqual(len(path), 2)

    def test_straight_path(self) -> None:
        """Simple forward straight path."""
        s1 = State(0.0, 0.0, 0.0)
        s2 = State(5.0, 0.0, 0.0)
        path = self.ss.get_path(s1, s2)
        self.assertGreater(len(path), 1)
        end = path[-1]
        self.assertAlmostEqual(end.x, s2.x, delta=EPS_DISTANCE)
        self.assertAlmostEqual(end.y, s2.y, delta=EPS_DISTANCE)

    def test_reverse_path(self) -> None:
        """RS can plan a path that involves reversing."""
        # goal is directly behind start → must reverse
        s1 = State(0.0, 0.0, 0.0)
        s2 = State(-3.0, 0.0, 0.0)
        d = self.ss.get_distance(s1, s2)
        self.assertGreater(d, 0.0)
        path = self.ss.get_path(s1, s2)
        self.assertGreater(len(path), 1)
        end = path[-1]
        self.assertAlmostEqual(end.x, s2.x, delta=EPS_DISTANCE)
        self.assertAlmostEqual(end.y, s2.y, delta=EPS_DISTANCE)

    def test_endpoint_accuracy(self) -> None:
        """Path endpoint should be close to the goal state."""
        for _ in range(SAMPLES):
            s1 = random_state(self.rng)
            s2 = random_state(self.rng)
            path = self.ss.get_path(s1, s2)
            if len(path) < 2:
                continue
            end = path[-1]
            dist = math.hypot(end.x - s2.x, end.y - s2.y)
            self.assertLess(
                dist,
                EPS_DISTANCE,
                msg=f"Endpoint distance {dist:.4f} too large for s1={s1}, s2={s2}",
            )
            angle_err = abs(
                (end.theta - s2.theta + math.pi) % (2 * math.pi) - math.pi
            )
            self.assertLess(
                angle_err,
                EPS_YAW,
                msg=f"Angle error {angle_err:.4f} too large for s1={s1}, s2={s2}",
            )

    def test_distance_non_negative(self) -> None:
        """All distances must be non-negative."""
        for _ in range(200):
            s1 = random_state(self.rng)
            s2 = random_state(self.rng)
            self.assertGreaterEqual(self.ss.get_distance(s1, s2), 0.0)

    def test_distance_nearly_symmetric(self) -> None:
        """RS distance is symmetric (both directions yield same length)."""
        for _ in range(200):
            s1 = random_state(self.rng)
            s2 = random_state(self.rng)
            d12 = self.ss.get_distance(s1, s2)
            d21 = self.ss.get_distance(s2, s1)
            self.assertAlmostEqual(
                d12, d21, delta=1e-6,
                msg=f"RS distance not symmetric: {d12} vs {d21}"
            )

    def test_all_controls_not_empty(self) -> None:
        """get_all_controls returns at least one path for non-trivial states."""
        s1 = State(0, 0, 0)
        s2 = State(2, 2, math.pi / 4)
        all_ctrl = self.ss.get_all_controls(s1, s2)
        self.assertGreater(len(all_ctrl), 0)

    def test_all_paths_endpoint(self) -> None:
        """All paths from get_all_paths should end at the goal."""
        rng = random.Random(7)
        for _ in range(50):
            s1 = random_state(rng)
            s2 = random_state(rng)
            all_paths = self.ss.get_all_paths(s1, s2)
            self.assertGreater(len(all_paths), 0)
            for path in all_paths:
                if len(path) < 2:
                    continue
                end = path[-1]
                dist = math.hypot(end.x - s2.x, end.y - s2.y)
                self.assertLess(dist, EPS_DISTANCE)

    def test_interpolate_endpoints(self) -> None:
        """Interpolate at t=0 returns start, at t=1 returns end."""
        s1 = State(0, 0, 0)
        s2 = State(3, 4, math.pi / 3)
        controls = self.ss.get_controls(s1, s2)
        if not controls:
            return
        state_at_0 = self.ss.interpolate(s1, controls, 0.0)
        state_at_1 = self.ss.interpolate(s1, controls, 1.0)
        self.assertAlmostEqual(state_at_0.x, s1.x, delta=1e-9)
        self.assertAlmostEqual(state_at_0.y, s1.y, delta=1e-9)
        self.assertAlmostEqual(state_at_1.x, s2.x, delta=EPS_DISTANCE)
        self.assertAlmostEqual(state_at_1.y, s2.y, delta=EPS_DISTANCE)

    def test_trivial_all_controls(self) -> None:
        """Same-state query returns the trivial single-control path."""
        s = State(1, 2, 0.5)
        all_ctrl = self.ss.get_all_controls(s, s)
        self.assertEqual(len(all_ctrl), 1)
        self.assertEqual(len(all_ctrl[0]), 1)
        self.assertEqual(all_ctrl[0][0].delta_s, 0.0)


if __name__ == "__main__":
    unittest.main()
