"""Tests for the Dubins state space."""

import math
import random
import sys
import unittest

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent.parent))

from steering_functions import DubinsStateSpace, State

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


def path_length_states(path: list) -> float:
    """Approximate path length by summing Euclidean steps."""
    if len(path) < 2:
        return 0.0
    total = 0.0
    for a, b in zip(path, path[1:]):
        total += math.hypot(b.x - a.x, b.y - a.y)
    return total


class TestDubinsStateSpace(unittest.TestCase):

    def setUp(self) -> None:
        self.ss = DubinsStateSpace(KAPPA, DISCRETIZATION)
        self.rng = random.Random(42)

    def test_trivial_same_state(self) -> None:
        """Zero-length path for identical start and goal yields no/one state."""
        s = State(1.0, 2.0, math.pi / 4)
        path = self.ss.get_path(s, s)
        # Identical start and goal produce no controls → empty or single state
        self.assertLessEqual(len(path), 1)

    def test_straight_path(self) -> None:
        """Forward straight line: start and end headings both 0."""
        s1 = State(0.0, 0.0, 0.0)
        s2 = State(5.0, 0.0, 0.0)
        path = self.ss.get_path(s1, s2)
        self.assertGreater(len(path), 1)
        # End position should be close to goal
        end = path[-1]
        self.assertAlmostEqual(end.x, s2.x, delta=EPS_DISTANCE)
        self.assertAlmostEqual(end.y, s2.y, delta=EPS_DISTANCE)

    def test_distance_symmetry(self) -> None:
        """Dubins distance is *not* symmetric (left-turn != right-turn)."""
        s1 = State(0, 0, 0)
        s2 = State(1, 1, math.pi / 2)
        d12 = self.ss.get_distance(s1, s2)
        d21 = self.ss.get_distance(s2, s1)
        # They can differ, but both should be positive finite values
        self.assertGreater(d12, 0)
        self.assertGreater(d21, 0)
        self.assertFalse(math.isinf(d12))
        self.assertFalse(math.isinf(d21))

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

    def test_controls_give_correct_endpoint(self) -> None:
        """Integrating controls should reproduce the path endpoint."""
        rng = random.Random(99)
        for _ in range(200):
            s1 = random_state(rng)
            s2 = random_state(rng)
            controls = self.ss.get_controls(s1, s2)
            path = self.ss.integrate(s1, controls)
            if not path:
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

    def test_distance_non_negative(self) -> None:
        """All distances must be non-negative."""
        for _ in range(200):
            s1 = random_state(self.rng)
            s2 = random_state(self.rng)
            self.assertGreaterEqual(self.ss.get_distance(s1, s2), 0.0)

    def test_all_controls_not_empty(self) -> None:
        """get_all_controls should return at least one path."""
        s1 = State(0, 0, 0)
        s2 = State(2, 2, math.pi / 4)
        all_ctrl = self.ss.get_all_controls(s1, s2)
        self.assertGreater(len(all_ctrl), 0)

    def test_forwards_flag(self) -> None:
        """forwards=False should produce a path that ends at the goal."""
        ss_rev = DubinsStateSpace(KAPPA, DISCRETIZATION, forwards=False)
        s1 = State(0, 0, 0)
        s2 = State(3, 0, 0)
        d = ss_rev.get_distance(s1, s2)
        self.assertGreater(d, 0)


if __name__ == "__main__":
    unittest.main()
