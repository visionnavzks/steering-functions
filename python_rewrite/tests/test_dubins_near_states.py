"""
Tests for DubinsStateSpace behaviour when start and goal are very near.

Key scenarios verified:
  1. Exact same state → distance = 0, empty path.
  2. Near-identical position with same heading → distance ≈ 0.
  3. Near-identical states whose headings straddle the 0 / 2π boundary
     (the wrap-around bug fixed in _dubins_shortest).
  4. Small but non-negligible offsets → algorithm produces a valid path
     whose endpoint matches the goal within floating-point tolerance.
  5. Normal (non-degenerate) cases still work.
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

# ---------------------------------------------------------------------------
# Make the python_rewrite package importable when running from any directory.
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_PY_REWRITE = _REPO_ROOT / "python_rewrite"
if str(_PY_REWRITE) not in sys.path:
    sys.path.insert(0, str(_PY_REWRITE))

from steering_functions.dubins_state_space import DubinsStateSpace  # noqa: E402
from steering_functions.state import State  # noqa: E402

_KAPPA = 1.0  # maximum curvature (turning radius = 1 m)
_DISC = 0.01  # discretisation step (small for path-endpoint accuracy)
_EPS = 1e-5   # tolerance for "approximately zero" distance comparisons


def _make_ss() -> DubinsStateSpace:
    return DubinsStateSpace(kappa=_KAPPA, discretization=_DISC)


def _path_endpoint_error(ss: DubinsStateSpace, s1: State, s2: State) -> float:
    """Return Euclidean error between integrated path endpoint and s2."""
    path = ss.get_path(s1, s2)
    if not path:
        # Empty path → endpoint is s1 itself
        return math.hypot(s1.x - s2.x, s1.y - s2.y)
    end = path[-1]
    return math.hypot(end.x - s2.x, end.y - s2.y)


class TestDubinsNearStates(unittest.TestCase):
    """DubinsStateSpace edge cases for near-identical start / goal states."""

    def setUp(self):
        self.ss = _make_ss()

    # ------------------------------------------------------------------
    # 1. Exact same state
    # ------------------------------------------------------------------
    def test_exact_same_state_distance_zero(self):
        s = State(1.5, -2.3, 0.7)
        self.assertAlmostEqual(self.ss.get_distance(s, s), 0.0, places=10)

    def test_exact_same_state_empty_controls(self):
        s = State(1.5, -2.3, 0.7)
        self.assertEqual(self.ss.get_controls(s, s), [])

    def test_exact_same_state_empty_path(self):
        s = State(1.5, -2.3, 0.7)
        self.assertEqual(self.ss.get_path(s, s), [])

    # ------------------------------------------------------------------
    # 2. Near-identical position, same heading
    # ------------------------------------------------------------------
    def test_near_position_same_heading(self):
        """Tiny position offset, same heading → distance ≈ 0."""
        s1 = State(0.0, 0.0, 1.0)
        for dx, dy in [(1e-7, 0), (0, 1e-7), (1e-8, 1e-8)]:
            s2 = State(dx, dy, 1.0)
            dist = self.ss.get_distance(s1, s2)
            self.assertLess(
                dist, _EPS,
                msg=f"Expected near-zero distance for tiny offset ({dx}, {dy}), got {dist}",
            )

    # ------------------------------------------------------------------
    # 3. Heading wrap-around bug (0 / 2π boundary)
    # ------------------------------------------------------------------
    def test_heading_straddle_zero_near_identical(self):
        """
        state1.theta slightly negative, state2.theta slightly positive.
        The actual heading difference is tiny (< _DUBINS_EPS), but before the
        fix abs(alpha - beta) computed ≈ 2π instead of ≈ dtheta, causing the
        near-identical special case to be skipped and returning distance ≈ 2π.
        """
        half_eps = 4e-7  # dtheta = 8e-7 < _DUBINS_EPS = 1e-6
        s1 = State(0.0, 0.0, -half_eps)
        s2 = State(0.0, 0.0,  half_eps)
        dist = self.ss.get_distance(s1, s2)
        self.assertLess(
            dist, _EPS,
            msg=f"Heading straddle 0: expected near-zero distance, got {dist:.6f}",
        )

    def test_heading_straddle_zero_tiny_position_offset(self):
        """Same as above but with a tiny position offset (d > 0 but < eps)."""
        half_eps = 4e-7
        s1 = State(0.0, 0.0, -half_eps)
        s2 = State(1e-7, 0.0,  half_eps)  # d = 1e-7 < _DUBINS_EPS
        dist = self.ss.get_distance(s1, s2)
        self.assertLess(
            dist, _EPS,
            msg=f"Heading straddle 0 with tiny offset: expected near-zero distance, got {dist:.6f}",
        )

    def test_heading_straddle_consistency(self):
        """
        Straddling-0 and non-straddling cases with the same heading difference
        should give the same distance (within floating-point tolerance).
        """
        dtheta = 5e-7  # well below _DUBINS_EPS = 1e-6
        # Straddling 0
        d_straddle = self.ss.get_distance(
            State(0.0, 0.0, -dtheta / 2),
            State(0.0, 0.0,  dtheta / 2),
        )
        # Not straddling (around θ = 1 rad)
        d_normal = self.ss.get_distance(
            State(0.0, 0.0, 1.0 - dtheta / 2),
            State(0.0, 0.0, 1.0 + dtheta / 2),
        )
        self.assertAlmostEqual(
            d_straddle, d_normal, places=8,
            msg=(
                f"Straddling-0 distance ({d_straddle:.6e}) differs from "
                f"non-straddling distance ({d_normal:.6e})"
            ),
        )

    # ------------------------------------------------------------------
    # 4. Slightly larger offset → path reaches goal
    # ------------------------------------------------------------------
    def test_path_reaches_goal_small_offset(self):
        """For a small-but-non-trivial offset the path must end at the goal."""
        s1 = State(0.0, 0.0, 0.0)
        s2 = State(0.05, 0.05, 0.1)
        err = _path_endpoint_error(self.ss, s1, s2)
        self.assertLess(err, 0.05, msg=f"Path endpoint error {err:.4f} too large")

    # ------------------------------------------------------------------
    # 5. Normal (non-degenerate) cases are unaffected
    # ------------------------------------------------------------------
    def test_normal_case_path_reaches_goal(self):
        """Standard case: path must end at the goal within discretisation error."""
        s1 = State(0.0, 0.0, 0.0)
        s2 = State(3.0, 2.0, math.pi / 4)
        err = _path_endpoint_error(self.ss, s1, s2)
        self.assertLess(err, 0.05, msg=f"Path endpoint error {err:.4f} too large")

    def test_normal_case_distance_positive(self):
        s1 = State(0.0, 0.0, 0.0)
        s2 = State(3.0, 2.0, math.pi / 4)
        dist = self.ss.get_distance(s1, s2)
        self.assertGreater(dist, 0.0)

    # ------------------------------------------------------------------
    # 6. No NaN / Inf produced for any near-state combo
    # ------------------------------------------------------------------
    def test_no_nan_or_inf_for_tiny_offsets(self):
        """Distance must be finite for tiny positional and heading offsets."""
        import random
        rng = random.Random(42)
        for _ in range(50):
            dx    = rng.uniform(-1e-6, 1e-6)
            dy    = rng.uniform(-1e-6, 1e-6)
            dtheta = rng.uniform(-1e-4, 1e-4)
            s1 = State(0.0, 0.0, 0.0)
            s2 = State(dx, dy, dtheta)
            dist = self.ss.get_distance(s1, s2)
            self.assertTrue(
                math.isfinite(dist),
                msg=f"Non-finite distance {dist} for dx={dx:.2e}, dy={dy:.2e}, dtheta={dtheta:.2e}",
            )


if __name__ == "__main__":
    unittest.main()
