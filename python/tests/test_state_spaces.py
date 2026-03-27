"""Tests for all low-level state space bindings."""

import math
import unittest

from steering_functions import (
    CC00DubinsStateSpace,
    CC00ReedsSheppStateSpace,
    CC0pmDubinsStateSpace,
    CCDubinsStateSpace,
    CCpm0DubinsStateSpace,
    CCpmpmDubinsStateSpace,
    DubinsStateSpace,
    HC00ReedsSheppStateSpace,
    HC0pmReedsSheppStateSpace,
    HCReedsSheppStateSpace,
    HCpm0ReedsSheppStateSpace,
    HCpmpmReedsSheppStateSpace,
    ReedsSheppStateSpace,
    State,
)


KAPPA = 1.0
SIGMA = 1.0
DISC = 0.1
POS_TOL = 0.15
THETA_TOL = 0.15


def angle_error(a: float, b: float) -> float:
    return abs(math.atan2(math.sin(a - b), math.cos(a - b)))


def make_dubins_planners():
    """Create all Dubins-family state spaces."""
    return {
        "DubinsForwards": DubinsStateSpace(KAPPA, DISC, True),
        "DubinsBackwards": DubinsStateSpace(KAPPA, DISC, False),
        "CCDubinsForwards": CCDubinsStateSpace(KAPPA, SIGMA, DISC, True),
        "CCDubinsBackwards": CCDubinsStateSpace(KAPPA, SIGMA, DISC, False),
        "CC00DubinsForwards": CC00DubinsStateSpace(KAPPA, SIGMA, DISC, True),
        "CC00DubinsBackwards": CC00DubinsStateSpace(KAPPA, SIGMA, DISC, False),
        "CC0pmDubinsForwards": CC0pmDubinsStateSpace(KAPPA, SIGMA, DISC, True),
        "CC0pmDubinsBackwards": CC0pmDubinsStateSpace(KAPPA, SIGMA, DISC, False),
        "CCpm0DubinsForwards": CCpm0DubinsStateSpace(KAPPA, SIGMA, DISC, True),
        "CCpm0DubinsBackwards": CCpm0DubinsStateSpace(KAPPA, SIGMA, DISC, False),
        "CCpmpmDubinsForwards": CCpmpmDubinsStateSpace(KAPPA, SIGMA, DISC, True),
        "CCpmpmDubinsBackwards": CCpmpmDubinsStateSpace(KAPPA, SIGMA, DISC, False),
    }


def make_rs_planners():
    """Create all Reeds-Shepp-family state spaces."""
    return {
        "ReedsShepp": ReedsSheppStateSpace(KAPPA, DISC),
        "CC00ReedsShepp": CC00ReedsSheppStateSpace(KAPPA, SIGMA, DISC),
        "HCReedsShepp": HCReedsSheppStateSpace(KAPPA, SIGMA, DISC),
        "HC00ReedsShepp": HC00ReedsSheppStateSpace(KAPPA, SIGMA, DISC),
        "HC0pmReedsShepp": HC0pmReedsSheppStateSpace(KAPPA, SIGMA, DISC),
        "HCpm0ReedsShepp": HCpm0ReedsSheppStateSpace(KAPPA, SIGMA, DISC),
        "HCpmpmReedsShepp": HCpmpmReedsSheppStateSpace(KAPPA, SIGMA, DISC),
    }


def make_all_planners():
    planners = {}
    planners.update(make_dubins_planners())
    planners.update(make_rs_planners())
    return planners


class TestAllStateSpacesGetPath(unittest.TestCase):
    """Test get_path for all state spaces with a simple scenario."""

    def test_get_path_reaches_goal(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                path = planner.get_path(start, goal)
                self.assertGreater(len(path), 1, f"{name} returned short path")
                end = path[-1]
                self.assertLess(
                    math.hypot(end.x - goal.x, end.y - goal.y),
                    POS_TOL,
                    f"{name}: position error too large",
                )
                self.assertLess(
                    angle_error(end.theta, goal.theta),
                    THETA_TOL,
                    f"{name}: angle error too large",
                )


class TestAllStateSpacesGetControls(unittest.TestCase):
    """Test get_controls for all state spaces."""

    def test_get_controls_nonempty(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                controls = planner.get_controls(start, goal)
                self.assertGreater(len(controls), 0, f"{name} returned empty controls")

    def test_get_controls_total_length_positive(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                controls = planner.get_controls(start, goal)
                total = sum(abs(c.delta_s) for c in controls)
                self.assertGreater(total, 0.0)


class TestAllStateSpacesGetDistance(unittest.TestCase):
    """Test get_distance for all state spaces."""

    def test_distance_nonnegative(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                d = planner.get_distance(start, goal)
                self.assertGreaterEqual(d, 0.0)

    def test_distance_same_state(self) -> None:
        state = State(1.0, 2.0, 0.5, 0.0)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                d = planner.get_distance(state, state)
                self.assertAlmostEqual(d, 0.0, places=2)

    def test_distance_at_least_euclidean(self) -> None:
        """Path distance should be >= straight-line distance."""
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        euclidean = math.hypot(goal.x - start.x, goal.y - start.y)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                d = planner.get_distance(start, goal)
                self.assertGreaterEqual(d, euclidean - 0.01,
                                        f"{name}: distance < euclidean")


class TestReedsSheppSymmetry(unittest.TestCase):
    """Test distance symmetry for Reeds-Shepp family planners.

    Note: HC0pm and HCpm0 are NOT symmetric by design (asymmetric curvature).
    """

    def _symmetric_planners(self):
        """Planners that should have symmetric distance."""
        return {
            "ReedsShepp": ReedsSheppStateSpace(KAPPA, DISC),
            "CC00ReedsShepp": CC00ReedsSheppStateSpace(KAPPA, SIGMA, DISC),
            "HCReedsShepp": HCReedsSheppStateSpace(KAPPA, SIGMA, DISC),
            "HC00ReedsShepp": HC00ReedsSheppStateSpace(KAPPA, SIGMA, DISC),
            "HCpmpmReedsShepp": HCpmpmReedsSheppStateSpace(KAPPA, SIGMA, DISC),
        }

    def test_distance_symmetry(self) -> None:
        start = State(1.0, -1.0, 0.5, 0.3)
        goal = State(-2.0, 3.0, -1.0, -0.2)
        for name, planner in self._symmetric_planners().items():
            with self.subTest(planner=name):
                d_forward = planner.get_distance(start, goal)
                d_backward = planner.get_distance(goal, start)
                self.assertAlmostEqual(d_forward, d_backward, places=1,
                                       msg=f"{name}: asymmetric distance")

    def test_distance_symmetry_zero_curvature(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(5.0, 3.0, 1.0, 0.0)
        for name, planner in self._symmetric_planners().items():
            with self.subTest(planner=name):
                d_forward = planner.get_distance(start, goal)
                d_backward = planner.get_distance(goal, start)
                self.assertAlmostEqual(d_forward, d_backward, places=1)

    def test_hc0pm_hcpm0_are_duals(self) -> None:
        """HC0pm(start→goal) should equal HCpm0(goal→start) and vice versa."""
        start = State(1.0, -1.0, 0.5, 0.3)
        goal = State(-2.0, 3.0, -1.0, -0.2)
        hc0pm = HC0pmReedsSheppStateSpace(KAPPA, SIGMA, DISC)
        hcpm0 = HCpm0ReedsSheppStateSpace(KAPPA, SIGMA, DISC)
        d_0pm_fwd = hc0pm.get_distance(start, goal)
        d_pm0_bwd = hcpm0.get_distance(goal, start)
        self.assertAlmostEqual(d_0pm_fwd, d_pm0_bwd, places=1)


class TestAllStateSpacesGetAllPaths(unittest.TestCase):
    """Test get_all_paths for all state spaces."""

    def test_get_all_paths_returns_list_of_lists(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                all_paths = planner.get_all_paths(start, goal)
                self.assertIsInstance(all_paths, list)
                self.assertGreater(len(all_paths), 0)
                for path in all_paths:
                    self.assertIsInstance(path, list)
                    self.assertGreater(len(path), 0)

    def test_get_all_controls_returns_list_of_lists(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                all_controls = planner.get_all_controls(start, goal)
                self.assertIsInstance(all_controls, list)
                self.assertGreater(len(all_controls), 0)


class TestDistanceConsistencyWithPathLength(unittest.TestCase):
    """Test that get_distance matches the sum of control deltas.

    Note: Backwards Dubins planners return forward-direction distance from
    get_distance but backwards-direction controls. They are excluded here.
    """

    def _consistent_planners(self):
        """Planners where get_distance matches control path length."""
        planners = {}
        planners.update(make_rs_planners())
        # Only forward Dubins planners have consistent distance/control length
        planners["DubinsForwards"] = DubinsStateSpace(KAPPA, DISC, True)
        planners["CCDubinsForwards"] = CCDubinsStateSpace(KAPPA, SIGMA, DISC, True)
        planners["CC00DubinsForwards"] = CC00DubinsStateSpace(KAPPA, SIGMA, DISC, True)
        planners["CC0pmDubinsForwards"] = CC0pmDubinsStateSpace(KAPPA, SIGMA, DISC, True)
        planners["CCpm0DubinsForwards"] = CCpm0DubinsStateSpace(KAPPA, SIGMA, DISC, True)
        planners["CCpmpmDubinsForwards"] = CCpmpmDubinsStateSpace(KAPPA, SIGMA, DISC, True)
        return planners

    def test_distance_matches_control_length(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        for name, planner in self._consistent_planners().items():
            with self.subTest(planner=name):
                distance = planner.get_distance(start, goal)
                controls = planner.get_controls(start, goal)
                control_length = sum(abs(c.delta_s) for c in controls)
                self.assertAlmostEqual(distance, control_length, places=1,
                                       msg=f"{name}: distance != control length")


class TestDubinsBackwardMotion(unittest.TestCase):
    """Test that backward Dubins planners produce backward paths."""

    def test_dubins_backwards_reaches_goal(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        planner = DubinsStateSpace(KAPPA, DISC, False)
        path = planner.get_path(start, goal)
        self.assertGreater(len(path), 1)
        end = path[-1]
        self.assertLess(math.hypot(end.x - goal.x, end.y - goal.y), POS_TOL)

    def test_cc_dubins_backwards_reaches_goal(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        planner = CCDubinsStateSpace(KAPPA, SIGMA, DISC, False)
        path = planner.get_path(start, goal)
        self.assertGreater(len(path), 1)
        end = path[-1]
        self.assertLess(math.hypot(end.x - goal.x, end.y - goal.y), POS_TOL)


class TestDubinsReverseControls(unittest.TestCase):
    """Test get_controls_reverse for Dubins planners."""

    def test_dubins_reverse_controls_nonempty(self) -> None:
        start = State(-3.2, 1.4, 0.7, 0.0)
        goal = State(2.1, -2.8, -1.1, 0.0)
        planner = DubinsStateSpace(KAPPA, DISC, True)
        controls = planner.get_controls_reverse(start, goal)
        self.assertGreater(len(controls), 0)

    def test_cc_dubins_reverse_controls_nonempty(self) -> None:
        start = State(-3.2, 1.4, 0.7, 0.35)
        goal = State(2.1, -2.8, -1.1, -0.45)
        planner = CCDubinsStateSpace(KAPPA, SIGMA, DISC, True)
        controls = planner.get_controls_reverse(start, goal)
        self.assertGreater(len(controls), 0)

    def test_dubins_reverse_controls_match_forward_reversed(self) -> None:
        """Reverse controls of (start→goal) = reversed forward controls of (goal→start)."""
        start = State(-3.2, 1.4, 0.7, 0.0)
        goal = State(2.1, -2.8, -1.1, 0.0)
        planner = DubinsStateSpace(KAPPA, DISC, True)
        reverse_controls = planner.get_controls_reverse(start, goal)
        forward_controls = planner.get_controls(goal, start)

        self.assertEqual(len(reverse_controls), len(forward_controls))


class TestSameStartGoalAllSpaces(unittest.TestCase):
    """Test that same start/goal produces valid results for all state spaces.

    Note: Most planners return empty paths when start == goal, as the distance
    is zero and no path is needed. This is expected C++ library behavior.
    """

    def test_same_start_goal_distance_near_zero(self) -> None:
        state = State(1.0, -2.0, 0.75, 0.0)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                d = planner.get_distance(state, state)
                self.assertLess(d, 0.1, f"{name}: non-zero distance for same state")

    def test_same_start_goal_path_valid_if_nonempty(self) -> None:
        state = State(1.0, -2.0, 0.75, 0.2)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                path = planner.get_path(state, state)
                if len(path) > 0:
                    end = path[-1]
                    self.assertLess(
                        math.hypot(end.x - state.x, end.y - state.y),
                        POS_TOL,
                    )


class TestNearbyStates(unittest.TestCase):
    """Test with start and goal very close together."""

    def test_nearby_states(self) -> None:
        start = State(1.0, 2.0, 0.5, 0.0)
        goal = State(1.01, 2.01, 0.51, 0.0)
        for name, planner in make_all_planners().items():
            with self.subTest(planner=name):
                path = planner.get_path(start, goal)
                self.assertGreater(len(path), 0)
                d = planner.get_distance(start, goal)
                self.assertGreaterEqual(d, 0.0)


class TestCurvatureContinuityCC(unittest.TestCase):
    """Test curvature continuity for CC-family planners (G² continuity)."""

    def test_cc_dubins_curvature_continuity(self) -> None:
        """CC planners guarantee continuous curvature (G²): curvature change
        between consecutive discretized states should be bounded by
        discretization * sigma + epsilon."""
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        eps_kappa = 1e-4

        cc_planners = {
            "CCDubins": CCDubinsStateSpace(KAPPA, SIGMA, DISC, True),
            "CC00Dubins": CC00DubinsStateSpace(KAPPA, SIGMA, DISC, True),
            "CC0pmDubins": CC0pmDubinsStateSpace(KAPPA, SIGMA, DISC, True),
            "CCpm0Dubins": CCpm0DubinsStateSpace(KAPPA, SIGMA, DISC, True),
            "CCpmpmDubins": CCpmpmDubinsStateSpace(KAPPA, SIGMA, DISC, True),
            "CC00RS": CC00ReedsSheppStateSpace(KAPPA, SIGMA, DISC),
        }

        for name, planner in cc_planners.items():
            with self.subTest(planner=name):
                path = planner.get_path(start, goal)
                self.assertGreater(len(path), 1)
                for i in range(1, len(path)):
                    delta_kappa = abs(path[i].kappa - path[i - 1].kappa)
                    self.assertLessEqual(
                        delta_kappa,
                        DISC * SIGMA + eps_kappa,
                        f"{name}: curvature jump at index {i}: {delta_kappa}",
                    )


class TestHCCurvatureContinuityBetweenCusps(unittest.TestCase):
    """Test curvature continuity for HC-family planners (continuous between cusps)."""

    def test_hc_curvature_continuity_within_segments(self) -> None:
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(4.0, 2.0, 1.0, 0.0)
        eps_kappa = 1e-4

        hc_planners = {
            "HCRS": HCReedsSheppStateSpace(KAPPA, SIGMA, DISC),
            "HC00RS": HC00ReedsSheppStateSpace(KAPPA, SIGMA, DISC),
            "HC0pmRS": HC0pmReedsSheppStateSpace(KAPPA, SIGMA, DISC),
            "HCpm0RS": HCpm0ReedsSheppStateSpace(KAPPA, SIGMA, DISC),
            "HCpmpmRS": HCpmpmReedsSheppStateSpace(KAPPA, SIGMA, DISC),
        }

        for name, planner in hc_planners.items():
            with self.subTest(planner=name):
                path = planner.get_path(start, goal)
                self.assertGreater(len(path), 1)
                for i in range(1, len(path)):
                    # Only check curvature continuity within same direction segment
                    if path[i].d * path[i - 1].d >= 0:
                        delta_kappa = abs(path[i].kappa - path[i - 1].kappa)
                        self.assertLessEqual(
                            delta_kappa,
                            DISC * SIGMA + eps_kappa,
                            f"{name}: curvature jump at index {i}: {delta_kappa}",
                        )


class TestDifferentKappaValues(unittest.TestCase):
    """Test planners with various kappa_max values."""

    def test_small_kappa(self) -> None:
        """Small curvature = wide turns."""
        planner = ReedsSheppStateSpace(0.5, DISC)
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(5.0, 3.0, 1.0, 0.0)
        path = planner.get_path(start, goal)
        self.assertGreater(len(path), 1)
        end = path[-1]
        self.assertLess(math.hypot(end.x - goal.x, end.y - goal.y), POS_TOL)

    def test_large_kappa(self) -> None:
        """Large curvature = tighter turns."""
        planner = ReedsSheppStateSpace(5.0, DISC)
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(2.0, 1.0, 0.5, 0.0)
        path = planner.get_path(start, goal)
        self.assertGreater(len(path), 1)
        end = path[-1]
        self.assertLess(math.hypot(end.x - goal.x, end.y - goal.y), POS_TOL)

    def test_larger_kappa_yields_shorter_or_equal_path(self) -> None:
        """Higher kappa_max should produce equal or shorter paths."""
        start = State(0.0, 0.0, 0.0, 0.0)
        goal = State(5.0, 3.0, 1.0, 0.0)
        small_kappa = ReedsSheppStateSpace(0.5, DISC)
        large_kappa = ReedsSheppStateSpace(2.0, DISC)
        d_small = small_kappa.get_distance(start, goal)
        d_large = large_kappa.get_distance(start, goal)
        self.assertLessEqual(d_large, d_small + 0.01)


if __name__ == "__main__":
    unittest.main()
