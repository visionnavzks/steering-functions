"""
Comprehensive consistency tests: Python steering functions vs C++ reference.

The C++ ``generate_test_data`` tool produces a JSON file with distances,
controls and path endpoints for 50 random start/goal pairs across **all**
19 state-space variants.  This test suite loads that data and checks that
the Python port produces matching results within floating-point tolerances.

Run from the repository root::

    uv run python -m pytest python_rewrite/tests/test_cpp_consistency.py -v

or::

    uv run python -m unittest python_rewrite.tests.test_cpp_consistency -v
"""

from __future__ import annotations

import json
import math
import os
import sys
import unittest
from pathlib import Path

# ---------------------------------------------------------------------------
# Make sure the Python steering_functions package is importable.
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_PY_ROOT = _REPO_ROOT / "python_rewrite"
if str(_PY_ROOT) not in sys.path:
    sys.path.insert(0, str(_PY_ROOT))

from steering_functions.state import State, Control  # noqa: E402
from steering_functions.dubins_state_space import DubinsStateSpace  # noqa: E402
from steering_functions.reeds_shepp_state_space import ReedsSheppStateSpace  # noqa: E402
from steering_functions.cc_dubins_state_space import (  # noqa: E402
    CC00_Dubins_State_Space,
    CC0pm_Dubins_State_Space,
    CCpm0_Dubins_State_Space,
    CCpmpm_Dubins_State_Space,
    CC_Dubins_State_Space,
)
from steering_functions.cc_reeds_shepp_state_space import (  # noqa: E402
    CC00_Reeds_Shepp_State_Space,
)
from steering_functions.hc_reeds_shepp_state_space import (  # noqa: E402
    HC00_Reeds_Shepp_State_Space,
    HC0pm_Reeds_Shepp_State_Space,
    HCpm0_Reeds_Shepp_State_Space,
    HCpmpm_Reeds_Shepp_State_Space,
    HC_Reeds_Shepp_State_Space,
)

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------
EPS_DISTANCE = 0.05       # positional error (m)
EPS_YAW = 0.05            # heading error (rad)
EPS_KAPPA_ABS = 1e-4      # absolute curvature error
EPS_CONTROL_LEN = 0.05    # total control-length mismatch (m)
EPS_DISTANCE_REL = 0.05   # 5 % relative distance mismatch

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
TEST_DATA_PATH = Path(__file__).resolve().parent / "test_data.json"


def _load_test_data() -> dict:
    """Load and cache the C++ reference data."""
    if not TEST_DATA_PATH.exists():
        raise FileNotFoundError(
            f"C++ reference data not found at {TEST_DATA_PATH}.\n"
            "Build and run generate_test_data first:\n"
            "  BUILD_TESTING=ON ./build.sh\n"
            "  ./build/src/steering_functions/generate_test_data python/tests/test_data.json"
        )
    with open(TEST_DATA_PATH) as f:
        return json.load(f)


def _make_state(d: dict) -> State:
    """Convert a JSON state dict to a Python State."""
    return State(
        x=d["x"],
        y=d["y"],
        theta=d["theta"],
        kappa=d.get("kappa", 0.0),
        d=d.get("d", 0.0),
    )


def _pos_distance(s1, s2) -> float:
    """Euclidean distance between two state-like objects."""
    dx = s1.x - s2.x if hasattr(s1, "x") else s1["x"] - s2["x"]
    dy = s1.y - s2.y if hasattr(s1, "y") else s1["y"] - s2["y"]
    return math.sqrt(dx * dx + dy * dy)


def _controls_length(controls) -> float:
    """Sum of |delta_s| over a list of Control objects."""
    return sum(abs(c.delta_s) for c in controls)


def _controls_length_json(controls_json) -> float:
    """Sum of |delta_s| from JSON control dicts."""
    return sum(abs(c["delta_s"]) for c in controls_json)


# ---------------------------------------------------------------------------
# Factory: name -> Python state space constructor
# ---------------------------------------------------------------------------
def _build_state_spaces(kappa: float, sigma: float, disc: float) -> dict:
    """Return a dict mapping C++ state-space name to Python instance."""
    return {
        "dubins_forwards": DubinsStateSpace(kappa, disc, True),
        "dubins_backwards": DubinsStateSpace(kappa, disc, False),
        "reeds_shepp": ReedsSheppStateSpace(kappa, disc),
        "cc00_dubins_forwards": CC00_Dubins_State_Space(kappa, sigma, disc, True),
        "cc00_dubins_backwards": CC00_Dubins_State_Space(kappa, sigma, disc, False),
        "cc0pm_dubins_forwards": CC0pm_Dubins_State_Space(kappa, sigma, disc, True),
        "cc0pm_dubins_backwards": CC0pm_Dubins_State_Space(kappa, sigma, disc, False),
        "ccpm0_dubins_forwards": CCpm0_Dubins_State_Space(kappa, sigma, disc, True),
        "ccpm0_dubins_backwards": CCpm0_Dubins_State_Space(kappa, sigma, disc, False),
        "ccpmpm_dubins_forwards": CCpmpm_Dubins_State_Space(kappa, sigma, disc, True),
        "ccpmpm_dubins_backwards": CCpmpm_Dubins_State_Space(kappa, sigma, disc, False),
        "cc_dubins_forwards": CC_Dubins_State_Space(kappa, sigma, disc, True),
        "cc_dubins_backwards": CC_Dubins_State_Space(kappa, sigma, disc, False),
        "cc00_reeds_shepp": CC00_Reeds_Shepp_State_Space(kappa, sigma, disc),
        "hc00_reeds_shepp": HC00_Reeds_Shepp_State_Space(kappa, sigma, disc),
        "hc0pm_reeds_shepp": HC0pm_Reeds_Shepp_State_Space(kappa, sigma, disc),
        "hcpm0_reeds_shepp": HCpm0_Reeds_Shepp_State_Space(kappa, sigma, disc),
        "hcpmpm_reeds_shepp": HCpmpm_Reeds_Shepp_State_Space(kappa, sigma, disc),
        "hc_reeds_shepp": HC_Reeds_Shepp_State_Space(kappa, sigma, disc),
    }


# ======================================================================
# Test class
# ======================================================================
class TestCppConsistency(unittest.TestCase):
    """Compare Python and C++ outputs across all state-space variants."""

    _data: dict | None = None
    _spaces: dict | None = None

    @classmethod
    def setUpClass(cls):
        cls._data = _load_test_data()
        kappa = cls._data["kappa"]
        sigma = cls._data["sigma"]
        disc = cls._data["discretization"]
        cls._spaces = _build_state_spaces(kappa, sigma, disc)

    # ------------------------------------------------------------------
    # 1. Distance consistency
    # ------------------------------------------------------------------
    def test_distance_consistency(self):
        """Python get_distance ≈ C++ get_distance for every sample."""
        failures = []
        for si, sample in enumerate(self._data["samples"]):
            start = _make_state(sample["start"])
            goal = _make_state(sample["goal"])
            for ss_name, cpp_result in sample["results"].items():
                py_ss = self._spaces.get(ss_name)
                if py_ss is None:
                    continue
                cpp_dist = cpp_result["distance"]
                try:
                    py_dist = py_ss.get_distance(start, goal)
                except Exception as exc:
                    failures.append(
                        f"  [{ss_name}] sample {si}: Python raised {exc!r}"
                    )
                    continue

                # Allow relative or absolute tolerance
                abs_err = abs(py_dist - cpp_dist)
                rel_err = abs_err / max(abs(cpp_dist), 1e-9)
                if abs_err > EPS_CONTROL_LEN and rel_err > EPS_DISTANCE_REL:
                    failures.append(
                        f"  [{ss_name}] sample {si}: "
                        f"C++={cpp_dist:.6f}  Py={py_dist:.6f}  "
                        f"abs_err={abs_err:.6f}  rel_err={rel_err:.4f}"
                    )

        if failures:
            # Print summary
            msg = (
                f"\n{len(failures)} distance mismatches out of "
                f"{len(self._data['samples']) * len(self._spaces)} checks:\n"
            )
            # Show first 30
            msg += "\n".join(failures[:30])
            if len(failures) > 30:
                msg += f"\n  ... and {len(failures) - 30} more"
            self.fail(msg)

    # ------------------------------------------------------------------
    # 2. Controls total-length consistency
    # ------------------------------------------------------------------
    def test_controls_length_consistency(self):
        """Sum(|delta_s|) of Python controls ≈ C++ controls length."""
        failures = []
        for si, sample in enumerate(self._data["samples"]):
            start = _make_state(sample["start"])
            goal = _make_state(sample["goal"])
            for ss_name, cpp_result in sample["results"].items():
                py_ss = self._spaces.get(ss_name)
                if py_ss is None:
                    continue
                cpp_len = _controls_length_json(cpp_result["controls"])
                try:
                    py_controls = py_ss.get_controls(start, goal)
                    py_len = _controls_length(py_controls)
                except Exception as exc:
                    failures.append(
                        f"  [{ss_name}] sample {si}: Python raised {exc!r}"
                    )
                    continue

                abs_err = abs(py_len - cpp_len)
                rel_err = abs_err / max(abs(cpp_len), 1e-9)
                if abs_err > EPS_CONTROL_LEN and rel_err > EPS_DISTANCE_REL:
                    failures.append(
                        f"  [{ss_name}] sample {si}: "
                        f"C++ ctrl_len={cpp_len:.6f}  Py ctrl_len={py_len:.6f}  "
                        f"abs_err={abs_err:.6f}"
                    )

        if failures:
            msg = (
                f"\n{len(failures)} controls-length mismatches:\n"
            )
            msg += "\n".join(failures[:30])
            if len(failures) > 30:
                msg += f"\n  ... and {len(failures) - 30} more"
            self.fail(msg)

    # ------------------------------------------------------------------
    # 3. Path endpoint consistency (goal reaching)
    # ------------------------------------------------------------------
    def test_path_endpoint_consistency(self):
        """Python path endpoint ≈ C++ path endpoint ≈ goal."""
        failures = []
        for si, sample in enumerate(self._data["samples"]):
            start = _make_state(sample["start"])
            goal = _make_state(sample["goal"])
            for ss_name, cpp_result in sample["results"].items():
                py_ss = self._spaces.get(ss_name)
                if py_ss is None:
                    continue
                cpp_end = cpp_result.get("path_end")
                if cpp_end is None:
                    continue

                try:
                    py_path = py_ss.get_path(start, goal)
                except Exception as exc:
                    failures.append(
                        f"  [{ss_name}] sample {si}: Python get_path raised {exc!r}"
                    )
                    continue

                if not py_path:
                    failures.append(
                        f"  [{ss_name}] sample {si}: Python returned empty path"
                    )
                    continue

                py_end = py_path[-1]

                # Check Python endpoint vs goal
                py_goal_dist = _pos_distance(py_end, goal)
                cpp_goal_dist = math.sqrt(
                    (cpp_end["x"] - goal.x) ** 2
                    + (cpp_end["y"] - goal.y) ** 2
                )

                if py_goal_dist > EPS_DISTANCE:
                    failures.append(
                        f"  [{ss_name}] sample {si}: "
                        f"Py endpoint {py_goal_dist:.6f}m from goal "
                        f"(C++ was {cpp_goal_dist:.6f}m from goal)"
                    )

        if failures:
            msg = (
                f"\n{len(failures)} path-endpoint mismatches:\n"
            )
            msg += "\n".join(failures[:30])
            if len(failures) > 30:
                msg += f"\n  ... and {len(failures) - 30} more"
            self.fail(msg)

    # ------------------------------------------------------------------
    # 4. Number of control segments consistency
    # ------------------------------------------------------------------
    def test_num_controls_consistency(self):
        """Python and C++ produce the same number of control segments."""
        failures = []
        for si, sample in enumerate(self._data["samples"]):
            start = _make_state(sample["start"])
            goal = _make_state(sample["goal"])
            for ss_name, cpp_result in sample["results"].items():
                py_ss = self._spaces.get(ss_name)
                if py_ss is None:
                    continue
                cpp_n = len(cpp_result["controls"])
                try:
                    py_controls = py_ss.get_controls(start, goal)
                    py_n = len(py_controls)
                except Exception as exc:
                    failures.append(
                        f"  [{ss_name}] sample {si}: Python raised {exc!r}"
                    )
                    continue

                if py_n != cpp_n:
                    failures.append(
                        f"  [{ss_name}] sample {si}: "
                        f"C++ has {cpp_n} segments, Python has {py_n}"
                    )

        if failures:
            msg = (
                f"\n{len(failures)} segment-count mismatches:\n"
            )
            msg += "\n".join(failures[:30])
            if len(failures) > 30:
                msg += f"\n  ... and {len(failures) - 30} more"
            self.fail(msg)

    # ------------------------------------------------------------------
    # 5. Per-control segment comparison
    # ------------------------------------------------------------------
    def test_per_control_values(self):
        """Each control segment's delta_s, kappa, sigma ≈ C++ values."""
        failures = []
        for si, sample in enumerate(self._data["samples"]):
            start = _make_state(sample["start"])
            goal = _make_state(sample["goal"])
            for ss_name, cpp_result in sample["results"].items():
                py_ss = self._spaces.get(ss_name)
                if py_ss is None:
                    continue
                cpp_ctrls = cpp_result["controls"]
                try:
                    py_ctrls = py_ss.get_controls(start, goal)
                except Exception:
                    continue  # already caught in other tests

                if len(py_ctrls) != len(cpp_ctrls):
                    continue  # already caught in test_num_controls

                for ci, (py_c, cpp_c) in enumerate(zip(py_ctrls, cpp_ctrls)):
                    ds_err = abs(py_c.delta_s - cpp_c["delta_s"])
                    k_err = abs(py_c.kappa - cpp_c["kappa"])
                    s_err = abs(py_c.sigma - cpp_c["sigma"])
                    if ds_err > EPS_CONTROL_LEN or k_err > 0.01 or s_err > 0.01:
                        failures.append(
                            f"  [{ss_name}] sample {si} ctrl {ci}: "
                            f"delta_s err={ds_err:.6f}  "
                            f"kappa err={k_err:.6f}  "
                            f"sigma err={s_err:.6f}"
                        )
                        break  # one failure per sample/ss is enough

        if failures:
            msg = (
                f"\n{len(failures)} per-control-value mismatches:\n"
            )
            msg += "\n".join(failures[:30])
            if len(failures) > 30:
                msg += f"\n  ... and {len(failures) - 30} more"
            self.fail(msg)

    # ------------------------------------------------------------------
    # 6. Curvature continuity (Python-only sanity)
    # ------------------------------------------------------------------
    def test_curvature_continuity(self):
        """Adjacent path states have curvature change ≤ disc*sigma + eps."""
        disc = self._data["discretization"]
        sigma = self._data["sigma"]
        max_dk = disc * sigma + 1e-4
        failures = []

        # Only test CC/HC variants (Dubins/RS don't guarantee curvature continuity)
        cc_hc_names = [
            n for n in self._spaces
            if n.startswith("cc") or n.startswith("hc")
        ]

        for si, sample in enumerate(self._data["samples"][:20]):  # first 20
            start = _make_state(sample["start"])
            goal = _make_state(sample["goal"])
            for ss_name in cc_hc_names:
                py_ss = self._spaces[ss_name]
                try:
                    path = py_ss.get_path(start, goal)
                except Exception:
                    continue
                if len(path) < 2:
                    continue
                for pi in range(1, len(path)):
                    if path[pi - 1].d * path[pi].d < 0:
                        continue
                    dk = abs(path[pi].kappa - path[pi - 1].kappa)
                    if dk > max_dk:
                        failures.append(
                            f"  [{ss_name}] sample {si} step {pi}: "
                            f"dk={dk:.6f} > max_dk={max_dk:.6f}"
                        )
                        break

        if failures:
            msg = (
                f"\n{len(failures)} curvature-continuity violations:\n"
            )
            msg += "\n".join(failures[:30])
            if len(failures) > 30:
                msg += f"\n  ... and {len(failures) - 30} more"
            self.fail(msg)

    # ------------------------------------------------------------------
    # 7. Path length consistency (distance ≈ sum of point-to-point)
    # ------------------------------------------------------------------
    def test_path_length_vs_distance(self):
        """The discretized path length ≈ the reported distance."""
        failures = []
        for si, sample in enumerate(self._data["samples"][:20]):
            start = _make_state(sample["start"])
            goal = _make_state(sample["goal"])
            for ss_name, py_ss in self._spaces.items():
                try:
                    dist = py_ss.get_distance(start, goal)
                    path = py_ss.get_path(start, goal)
                except Exception:
                    continue
                if len(path) < 2:
                    continue
                # Compute path length by summing segment lengths
                path_len = 0.0
                for pi in range(1, len(path)):
                    path_len += _pos_distance(path[pi], path[pi - 1])

                err = abs(dist - path_len)
                if err > EPS_DISTANCE and err / max(dist, 1e-9) > 0.05:
                    failures.append(
                        f"  [{ss_name}] sample {si}: "
                        f"distance={dist:.6f}  path_len={path_len:.6f}  "
                        f"err={err:.6f}"
                    )

        if failures:
            msg = (
                f"\n{len(failures)} path-length vs distance mismatches:\n"
            )
            msg += "\n".join(failures[:30])
            if len(failures) > 30:
                msg += f"\n  ... and {len(failures) - 30} more"
            self.fail(msg)

    # ------------------------------------------------------------------
    # 8. Symmetry: get_distance(a,b) for RS-type should ≈ get_distance(b,a)
    # ------------------------------------------------------------------
    def test_distance_symmetry(self):
        """RS-type state spaces: distance(a,b) ≈ distance(b,a)."""
        symmetric_names = [
            "reeds_shepp",
            "cc00_reeds_shepp",
            "hc00_reeds_shepp",
            "hc_reeds_shepp",
        ]
        failures = []
        for si, sample in enumerate(self._data["samples"][:20]):
            start = _make_state(sample["start"])
            goal = _make_state(sample["goal"])
            for ss_name in symmetric_names:
                py_ss = self._spaces.get(ss_name)
                if py_ss is None:
                    continue
                try:
                    d_ab = py_ss.get_distance(start, goal)
                    d_ba = py_ss.get_distance(goal, start)
                except Exception:
                    continue
                err = abs(d_ab - d_ba)
                if err > EPS_CONTROL_LEN:
                    failures.append(
                        f"  [{ss_name}] sample {si}: "
                        f"d(a,b)={d_ab:.6f}  d(b,a)={d_ba:.6f}  err={err:.6f}"
                    )

        if failures:
            msg = f"\n{len(failures)} symmetry violations:\n"
            msg += "\n".join(failures[:30])
            self.fail(msg)


# ======================================================================
# Per-state-space detailed tests (individual subTest for granular output)
# ======================================================================
class TestPerStateSpace(unittest.TestCase):
    """Run per-state-space distance checks with subTest for clarity."""

    _data: dict | None = None
    _spaces: dict | None = None

    @classmethod
    def setUpClass(cls):
        cls._data = _load_test_data()
        kappa = cls._data["kappa"]
        sigma = cls._data["sigma"]
        disc = cls._data["discretization"]
        cls._spaces = _build_state_spaces(kappa, sigma, disc)

    def _run_for_ss(self, ss_name: str):
        py_ss = self._spaces[ss_name]
        n_ok = 0
        n_fail = 0
        for si, sample in enumerate(self._data["samples"]):
            start = _make_state(sample["start"])
            goal = _make_state(sample["goal"])
            cpp_result = sample["results"].get(ss_name)
            if cpp_result is None:
                continue
            cpp_dist = cpp_result["distance"]
            try:
                py_dist = py_ss.get_distance(start, goal)
            except Exception:
                n_fail += 1
                continue
            abs_err = abs(py_dist - cpp_dist)
            rel_err = abs_err / max(abs(cpp_dist), 1e-9)
            if abs_err > EPS_CONTROL_LEN and rel_err > EPS_DISTANCE_REL:
                n_fail += 1
            else:
                n_ok += 1
        return n_ok, n_fail

    def test_dubins_forwards(self):
        ok, fail = self._run_for_ss("dubins_forwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_dubins_backwards(self):
        ok, fail = self._run_for_ss("dubins_backwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_reeds_shepp(self):
        ok, fail = self._run_for_ss("reeds_shepp")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_cc00_dubins_forwards(self):
        ok, fail = self._run_for_ss("cc00_dubins_forwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_cc00_dubins_backwards(self):
        ok, fail = self._run_for_ss("cc00_dubins_backwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_cc0pm_dubins_forwards(self):
        ok, fail = self._run_for_ss("cc0pm_dubins_forwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_cc0pm_dubins_backwards(self):
        ok, fail = self._run_for_ss("cc0pm_dubins_backwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_ccpm0_dubins_forwards(self):
        ok, fail = self._run_for_ss("ccpm0_dubins_forwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_ccpm0_dubins_backwards(self):
        ok, fail = self._run_for_ss("ccpm0_dubins_backwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_ccpmpm_dubins_forwards(self):
        ok, fail = self._run_for_ss("ccpmpm_dubins_forwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_ccpmpm_dubins_backwards(self):
        ok, fail = self._run_for_ss("ccpmpm_dubins_backwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_cc_dubins_forwards(self):
        ok, fail = self._run_for_ss("cc_dubins_forwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_cc_dubins_backwards(self):
        ok, fail = self._run_for_ss("cc_dubins_backwards")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_cc00_reeds_shepp(self):
        ok, fail = self._run_for_ss("cc00_reeds_shepp")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_hc00_reeds_shepp(self):
        ok, fail = self._run_for_ss("hc00_reeds_shepp")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_hc0pm_reeds_shepp(self):
        ok, fail = self._run_for_ss("hc0pm_reeds_shepp")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_hcpm0_reeds_shepp(self):
        ok, fail = self._run_for_ss("hcpm0_reeds_shepp")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_hcpmpm_reeds_shepp(self):
        ok, fail = self._run_for_ss("hcpmpm_reeds_shepp")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")

    def test_hc_reeds_shepp(self):
        ok, fail = self._run_for_ss("hc_reeds_shepp")
        self.assertEqual(fail, 0, f"{fail}/{ok+fail} samples failed")


if __name__ == "__main__":
    unittest.main()
