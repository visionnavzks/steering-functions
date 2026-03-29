#!/usr/bin/env python3
"""
Interactive visualization of steering function paths.

Usage::

    python3 python/visualize.py

Interaction
-----------
* **Right panel – Method**: choose the steering algorithm.
* **Right panel – Mode**: *Best* (shortest path only) or *All* (every candidate path).
* **kappa_max / sigma_max sliders**: tune robot constraints.
* **State text boxes**: edit start/goal x, y, and heading directly.

Keyboard shortcuts
------------------
* ``r`` – reset to default start/goal
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import List

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, RadioButtons, Slider, TextBox

# ---------------------------------------------------------------------------
# Make steering_functions importable when run from repo root
# ---------------------------------------------------------------------------
_PY_ROOT = Path(__file__).resolve().parent
if str(_PY_ROOT) not in sys.path:
    sys.path.insert(0, str(_PY_ROOT))

from steering_functions.state import State  # noqa: E402
from steering_functions.steering_path import PathType, SteeringPath  # noqa: E402

# ---------------------------------------------------------------------------
# Method registry:  display name -> (PathType, needs_sigma)
# ---------------------------------------------------------------------------
METHODS: dict[str, tuple[PathType, bool]] = {
    "Dubins":       (PathType.DUBINS,      False),
    "CC±±-Dubins":  (PathType.CCPMPM_DUBINS, True),
    "CC±0-Dubins":  (PathType.CCPM0_DUBINS,  True),
    "CC0±-Dubins":  (PathType.CC0PM_DUBINS,  True),
    "CC00-Dubins":  (PathType.CC00_DUBINS,   True),
    "CC-Dubins":    (PathType.CC_DUBINS,     True),
    "RS":           (PathType.RS,            False),
    "CC00-RS":      (PathType.CC00_RS,       True),
    "HC±±-RS":      (PathType.HCPMPM_RS,     True),
    "HC±0-RS":      (PathType.HCPM0_RS,      True),
    "HC0±-RS":      (PathType.HC0PM_RS,      True),
    "HC00-RS":      (PathType.HC00_RS,       True),
    "HC-RS":        (PathType.HC_RS,         True),
}

METHOD_NAMES = list(METHODS.keys())

# Colors for "All" paths mode (cycled if more paths than colors)
_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
]

# Default robot / algorithm parameters
_DEFAULT_KAPPA_MAX = 1.0   # 1/m
_DEFAULT_SIGMA_MAX = 1.0   # 1/m²
_DISCRETIZATION   = 0.02   # m


# ---------------------------------------------------------------------------
# Helper: extract (x, y) arrays from a list of State objects
# ---------------------------------------------------------------------------
def _path_xy(path: List[State]):
    xs = [s.x for s in path]
    ys = [s.y for s in path]
    return xs, ys


def _total_length(controls) -> float:
    """Compute total arc length from a control sequence."""
    return sum(abs(c.delta_s) for c in controls)


def _state_arrow(ax, state: State, color: str, length: float = 0.25, zorder: int = 5):
    """Draw a filled arrow showing position and heading of a state."""
    dx = length * math.cos(state.theta)
    dy = length * math.sin(state.theta)
    ax.annotate(
        "",
        xy=(state.x + dx, state.y + dy),
        xytext=(state.x, state.y),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=2.0),
        zorder=zorder,
    )
    ax.plot(state.x, state.y, "o", color=color, markersize=9, zorder=zorder)


# ---------------------------------------------------------------------------
# Main visualizer class
# ---------------------------------------------------------------------------
class SteeringVisualizer:
    def __init__(self) -> None:
        # Default start / goal states
        self.start = State(x=-3.0, y=0.0,  theta=0.0)
        self.goal  = State(x= 3.0, y=2.0,  theta=math.pi / 4)

        self.method     : str   = METHOD_NAMES[0]   # "Dubins"
        self.show_all   : bool  = False              # best vs all
        self.kappa_max  : float = _DEFAULT_KAPPA_MAX
        self.sigma_max  : float = _DEFAULT_SIGMA_MAX

        self._updating_textboxes: bool = False  # guard against recursive TextBox callbacks

        self._setup_figure()
        self._connect_events()
        self._compute_and_draw()

    # ------------------------------------------------------------------
    # Figure / axes / widget construction
    # ------------------------------------------------------------------
    def _setup_figure(self) -> None:
        self.fig = plt.figure(figsize=(16, 9))
        self.fig.patch.set_facecolor("#f8f8f8")
        self.fig.canvas.manager.set_window_title("Steering Functions Visualizer")

        # ── Main plot axes ─────────────────────────────────────────────
        self.ax = self.fig.add_axes([0.03, 0.12, 0.60, 0.85])
        self.ax.set_aspect("equal")
        self.ax.set_facecolor("#ffffff")
        self.ax.grid(True, linestyle="--", alpha=0.4)
        self.ax.set_title("Steering Function Path Visualizer", fontsize=13, pad=8)
        self.ax.set_xlabel("x [m]")
        self.ax.set_ylabel("y [m]")

        # ── Method radio buttons ───────────────────────────────────────
        ax_radio = self.fig.add_axes(
            [0.66, 0.30, 0.31, 0.65],
            facecolor="#f0f0f0",
        )
        self.radio_method = RadioButtons(ax_radio, METHOD_NAMES, active=0)
        ax_radio.set_title("Method", fontsize=10, pad=4)
        for lbl in self.radio_method.labels:
            lbl.set_fontsize(9)
        self.radio_method.on_clicked(self._on_method_change)

        # ── Mode radio (Best / All) ────────────────────────────────────
        ax_mode = self.fig.add_axes([0.66, 0.21, 0.31, 0.07], facecolor="#f0f0f0")
        self.radio_mode = RadioButtons(ax_mode, ["Best path", "All paths"], active=0)
        ax_mode.set_title("Display mode", fontsize=10, pad=2)
        for lbl in self.radio_mode.labels:
            lbl.set_fontsize(9)
        self.radio_mode.on_clicked(self._on_mode_change)

        # ── kappa_max slider ──────────────────────────────────────────
        ax_kappa = self.fig.add_axes([0.66, 0.145, 0.31, 0.025])
        self.slider_kappa = Slider(
            ax_kappa, "κ_max", 0.1, 5.0,
            valinit=_DEFAULT_KAPPA_MAX, valstep=0.05, color="#4a90d9",
        )
        self.slider_kappa.label.set_fontsize(9)
        self.slider_kappa.on_changed(self._on_kappa_change)

        # ── sigma_max slider ──────────────────────────────────────────
        ax_sigma = self.fig.add_axes([0.66, 0.105, 0.31, 0.025])
        self.slider_sigma = Slider(
            ax_sigma, "σ_max", 0.1, 5.0,
            valinit=_DEFAULT_SIGMA_MAX, valstep=0.05, color="#4ab04a",
        )
        self.slider_sigma.label.set_fontsize(9)
        self.slider_sigma.on_changed(self._on_sigma_change)

        # ── Reset button ───────────────────────────────────────────────
        ax_reset = self.fig.add_axes([0.66, 0.008, 0.145, 0.035])
        self.btn_reset = Button(ax_reset, "↺ Reset", color="#e8e8e8", hovercolor="#d0d0d0")
        self.btn_reset.label.set_fontsize(9)
        self.btn_reset.on_clicked(self._on_reset)

        # ── Info text at bottom of main axes ──────────────────────────
        self._info_text = self.ax.text(
            0.01, 0.01, "",
            transform=self.ax.transAxes,
            fontsize=8, verticalalignment="bottom",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8),
        )

        # ── State input text boxes ─────────────────────────────────────
        self.fig.text(0.040, 0.079, "Start:", fontsize=9, fontweight="bold",
                      color="#2ca02c")
        self.fig.text(0.290, 0.079, "Goal:",  fontsize=9, fontweight="bold",
                      color="#d62728")
        self.fig.text(
            0.505, 0.079,
            "Edit state values in the text boxes  "
            "|  ↵ Enter to apply  |  Keys: [r] Reset  [R] Randomize",
            fontsize=7.5, color="#666666",
        )

        # Start TextBoxes  (x, y, θ in degrees)
        ax_sx = self.fig.add_axes([0.040, 0.022, 0.058, 0.036])
        ax_sy = self.fig.add_axes([0.106, 0.022, 0.058, 0.036])
        ax_st = self.fig.add_axes([0.172, 0.022, 0.066, 0.036])
        self.tb_sx = TextBox(ax_sx, "x",  initial=f"{self.start.x:.3f}")
        self.tb_sy = TextBox(ax_sy, "y",  initial=f"{self.start.y:.3f}")
        self.tb_st = TextBox(ax_st, "θ°", initial=f"{math.degrees(self.start.theta):.1f}")
        for tb in (self.tb_sx, self.tb_sy, self.tb_st):
            tb.label.set_fontsize(9)

        # Goal TextBoxes
        ax_gx = self.fig.add_axes([0.290, 0.022, 0.058, 0.036])
        ax_gy = self.fig.add_axes([0.356, 0.022, 0.058, 0.036])
        ax_gt = self.fig.add_axes([0.422, 0.022, 0.066, 0.036])
        self.tb_gx = TextBox(ax_gx, "x",  initial=f"{self.goal.x:.3f}")
        self.tb_gy = TextBox(ax_gy, "y",  initial=f"{self.goal.y:.3f}")
        self.tb_gt = TextBox(ax_gt, "θ°", initial=f"{math.degrees(self.goal.theta):.1f}")
        for tb in (self.tb_gx, self.tb_gy, self.tb_gt):
            tb.label.set_fontsize(9)

        # Randomize button
        ax_rand = self.fig.add_axes([0.505, 0.022, 0.110, 0.036])
        self.btn_rand = Button(ax_rand, "⚄ Randomize", color="#d4e8ff", hovercolor="#aaccff")
        self.btn_rand.label.set_fontsize(9)
        self.btn_rand.on_clicked(self._on_randomize)

        # Connect TextBox submit callbacks (fired on Enter)
        self.tb_sx.on_submit(lambda v: self._on_state_text("start", "x",     v))
        self.tb_sy.on_submit(lambda v: self._on_state_text("start", "y",     v))
        self.tb_st.on_submit(lambda v: self._on_state_text("start", "theta", v))
        self.tb_gx.on_submit(lambda v: self._on_state_text("goal",  "x",     v))
        self.tb_gy.on_submit(lambda v: self._on_state_text("goal",  "y",     v))
        self.tb_gt.on_submit(lambda v: self._on_state_text("goal",  "theta", v))

    # ------------------------------------------------------------------
    # Event connection
    # ------------------------------------------------------------------
    def _connect_events(self) -> None:
        self.fig.canvas.mpl_connect("key_press_event",      self._on_key)

    # ------------------------------------------------------------------
    # Widget callbacks
    # ------------------------------------------------------------------
    def _on_method_change(self, label: str) -> None:
        self.method = label
        self._compute_and_draw()

    def _on_mode_change(self, label: str) -> None:
        self.show_all = (label == "All paths")
        self._compute_and_draw()

    def _on_kappa_change(self, val: float) -> None:
        self.kappa_max = val
        self._compute_and_draw()

    def _on_sigma_change(self, val: float) -> None:
        self.sigma_max = val
        self._compute_and_draw()

    def _on_reset(self, _event) -> None:
        self.start = State(x=-3.0, y=0.0, theta=0.0)
        self.goal  = State(x= 3.0, y=2.0, theta=math.pi / 4)
        self._sync_textboxes()
        self._compute_and_draw()

    # ------------------------------------------------------------------
    # Keyboard handlers
    # ------------------------------------------------------------------
    def _on_key(self, event) -> None:
        if event.key == "r":
            self._on_reset(None)
        elif event.key == "R":
            self._on_randomize(None)

    # ------------------------------------------------------------------
    # Core: compute paths and render
    # ------------------------------------------------------------------
    def _compute_and_draw(self) -> None:
        self.ax.cla()
        self.ax.set_aspect("equal")
        self.ax.set_facecolor("#ffffff")
        self.ax.grid(True, linestyle="--", alpha=0.4)
        self.ax.set_xlabel("x [m]")
        self.ax.set_ylabel("y [m]")

        path_type, needs_sigma = METHODS[self.method]
        sigma = self.sigma_max if needs_sigma else max(self.sigma_max, 1e-3)

        info_lines: list[str] = []
        info_lines.append(
            f"Method: {self.method}   "
            f"κ_max={self.kappa_max:.2f} m⁻¹   "
            f"σ_max={sigma:.2f} m⁻²"
        )

        try:
            planner = SteeringPath(
                path_type,
                kappa_max=self.kappa_max,
                sigma_max=sigma,
                discretization=_DISCRETIZATION,
            )

            if self.show_all:
                self._draw_all_paths(planner, info_lines)
            else:
                self._draw_best_path(planner, info_lines)

        except Exception as exc:  # noqa: BLE001
            self.ax.set_title(f"⚠  {exc}", color="red", fontsize=10)
            info_lines.append(f"Error: {exc}")

        # Draw start and goal states on top
        _state_arrow(self.ax, self.start, color="#2ca02c", length=0.30, zorder=10)
        _state_arrow(self.ax, self.goal,  color="#d62728", length=0.30, zorder=10)

        # Legend entries for start / goal
        start_patch = mpatches.Patch(color="#2ca02c", label="Start")
        goal_patch  = mpatches.Patch(color="#d62728", label="Goal")
        handles     = [start_patch, goal_patch]

        # Title and legend
        mode_str = "All paths" if self.show_all else "Best path"
        self.ax.set_title(
            f"{self.method}  –  {mode_str}",
            fontsize=12, pad=6,
        )

        self.ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.8)

        # Auto-zoom with some margin around start/goal
        self._auto_zoom()

        # Info text
        self._info_text = self.ax.text(
            0.01, 0.01, "\n".join(info_lines),
            transform=self.ax.transAxes,
            fontsize=8, verticalalignment="bottom",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.85),
        )

        self.fig.canvas.draw_idle()

    # ------------------------------------------------------------------
    def _draw_best_path(self, planner: SteeringPath, info_lines: list[str]) -> None:
        controls = planner.compute_shortest_control_sequence(self.start, self.goal)
        path     = planner.compute_shortest_path(self.start, self.goal)

        if not path:
            info_lines.append("No path found.")
            return

        xs, ys = _path_xy(path)
        length = _total_length(controls)

        self.ax.plot(xs, ys, color="#1f77b4", lw=2.2, zorder=3, label=f"Path  {length:.3f} m")
        # Add length annotation
        mid = len(path) // 2
        self.ax.annotate(
            f"{length:.2f} m",
            xy=(path[mid].x, path[mid].y),
            xytext=(10, 10), textcoords="offset points",
            fontsize=8, color="#1f77b4",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7),
        )
        info_lines.append(f"Path length: {length:.4f} m")

    # ------------------------------------------------------------------
    def _draw_all_paths(self, planner: SteeringPath, info_lines: list[str]) -> None:
        all_controls = planner.compute_all_control_sequences(self.start, self.goal)
        all_paths    = planner.compute_all_paths(self.start, self.goal)

        if not all_paths:
            info_lines.append("No paths found.")
            return

        # Sort by length (shortest first) so we can highlight the best
        lengths = [_total_length(cs) for cs in all_controls]
        order   = sorted(range(len(lengths)), key=lambda i: lengths[i])

        best_idx   = order[0]
        best_len   = lengths[best_idx]
        info_lines.append(f"Paths found: {len(all_paths)}   |   Best: {best_len:.4f} m")

        for rank, idx in enumerate(order):
            path = all_paths[idx]
            if not path:
                continue
            xs, ys = _path_xy(path)
            color  = _PALETTE[rank % len(_PALETTE)]
            lw     = 2.5 if idx == best_idx else 1.2
            alpha  = 1.0 if idx == best_idx else 0.55
            zorder = 5 if idx == best_idx else 2
            label  = f"#{rank+1}  {lengths[idx]:.2f} m" + (" ★" if idx == best_idx else "")
            self.ax.plot(xs, ys, color=color, lw=lw, alpha=alpha, zorder=zorder, label=label)

        self.ax.legend(
            loc="upper right", fontsize=7, framealpha=0.85,
            ncol=2 if len(all_paths) > 8 else 1,
        )

    # ------------------------------------------------------------------
    def _auto_zoom(self) -> None:
        """Set axis limits to comfortably show both states plus any drawn paths."""
        # Collect all line data to determine bounds
        all_x, all_y = [], []
        for line in self.ax.lines:
            xd, yd = line.get_xdata(), line.get_ydata()
            all_x.extend(xd)
            all_y.extend(yd)

        # Always include start and goal positions
        for s in (self.start, self.goal):
            all_x.append(s.x)
            all_y.append(s.y)

        if not all_x:
            return

        margin = 1.5
        x_min, x_max = min(all_x) - margin, max(all_x) + margin
        y_min, y_max = min(all_y) - margin, max(all_y) + margin

        # Keep aspect ratio equal: expand the narrower axis
        cx, cy = (x_min + x_max) / 2, (y_min + y_max) / 2
        half   = max((x_max - x_min) / 2, (y_max - y_min) / 2)
        self.ax.set_xlim(cx - half, cx + half)
        self.ax.set_ylim(cy - half, cy + half)

    # ------------------------------------------------------------------
    def _sync_textboxes(self) -> None:
        """Update all TextBox values to reflect the current start/goal states."""
        self._updating_textboxes = True
        self.tb_sx.set_val(f"{self.start.x:.3f}")
        self.tb_sy.set_val(f"{self.start.y:.3f}")
        self.tb_st.set_val(f"{math.degrees(self.start.theta):.1f}")
        self.tb_gx.set_val(f"{self.goal.x:.3f}")
        self.tb_gy.set_val(f"{self.goal.y:.3f}")
        self.tb_gt.set_val(f"{math.degrees(self.goal.theta):.1f}")
        self._updating_textboxes = False

    def _on_state_text(self, state_name: str, field: str, value: str) -> None:
        """Apply a typed value from a TextBox (triggered by Enter)."""
        if self._updating_textboxes:
            return
        try:
            val = float(value)
        except ValueError:
            self._sync_textboxes()   # revert display
            return
        if state_name == "start":
            s = State(x=self.start.x, y=self.start.y, theta=self.start.theta)
        else:
            s = State(x=self.goal.x,  y=self.goal.y,  theta=self.goal.theta)
        if field == "x":
            s.x = val
        elif field == "y":
            s.y = val
        else:   # degrees → radians
            s.theta = math.radians(val)
        if state_name == "start":
            self.start = s
        else:
            self.goal = s
        self._compute_and_draw()

    def _on_randomize(self, _event) -> None:
        """Randomly sample start and goal states within ±5 m."""
        rng = np.random.default_rng()
        R = 5.0
        self.start = State(
            x=float(rng.uniform(-R, R)),
            y=float(rng.uniform(-R, R)),
            theta=float(rng.uniform(-math.pi, math.pi)),
        )
        self.goal = State(
            x=float(rng.uniform(-R, R)),
            y=float(rng.uniform(-R, R)),
            theta=float(rng.uniform(-math.pi, math.pi)),
        )
        self._sync_textboxes()
        self._compute_and_draw()

    # ------------------------------------------------------------------
    def show(self) -> None:
        plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    viz = SteeringVisualizer()
    viz.show()
