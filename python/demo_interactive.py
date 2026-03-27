"""Interactive demo for steering_functions.

Features:
- Dropdown to select the path planning algorithm (PathType)
- "Randomize" button to generate new random start/goal and recompute the path
- The planned path, start pose, and goal pose are drawn on a 2D canvas
"""

import math
import random
import matplotlib

matplotlib.use("TkAgg")  # change to "Qt5Agg" or "Agg" if TkAgg is unavailable

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons

import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from steering_functions import PathType, State, SteeringPath

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
KAPPA_MAX = 1.0       # max curvature  [1/m]
SIGMA_MAX = 1.0       # max curvature rate [1/m²]
DISCRETIZATION = 0.1  # arc-length step [m]

WORLD_RANGE = (-5.0, 5.0)  # x and y range for random poses

ARROW_LEN = 0.4       # length of pose arrows

# PathType choices shown in the radio button (name -> PathType)
METHOD_OPTIONS = {
    "Dubins":       PathType.DUBINS,
    "Reeds-Shepp":  PathType.RS,
    "CC-Dubins":    PathType.CC_DUBINS,
    "CC00-Dubins":  PathType.CC00_DUBINS,
    "CC00-RS":      PathType.CC00_RS,
    "HC-RS":        PathType.HC_RS,
    "HC00-RS":      PathType.HC00_RS,
    "HCPMPM-RS":    PathType.HCPMPM_RS,
}
METHOD_LABELS = list(METHOD_OPTIONS.keys())

DISPLAY_OPTIONS = {
    "Shortest": "shortest",
    "All": "all",
}
DISPLAY_LABELS = list(DISPLAY_OPTIONS.keys())

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def random_state():
    lo, hi = WORLD_RANGE
    x = random.uniform(lo, hi)
    y = random.uniform(lo, hi)
    theta = random.uniform(-math.pi, math.pi)
    return State(x=x, y=y, theta=theta)


def draw_pose(ax, state, color, label):
    """Draw a filled circle + heading arrow for a state."""
    ax.plot(state.x, state.y, "o", color=color, markersize=8, zorder=5)
    dx = ARROW_LEN * math.cos(state.theta)
    dy = ARROW_LEN * math.sin(state.theta)
    ax.annotate(
        "", xy=(state.x + dx, state.y + dy), xytext=(state.x, state.y),
        arrowprops=dict(arrowstyle="->", color=color, lw=2),
        zorder=6,
    )
    ax.text(state.x + 0.1, state.y + 0.15, label, color=color, fontsize=9, zorder=7)


def path_length(states):
    if not states:
        return float("nan")
    return states[-1].s if hasattr(states[-1], "s") else float("nan")


def draw_path(ax, states, color, label=None, linewidth=2, alpha=1.0, zorder=3):
    xs = [state.x for state in states]
    ys = [state.y for state in states]
    ax.plot(xs, ys, "-", color=color, linewidth=linewidth, alpha=alpha, label=label, zorder=zorder)


def compute_and_draw(ax, path_type, start, goal, display_mode):
    """Run the planner and render the result on *ax*."""
    ax.cla()
    ax.set_xlim(WORLD_RANGE[0] - 1, WORLD_RANGE[1] + 1)
    ax.set_ylim(WORLD_RANGE[0] - 1, WORLD_RANGE[1] + 1)
    ax.set_aspect("equal")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.set_title(f"Method: {path_type.name}", fontsize=11)

    # Draw start and goal poses
    draw_pose(ax, start, "green", "Start")
    draw_pose(ax, goal,  "red",   "Goal")

    try:
        planner = SteeringPath(path_type, KAPPA_MAX, SIGMA_MAX, DISCRETIZATION)
        if display_mode == "all":
            all_paths = planner.computeAllPaths(start, goal)

            if all_paths:
                colors = [plt.cm.tab10(index % 10) for index in range(max(len(all_paths), 1))]
                shortest_index = 0
                shortest_length = float("inf")

                for index, states in enumerate(all_paths):
                    if not states:
                        continue
                    length = path_length(states)
                    if not math.isnan(length) and length < shortest_length:
                        shortest_index = index
                        shortest_length = length

                for index, states in enumerate(all_paths):
                    if not states:
                        continue
                    is_shortest = index == shortest_index
                    draw_path(
                        ax,
                        states,
                        colors[index % len(colors)],
                        label="Shortest" if is_shortest else None,
                        linewidth=2.4 if is_shortest else 1.6,
                        alpha=0.95 if is_shortest else 0.45,
                        zorder=4 if is_shortest else 3,
                    )

                summary = f"Paths: {len([states for states in all_paths if states])}"
                if shortest_length != float("inf"):
                    summary += f"\nShortest: {shortest_length:.2f} m"
                ax.text(
                    0.02, 0.97, summary,
                    transform=ax.transAxes, va="top", fontsize=9, color="steelblue",
                )
            else:
                ax.text(0.5, 0.5, "No path found", transform=ax.transAxes,
                        ha="center", va="center", fontsize=12, color="gray")
        else:
            path_states = planner.computeShortestPath(start, goal)

            if path_states:
                draw_path(ax, path_states, "steelblue", label="Path")
                length = path_length(path_states)
                ax.text(
                    0.02, 0.97, f"Length: {length:.2f} m",
                    transform=ax.transAxes, va="top", fontsize=9, color="steelblue",
                )
            else:
                ax.text(0.5, 0.5, "No path found", transform=ax.transAxes,
                        ha="center", va="center", fontsize=12, color="gray")
    except Exception as exc:
        ax.text(0.5, 0.5, f"Error:\n{exc}", transform=ax.transAxes,
                ha="center", va="center", fontsize=9, color="red", wrap=True)

    ax.figure.canvas.draw_idle()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    fig = plt.figure(figsize=(11, 7))
    fig.suptitle("Steering Functions Interactive Demo", fontsize=13, fontweight="bold")

    # Layout: left panel (controls) + right panel (path canvas)
    ax_radio = fig.add_axes([0.01, 0.48, 0.16, 0.42])
    ax_display = fig.add_axes([0.01, 0.26, 0.16, 0.16])
    ax_btn   = fig.add_axes([0.03, 0.10, 0.12, 0.12])
    ax_main  = fig.add_axes([0.22, 0.05, 0.76, 0.88])

    # --- Radio buttons for method selection ---
    radio = RadioButtons(ax_radio, METHOD_LABELS, active=0)
    ax_radio.set_title("Method", pad=6, fontsize=10)

    # --- Radio buttons for display mode ---
    display_radio = RadioButtons(ax_display, DISPLAY_LABELS, active=0)
    ax_display.set_title("Show", pad=6, fontsize=10)

    # --- Randomize button ---
    btn = Button(ax_btn, "Randomize\nStart / Goal", color="lightyellow", hovercolor="gold")

    # --- Mutable state ---
    state = {
        "start": random_state(),
        "goal":  random_state(),
        "path_type": METHOD_OPTIONS[METHOD_LABELS[0]],
        "display_mode": DISPLAY_OPTIONS[DISPLAY_LABELS[0]],
    }

    def refresh():
        compute_and_draw(
            ax_main,
            state["path_type"],
            state["start"],
            state["goal"],
            state["display_mode"],
        )

    def on_method_change(label):
        state["path_type"] = METHOD_OPTIONS[label]
        refresh()

    def on_display_change(label):
        state["display_mode"] = DISPLAY_OPTIONS[label]
        refresh()

    def on_randomize(_event):
        state["start"] = random_state()
        state["goal"]  = random_state()
        refresh()

    radio.on_clicked(on_method_change)
    display_radio.on_clicked(on_display_change)
    btn.on_clicked(on_randomize)

    # Initial draw
    refresh()
    plt.show()


if __name__ == "__main__":
    main()
