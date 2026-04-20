from __future__ import annotations

import math
from typing import Optional

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons, Slider

from steering_functions_rust import PathType, State, SteeringPath, supported_path_types


METHODS = {
    "Dubins": (PathType.DUBINS, False),
    "CC-Dubins": (PathType.CC_DUBINS, True),
    "CC00-Dubins": (PathType.CC00_DUBINS, True),
    "CC0pm-Dubins": (PathType.CC0PM_DUBINS, True),
    "CCpm0-Dubins": (PathType.CCPM0_DUBINS, True),
    "CCpmpm-Dubins": (PathType.CCPMPM_DUBINS, True),
    "RS": (PathType.RS, False),
    "CC00-RS": (PathType.CC00_RS, True),
    "HC-RS": (PathType.HC_RS, True),
    "HC00-RS": (PathType.HC00_RS, True),
    "HC0pm-RS": (PathType.HC0PM_RS, True),
    "HCpm0-RS": (PathType.HCPM0_RS, True),
    "HCpmpm-RS": (PathType.HCPMPM_RS, True),
}

SUPPORTED = list(supported_path_types())
METHODS = {
    name: config
    for name, config in METHODS.items()
    if any(config[0] == supported for supported in SUPPORTED)
}
METHOD_NAMES = list(METHODS)


def _path_xy(path):
    return [state.x for state in path], [state.y for state in path]


def _total_length(controls):
    return sum(abs(control.delta_s) for control in controls)


def _draw_pose(ax, state: State, color: str, label: str):
    ax.plot(state.x, state.y, "o", color=color, ms=8, label=label)
    dx = 0.35 * math.cos(state.theta)
    dy = 0.35 * math.sin(state.theta)
    ax.arrow(
        state.x,
        state.y,
        dx,
        dy,
        color=color,
        width=0.02,
        length_includes_head=True,
        head_width=0.14,
        head_length=0.18,
        zorder=5,
    )


class SteeringFunctionsRustVisualizer:
    def __init__(self):
        self.start = State(x=-3.0, y=0.0, theta=0.0)
        self.goal = State(x=3.0, y=2.0, theta=math.pi / 4)
        self.kappa_max = 1.0
        self.sigma_max = 1.0
        self.discretization = 0.02
        self.method_name = METHOD_NAMES[0]
        self.show_all = False
        self.active_target = "start"
        self._press_data: Optional[tuple[float, float]] = None

        self.fig = plt.figure(figsize=(14, 8))
        self.ax = self.fig.add_axes([0.05, 0.12, 0.58, 0.82])
        self.ax_method = self.fig.add_axes([0.68, 0.46, 0.28, 0.46], facecolor="#f2f2f2")
        self.ax_mode = self.fig.add_axes([0.68, 0.35, 0.28, 0.08], facecolor="#f2f2f2")
        self.ax_kappa = self.fig.add_axes([0.68, 0.24, 0.28, 0.03])
        self.ax_sigma = self.fig.add_axes([0.68, 0.18, 0.28, 0.03])
        self.ax_start = self.fig.add_axes([0.68, 0.10, 0.13, 0.05])
        self.ax_goal = self.fig.add_axes([0.83, 0.10, 0.13, 0.05])
        self.ax_reset = self.fig.add_axes([0.68, 0.03, 0.28, 0.05])

        self.ax.set_title("Rust steering-functions visualizer")
        self.ax.set_xlabel("x [m]")
        self.ax.set_ylabel("y [m]")
        self.ax.set_aspect("equal")
        self.ax.grid(True, ls="--", alpha=0.35)

        self.method_radio = RadioButtons(self.ax_method, METHOD_NAMES, active=0)
        self.mode_radio = RadioButtons(self.ax_mode, ["Best", "All"], active=0)
        self.kappa_slider = Slider(self.ax_kappa, "kappa_max", 0.1, 5.0, valinit=1.0, valstep=0.05)
        self.sigma_slider = Slider(self.ax_sigma, "sigma_max", 0.1, 5.0, valinit=1.0, valstep=0.05)
        self.start_button = Button(self.ax_start, "Set Start", color="#d8f0d8")
        self.goal_button = Button(self.ax_goal, "Set Goal", color="#f4d6d6")
        self.reset_button = Button(self.ax_reset, "Reset")

        self.info_text = self.ax.text(
            0.01,
            0.01,
            "",
            transform=self.ax.transAxes,
            va="bottom",
            fontsize=9,
            bbox={"boxstyle": "round", "facecolor": "#fff6bf", "alpha": 0.9},
        )

        self.method_radio.on_clicked(self._on_method_change)
        self.mode_radio.on_clicked(self._on_mode_change)
        self.kappa_slider.on_changed(self._on_kappa_change)
        self.sigma_slider.on_changed(self._on_sigma_change)
        self.start_button.on_clicked(lambda _event: self._set_target("start"))
        self.goal_button.on_clicked(lambda _event: self._set_target("goal"))
        self.reset_button.on_clicked(self._on_reset)

        self.fig.canvas.mpl_connect("button_press_event", self._on_press)
        self.fig.canvas.mpl_connect("button_release_event", self._on_release)

        self._draw()

    def _planner(self) -> SteeringPath:
        path_type, needs_sigma = METHODS[self.method_name]
        sigma_max = self.sigma_max if needs_sigma else 1.0
        return SteeringPath(path_type, self.kappa_max, sigma_max, self.discretization)

    def _set_target(self, target: str):
        self.active_target = target

    def _on_method_change(self, label: str):
        self.method_name = label
        self._draw()

    def _on_mode_change(self, label: str):
        self.show_all = label == "All"
        self._draw()

    def _on_kappa_change(self, value: float):
        self.kappa_max = value
        self._draw()

    def _on_sigma_change(self, value: float):
        self.sigma_max = value
        self._draw()

    def _on_reset(self, _event):
        self.start = State(x=-3.0, y=0.0, theta=0.0)
        self.goal = State(x=3.0, y=2.0, theta=math.pi / 4)
        self.kappa_slider.set_val(1.0)
        self.sigma_slider.set_val(1.0)
        self._draw()

    def _on_press(self, event):
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return
        self._press_data = (event.xdata, event.ydata)

    def _on_release(self, event):
        if self._press_data is None or event.inaxes != self.ax:
            self._press_data = None
            return

        x0, y0 = self._press_data
        x1 = event.xdata if event.xdata is not None else x0
        y1 = event.ydata if event.ydata is not None else y0
        theta = math.atan2(y1 - y0, x1 - x0) if abs(x1 - x0) + abs(y1 - y0) > 1e-9 else 0.0

        target = self.start if self.active_target == "start" else self.goal
        target.x = x0
        target.y = y0
        target.theta = theta
        self._press_data = None
        self._draw()

    def _draw(self):
        self.ax.clear()
        self.ax.set_title("Rust steering-functions visualizer")
        self.ax.set_xlabel("x [m]")
        self.ax.set_ylabel("y [m]")
        self.ax.set_aspect("equal")
        self.ax.grid(True, ls="--", alpha=0.35)

        planner = self._planner()
        try:
            if self.show_all:
                paths = planner.compute_all_paths(self.start, self.goal)
                for index, path in enumerate(paths):
                    xs, ys = _path_xy(path)
                    self.ax.plot(xs, ys, lw=2.0, alpha=0.75, label=f"path {index + 1}")
                controls = planner.compute_shortest_control_sequence(self.start, self.goal)
                info = f"{self.method_name} | {len(paths)} paths | best length={_total_length(controls):.3f} m"
            else:
                path = planner.compute_shortest_path(self.start, self.goal)
                controls = planner.compute_shortest_control_sequence(self.start, self.goal)
                xs, ys = _path_xy(path)
                self.ax.plot(xs, ys, color="#1f77b4", lw=2.5, label="shortest path")
                info = f"{self.method_name} | segments={len(controls)} | length={_total_length(controls):.3f} m"
        except NotImplementedError as exc:
            info = str(exc)
        except Exception as exc:
            info = f"Planner error: {exc}"

        _draw_pose(self.ax, self.start, "#2ca02c", "start")
        _draw_pose(self.ax, self.goal, "#d62728", "goal")
        self.info_text = self.ax.text(
            0.01,
            0.01,
            info,
            transform=self.ax.transAxes,
            va="bottom",
            fontsize=9,
            bbox={"boxstyle": "round", "facecolor": "#fff6bf", "alpha": 0.9},
        )
        self.ax.legend(loc="upper left")
        self.fig.canvas.draw_idle()


def main():
    SteeringFunctionsRustVisualizer()
    plt.show()


if __name__ == "__main__":
    main()