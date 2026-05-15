use eframe::egui::{self, Color32, DragValue, RichText, ScrollArea};
use egui_plot::{Legend, Line, Plot, PlotBounds, PlotPoints, Points};

use steering_functions_lite::{
    Control, DubinsDirectionMode, PathType, State, SteeringPath,
};

#[derive(Default)]
struct PlanningViewData {
    summary: String,
    shortest_path: Option<Vec<State>>,
    shortest_controls: Vec<Control>,
    all_paths: Vec<Vec<State>>,
    all_controls: Vec<Vec<Control>>,
}

struct VisualizerApp {
    path_type: PathType,
    dubins_mode: DubinsDirectionMode,
    show_all: bool,
    kappa_max: f64,
    discretization: f64,
    start: State,
    goal: State,
    fit_once: bool,
}

impl Default for VisualizerApp {
    fn default() -> Self {
        Self {
            path_type: PathType::Dubins,
            dubins_mode: DubinsDirectionMode::ForwardOnly,
            show_all: false,
            kappa_max: 1.0,
            discretization: 0.05,
            start: State {
                x: -3.0,
                y: 0.0,
                theta: 0.0,
                kappa: 0.0,
                ..State::default()
            },
            goal: State {
                x: 3.0,
                y: 2.0,
                theta: std::f64::consts::FRAC_PI_4,
                kappa: 0.0,
                ..State::default()
            },
            fit_once: true,
        }
    }
}

impl VisualizerApp {
    fn reset(&mut self) {
        *self = Self::default();
    }

    fn path_type_label(path_type: PathType) -> &'static str {
        match path_type {
            PathType::Dubins => "Dubins",
            PathType::Rs => "Reeds-Shepp",
        }
    }

    fn dubins_mode_label(mode: DubinsDirectionMode) -> &'static str {
        match mode {
            DubinsDirectionMode::ForwardOnly => "Forward only",
            DubinsDirectionMode::ReverseOnly => "Reverse only",
            DubinsDirectionMode::ForwardOrReverse => "Forward or reverse",
        }
    }

    fn planner(&self) -> SteeringPath {
        SteeringPath::new(self.path_type, self.kappa_max, self.discretization)
            .with_dubins_direction_mode(self.dubins_mode)
    }

    fn arrow_points(state: &State, length: f64) -> Vec<[f64; 2]> {
        vec![
            [state.x, state.y],
            [state.x + length * state.theta.cos(), state.y + length * state.theta.sin()],
        ]
    }

    fn path_plot_points(path: &[State]) -> PlotPoints {
        PlotPoints::from_iter(path.iter().map(|state| [state.x, state.y]))
    }

    fn edit_state(ui: &mut egui::Ui, label: &str, state: &mut State, kappa_max: f64) {
        ui.label(RichText::new(label).strong());
        ui.horizontal(|ui| {
            ui.label("x");
            ui.add(DragValue::new(&mut state.x).speed(0.05));
            ui.label("y");
            ui.add(DragValue::new(&mut state.y).speed(0.05));
        });
        ui.horizontal(|ui| {
            ui.label("theta");
            ui.add(DragValue::new(&mut state.theta).speed(0.02));
            ui.label("kappa");
            ui.add(DragValue::new(&mut state.kappa).speed(0.02).range(-kappa_max..=kappa_max));
        });
    }

    fn format_control(control: &Control, index: usize) -> String {
        format!(
            "{:02}: ds={:+.3}, kappa={:+.3}, sigma={:+.3}",
            index,
            control.delta_s,
            control.kappa,
            control.sigma
        )
    }

    fn path_length(controls: &[Control]) -> f64 {
        controls.iter().map(|c| c.delta_s.abs()).sum()
    }

    fn compute_view_data(&self) -> PlanningViewData {
        let planner = self.planner();

        if self.show_all {
            let all_paths = planner.compute_all_paths(&self.start, &self.goal);
            let all_controls = planner.compute_all_control_sequences(&self.start, &self.goal);
            let shortest_controls = planner.compute_shortest_control_sequence(&self.start, &self.goal);
            let summary = format!(
                "{} candidate paths, best length {:.3} m",
                all_paths.len(),
                Self::path_length(&shortest_controls)
            );

            PlanningViewData {
                summary,
                shortest_controls,
                all_paths,
                all_controls,
                ..PlanningViewData::default()
            }
        } else {
            let shortest_path = planner.compute_shortest_path(&self.start, &self.goal);
            let shortest_controls = planner.compute_shortest_control_sequence(&self.start, &self.goal);
            let summary = format!(
                "{} segments, length {:.3} m",
                shortest_controls.len(),
                Self::path_length(&shortest_controls)
            );

            PlanningViewData {
                summary,
                shortest_path: Some(shortest_path),
                shortest_controls,
                ..PlanningViewData::default()
            }
        }
    }
}

impl eframe::App for VisualizerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let mut follow_start_with_mouse = false;
        let mut follow_goal_with_mouse = false;

        ctx.input(|input| {
            follow_start_with_mouse = input.key_down(egui::Key::S);
            follow_goal_with_mouse = input.key_down(egui::Key::G);
        });

        let planning_data = self.compute_view_data();

        egui::SidePanel::right("controls")
            .min_width(300.0)
            .show(ctx, |ui| {
                ui.heading("Steering Functions Lite");
                ui.separator();

                egui::ComboBox::from_label("Method")
                    .selected_text(Self::path_type_label(self.path_type))
                    .show_ui(ui, |ui| {
                        for path_type in SteeringPath::supported_path_types() {
                            ui.selectable_value(
                                &mut self.path_type,
                                path_type,
                                Self::path_type_label(path_type),
                            );
                        }
                    });

                if self.path_type == PathType::Dubins {
                    egui::ComboBox::from_label("Dubins mode")
                        .selected_text(Self::dubins_mode_label(self.dubins_mode))
                        .show_ui(ui, |ui| {
                            for mode in [
                                DubinsDirectionMode::ForwardOnly,
                                DubinsDirectionMode::ReverseOnly,
                                DubinsDirectionMode::ForwardOrReverse,
                            ] {
                                ui.selectable_value(
                                    &mut self.dubins_mode,
                                    mode,
                                    Self::dubins_mode_label(mode),
                                );
                            }
                        });
                }

                ui.checkbox(&mut self.show_all, "Show all candidate paths");
                ui.add(egui::Slider::new(&mut self.kappa_max, 0.1..=5.0).text("kappa_max"));
                ui.add(
                    egui::Slider::new(&mut self.discretization, 0.01..=0.2)
                        .text("discretization"),
                );

                ui.separator();
                Self::edit_state(ui, "Start", &mut self.start, self.kappa_max);
                ui.separator();
                Self::edit_state(ui, "Goal", &mut self.goal, self.kappa_max);
                ui.separator();
                ui.label("Hotkeys: hold S to move Start with the mouse, hold G to move Goal.");

                if ui.button("Reset").clicked() {
                    self.reset();
                }

                ui.separator();
                ui.label(RichText::new("Control Commands").strong());
                ScrollArea::vertical().max_height(280.0).show(ui, |ui| {
                    if self.show_all {
                        for (sequence_index, controls) in planning_data.all_controls.iter().enumerate() {
                            let total_length = Self::path_length(controls);
                            ui.label(
                                RichText::new(format!(
                                    "Path {} ({:.3} m)",
                                    sequence_index + 1,
                                    total_length
                                ))
                                .strong(),
                            );
                            if controls.is_empty() {
                                ui.monospace("  <empty>");
                            } else {
                                for (control_index, control) in controls.iter().enumerate() {
                                    ui.monospace(Self::format_control(control, control_index));
                                }
                            }
                            ui.add_space(6.0);
                        }
                    } else if planning_data.shortest_controls.is_empty() {
                        ui.monospace("<empty>");
                    } else {
                        for (control_index, control) in planning_data.shortest_controls.iter().enumerate() {
                            ui.monospace(Self::format_control(control, control_index));
                        }
                    }
                });
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            Plot::new("steering_plot")
                .legend(Legend::default())
                .data_aspect(1.0)
                .allow_scroll(false)
                .allow_zoom(false)
                .show(ui, |plot_ui| {
                    if self.fit_once {
                        plot_ui.set_plot_bounds(PlotBounds::from_min_max([-6.0, -6.0], [6.0, 6.0]));
                        self.fit_once = false;
                    }

                    if plot_ui.response().hovered() {
                        let scroll_delta_y = ctx.input(|input| input.smooth_scroll_delta.y);
                        if scroll_delta_y.abs() > f32::EPSILON {
                            let zoom_factor = (scroll_delta_y / 200.0).exp();
                            plot_ui.set_auto_bounds(false.into());
                            plot_ui.zoom_bounds_around_hovered(egui::Vec2::splat(zoom_factor));
                        }
                    }

                    if let Some(position) = plot_ui.pointer_coordinate() {
                        if follow_start_with_mouse {
                            self.start.x = position.x;
                            self.start.y = position.y;
                        }
                        if follow_goal_with_mouse {
                            self.goal.x = position.x;
                            self.goal.y = position.y;
                        }
                    }

                    if self.show_all {
                        for (index, path) in planning_data.all_paths.iter().enumerate() {
                            let color = Color32::from_rgb(
                                (40 + (index * 45 % 180)) as u8,
                                (120 + (index * 35 % 100)) as u8,
                                (200 - (index * 25 % 120)) as u8,
                            );
                            plot_ui.line(
                                Line::new(Self::path_plot_points(path))
                                    .name(format!("path {}", index + 1))
                                    .color(color),
                            );
                        }
                    } else if let Some(path) = &planning_data.shortest_path {
                        plot_ui.line(
                            Line::new(Self::path_plot_points(path))
                                .name("shortest path")
                                .color(Color32::LIGHT_BLUE)
                                .width(2.5),
                        );
                    }

                    plot_ui.points(
                        Points::new(vec![[self.start.x, self.start.y]])
                            .name("start")
                            .color(Color32::GREEN)
                            .radius(6.0),
                    );
                    plot_ui.points(
                        Points::new(vec![[self.goal.x, self.goal.y]])
                            .name("goal")
                            .color(Color32::RED)
                            .radius(6.0),
                    );
                    plot_ui.line(
                        Line::new(Self::arrow_points(&self.start, 0.4))
                            .name("start heading")
                            .color(Color32::GREEN),
                    );
                    plot_ui.line(
                        Line::new(Self::arrow_points(&self.goal, 0.4))
                            .name("goal heading")
                            .color(Color32::RED),
                    );
                });

            ui.separator();
            ui.label(&planning_data.summary);
            ui.label("Hold S or G to let Start/Goal follow the mouse.");
        });
    }
}

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions::default();
    eframe::run_native(
        "steering-functions-lite GUI",
        options,
        Box::new(|_cc| Ok(Box::new(VisualizerApp::default()))),
    )
}
