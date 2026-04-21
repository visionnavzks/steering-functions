use eframe::egui::{self, Color32, DragValue, RichText, ScrollArea};
use egui_plot::{Legend, Line, Plot, PlotBounds, PlotPoints, Points};

use steering_functions::{Control, PathType, State, SteeringPath};

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
    show_all: bool,
    kappa_max: f64,
    sigma_max: f64,
    discretization: f64,
    start: State,
    goal: State,
    drag_target: Option<Target>,
    fit_once: bool,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Target {
    Start,
    Goal,
}

impl Default for VisualizerApp {
    fn default() -> Self {
        Self {
            path_type: PathType::Dubins,
            show_all: false,
            kappa_max: 1.0,
            sigma_max: 1.0,
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
            drag_target: None,
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
            PathType::None => "None",
            PathType::CcDubins => "CC-Dubins",
            PathType::Cc00Dubins => "CC00-Dubins",
            PathType::Cc0pmDubins => "CC0pm-Dubins",
            PathType::Ccpm0Dubins => "CCpm0-Dubins",
            PathType::CcpmpmDubins => "CCpmpm-Dubins",
            PathType::Dubins => "Dubins",
            PathType::Cc00Rs => "CC00-RS",
            PathType::HcRs => "HC-RS",
            PathType::Hc00Rs => "HC00-RS",
            PathType::Hc0pmRs => "HC0pm-RS",
            PathType::Hcpm0Rs => "HCpm0-RS",
            PathType::HcpmpmRs => "HCpmpm-RS",
            PathType::Rs => "RS",
        }
    }

    fn planner(&self) -> Result<SteeringPath, String> {
        SteeringPath::try_new(
            self.path_type,
            self.kappa_max,
            self.sigma_max,
            self.discretization,
        )
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

    fn position_distance_squared(state: &State, x: f64, y: f64) -> f64 {
        let dx = state.x - x;
        let dy = state.y - y;
        dx * dx + dy * dy
    }

    fn pick_drag_target(&self, x: f64, y: f64) -> Option<Target> {
        let drag_radius_squared = 0.35_f64.powi(2);
        let start_distance = Self::position_distance_squared(&self.start, x, y);
        let goal_distance = Self::position_distance_squared(&self.goal, x, y);

        if start_distance <= drag_radius_squared && start_distance <= goal_distance {
            Some(Target::Start)
        } else if goal_distance <= drag_radius_squared {
            Some(Target::Goal)
        } else {
            None
        }
    }

    fn target_state_mut(&mut self, target: Target) -> &mut State {
        match target {
            Target::Start => &mut self.start,
            Target::Goal => &mut self.goal,
        }
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

    fn compute_view_data(&self) -> Result<PlanningViewData, String> {
        let planner = self.planner()?;

        if self.show_all {
            let all_paths = planner.compute_all_paths(&self.start, &self.goal)?;
            let all_controls = planner.compute_all_control_sequences(&self.start, &self.goal)?;
            let shortest_controls = planner.compute_shortest_control_sequence(&self.start, &self.goal)?;
            let summary = format!(
                "{} candidate paths, best length {:.3} m",
                all_paths.len(),
                shortest_controls.iter().map(|c| c.delta_s.abs()).sum::<f64>()
            );

            Ok(PlanningViewData {
                summary,
                shortest_controls,
                all_paths,
                all_controls,
                ..PlanningViewData::default()
            })
        } else {
            let shortest_path = planner.compute_shortest_path(&self.start, &self.goal)?;
            let shortest_controls = planner.compute_shortest_control_sequence(&self.start, &self.goal)?;
            let summary = format!(
                "{} segments, length {:.3} m",
                shortest_controls.len(),
                shortest_controls.iter().map(|c| c.delta_s.abs()).sum::<f64>()
            );

            Ok(PlanningViewData {
                summary,
                shortest_path: Some(shortest_path),
                shortest_controls,
                ..PlanningViewData::default()
            })
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
            .min_width(280.0)
            .show(ctx, |ui| {
                ui.heading("Steering Functions");
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

                ui.checkbox(&mut self.show_all, "Show all candidate paths");
                ui.add(egui::Slider::new(&mut self.kappa_max, 0.1..=5.0).text("kappa_max"));
                ui.add(egui::Slider::new(&mut self.sigma_max, 0.1..=5.0).text("sigma_max"));
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

                match &planning_data {
                    Ok(data) => {
                        ScrollArea::vertical().max_height(260.0).show(ui, |ui| {
                            if self.show_all {
                                for (sequence_index, controls) in data.all_controls.iter().enumerate() {
                                    let total_length: f64 = controls.iter().map(|c| c.delta_s.abs()).sum();
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
                            } else if data.shortest_controls.is_empty() {
                                ui.monospace("<empty>");
                            } else {
                                for (control_index, control) in data.shortest_controls.iter().enumerate() {
                                    ui.monospace(Self::format_control(control, control_index));
                                }
                            }
                        });
                    }
                    Err(error) => {
                        ui.colored_label(Color32::LIGHT_RED, error);
                    }
                }
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            let mut summary = String::new();

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

                    let pointer_down = ctx.input(|input| input.pointer.primary_down());
                    let pointer_position = plot_ui.pointer_coordinate();

                    if let Some(position) = pointer_position {
                        if follow_start_with_mouse {
                            self.start.x = position.x;
                            self.start.y = position.y;
                        }
                        if follow_goal_with_mouse {
                            self.goal.x = position.x;
                            self.goal.y = position.y;
                        }
                    }

                    if pointer_down {
                        if self.drag_target.is_none() {
                            if let Some(position) = pointer_position {
                                self.drag_target = self.pick_drag_target(position.x, position.y);
                            }
                        }

                        if let (Some(drag_target), Some(position)) = (self.drag_target, pointer_position) {
                            let target = self.target_state_mut(drag_target);
                            target.x = position.x;
                            target.y = position.y;
                        }
                    } else {
                        self.drag_target = None;
                    }

                    match &planning_data {
                        Ok(data) => {
                            summary = data.summary.clone();
                            if self.show_all {
                                for (index, path) in data.all_paths.iter().enumerate() {
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
                            } else if let Some(path) = &data.shortest_path {
                                plot_ui.line(
                                    Line::new(Self::path_plot_points(path))
                                        .name("shortest path")
                                        .color(Color32::LIGHT_BLUE)
                                        .width(2.5),
                                );
                            }
                        }
                        Err(error) => summary = error.clone(),
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
            ui.label(summary);
            ui.label("Hold S or G to let Start/Goal follow the mouse, or drag the green/red markers directly.");
        });
    }
}

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions::default();
    eframe::run_native(
        "steering-functions Rust visualizer",
        options,
        Box::new(|_cc| Ok(Box::new(VisualizerApp::default()))),
    )
}
