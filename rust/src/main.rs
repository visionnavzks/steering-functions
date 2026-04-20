use eframe::egui::{self, Color32, DragValue, RichText};
use egui_plot::{Legend, Line, Plot, PlotBounds, PlotPoints, Points};

use steering_functions::{PathType, State, SteeringPath};

struct VisualizerApp {
    path_type: PathType,
    show_all: bool,
    kappa_max: f64,
    sigma_max: f64,
    discretization: f64,
    start: State,
    goal: State,
    set_target: Target,
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
            set_target: Target::Start,
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
}

impl eframe::App for VisualizerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
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

                ui.horizontal(|ui| {
                    ui.selectable_value(&mut self.set_target, Target::Start, "Click sets Start");
                    ui.selectable_value(&mut self.set_target, Target::Goal, "Click sets Goal");
                });

                if ui.button("Reset").clicked() {
                    self.reset();
                }
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            let planner = self.planner();
            let mut summary = String::new();

            Plot::new("steering_plot")
                .legend(Legend::default())
                .data_aspect(1.0)
                .show(ui, |plot_ui| {
                    if self.fit_once {
                        plot_ui.set_plot_bounds(PlotBounds::from_min_max([-6.0, -6.0], [6.0, 6.0]));
                        self.fit_once = false;
                    }

                    if plot_ui.response().clicked() {
                        if let Some(position) = plot_ui.pointer_coordinate() {
                            let target = match self.set_target {
                                Target::Start => &mut self.start,
                                Target::Goal => &mut self.goal,
                            };
                            target.x = position.x;
                            target.y = position.y;
                        }
                    }

                    match planner {
                        Ok(ref planner) => {
                            if self.show_all {
                                match planner.compute_all_paths(&self.start, &self.goal) {
                                    Ok(paths) => {
                                        for (index, path) in paths.iter().enumerate() {
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
                                        match planner.compute_shortest_control_sequence(&self.start, &self.goal) {
                                            Ok(controls) => {
                                                summary = format!(
                                                    "{} candidate paths, best length {:.3} m",
                                                    paths.len(),
                                                    controls.iter().map(|c| c.delta_s.abs()).sum::<f64>()
                                                );
                                            }
                                            Err(error) => summary = error,
                                        }
                                    }
                                    Err(error) => summary = error,
                                }
                            } else {
                                match planner.compute_shortest_path(&self.start, &self.goal) {
                                    Ok(path) => {
                                        plot_ui.line(
                                            Line::new(Self::path_plot_points(&path))
                                                .name("shortest path")
                                                .color(Color32::LIGHT_BLUE)
                                                .width(2.5),
                                        );
                                        match planner.compute_shortest_control_sequence(&self.start, &self.goal) {
                                            Ok(controls) => {
                                                summary = format!(
                                                    "{} segments, length {:.3} m",
                                                    controls.len(),
                                                    controls.iter().map(|c| c.delta_s.abs()).sum::<f64>()
                                                );
                                            }
                                            Err(error) => summary = error,
                                        }
                                    }
                                    Err(error) => summary = error,
                                }
                            }
                        }
                        Err(error) => summary = error,
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
            ui.label("Click inside the plot to move the selected target point.");
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
