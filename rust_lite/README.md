# steering-functions-lite

A minimal Rust extraction of steering functions that only includes:

- Dubins
- Reeds-Shepp

## Build

```bash
cargo check
```

## GUI Example

Run the interactive planner visualizer:

```bash
cargo run --example gui
```

The GUI supports:

- Dubins and Reeds-Shepp
- Dubins direction mode: forward only, reverse only, forward or reverse
- shortest path and all candidate paths display

## Usage

```rust
use steering_functions_lite::{PathType, State, SteeringPath};

let planner = SteeringPath::new(PathType::Dubins, 1.0, 0.05);
let start = State {
    x: 0.0,
    y: 0.0,
    theta: 0.0,
    ..State::default()
};
let goal = State {
    x: 3.0,
    y: 2.0,
    theta: 1.0,
    ..State::default()
};

let path = planner.compute_shortest_path(&start, &goal);
println!("path samples: {}", path.len());
```
