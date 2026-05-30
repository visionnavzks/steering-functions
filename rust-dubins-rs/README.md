# steering-functions-dubins-rs

A minimal Rust crate for steering functions including:

- **Dubins** curves (forward-only, reverse-only, or bidirectional)
- **Reeds-Shepp** curves (bidirectional)

## Build

```bash
cargo check
```

## Test

```bash
cargo test
```

## GUI Example

Run the interactive planner visualizer:

```bash
cargo run --example gui
```

The GUI supports:

- Dubins and Reeds-Shepp path types
- Dubins direction mode: forward only, reverse only, forward or reverse
- Shortest path and all candidate paths display
- Control commands with segment type labels

## Usage

```rust
use steering_functions_dubins_rs::{PathType, State, SteeringPath};

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

## License

Apache-2.0
