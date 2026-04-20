# steering-functions Rust

Rust port of the steering-functions planners, with optional Python bindings built through PyO3 and maturin.

## Native Rust

Compile the crate:

```bash
cargo check
```

## Python Bindings

From the `rust/` directory, install the extension into the active environment:

```bash
uv run --with maturin maturin develop
```

Then import it from Python:

```python
from steering_functions_rust import PathType, State, SteeringPath

planner = SteeringPath(PathType.DUBINS, 1.0, 1.0, 0.05)
start = State(x=0.0, y=0.0, theta=0.0)
goal = State(x=3.0, y=2.0, theta=1.0)
path = planner.compute_shortest_path(start, goal)
```

## Visualizer

Install the optional plotting dependency and launch the Python visualizer:

```bash
uv run --with maturin --with matplotlib maturin develop
uv run --with matplotlib python -m steering_functions_rust.visualize
```

## Supported Path Types

The current Rust port exposes these planners end-to-end:

- `DUBINS`
- `CC_DUBINS`
- `CC00_DUBINS`
- `CC0PM_DUBINS`
- `CCPM0_DUBINS`
- `CCPMPM_DUBINS`
- `RS`
- `CC00_RS`
- `HC00_RS`
- `HC0PM_RS`

The `HC_RS`, `HCPM0_RS`, and `HCPMPM_RS` enum variants are present for API compatibility but are not implemented in the Rust port yet.