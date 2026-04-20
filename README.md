# steering-functions

Standalone workspace wrapper around the upstream steering functions library and a higher-level steering path library.

## Build

Build both libraries from the repository root:

```bash
./build.sh
```

Build and run the unit tests:

```bash
BUILD_TESTING=ON RUN_TESTS=1 ./build.sh
```

Useful environment variables:

- `BUILD_DIR` to change the build directory. Default: `build`
- `BUILD_TYPE` to set the CMake build type. Default: `Release`
- `BUILD_TESTING` to enable or disable tests. Default: `OFF`
- `RUN_TESTS` to run `ctest` after the build. Default: `0`
- `BUILD_PYTHON_BINDINGS` to include the pybind11 extension in the root CMake build. Default: `OFF`

Useful CMake options:

- `-DSTEERING_FUNCTIONS_ENABLE_EXTENDED_TESTS=ON` to include the larger upstream integration test executable in `ctest`

## Python

The repository now also exposes the complete high-level planning API as a Python package backed by the existing implementation.

There are now two supported build flows.

Build the extension through CMake and copy it into `python/steering_functions/`:

```bash
bash ./build_python.sh
```

Build, then run the Python unit tests:

```bash
export RUN_PYTHON_TESTS=1
bash ./build_python.sh
```

Build, then launch the interactive demo:

```bash
export RUN_PYTHON_DEMO=1
bash ./build_python.sh
```

You can also include the bindings in a normal root CMake build:

```bash
./build.sh -DBUILD_PYTHON_BINDINGS=ON
```

Install it from the repository root:

```bash
python3 -m pip install .
```

Run the focused Python binding tests:

```bash
python3 -m unittest discover -s python/tests -v
```

Run the interactive demo manually:

```bash
export PYTHONPATH="$PWD/python${PYTHONPATH:+:$PYTHONPATH}"
python3 python/demo_interactive.py
```

Example usage:

```python
from steering_functions import PathType, State, SteeringPath

planner = SteeringPath(PathType.RS, 1.0, 1.0, 0.1)
start = State(0.0, 0.0, 0.0)
goal = State(4.0, 2.0, 1.0)
path = planner.computeShortestPath(start, goal)
```

## Rust

The `rust/` crate contains the Rust port. It now also exposes Python bindings through PyO3/maturin.

Build-check the Rust crate:

```bash
cd rust
cargo check
```

Install the Rust Python bindings into the active environment:

```bash
cd rust
uv run --with maturin maturin develop
```

Launch the Rust-backed visualizer:

```bash
cd rust
uv run --with matplotlib python -m steering_functions_rust.visualize
```

Launch the native Rust egui visualizer:

```bash
cd rust
cargo run --bin visualizer
```

Current Rust binding coverage includes Dubins, CC-Dubins variants, RS, CC00-RS, HC-RS, HC00-RS, HC0pm-RS, HCpm0-RS, and HCpmpm-RS.

## Layout

- `include/steering_state/state.h`: shared public `State` type used by both libraries
- `src/steering_functions`: steering function implementations and unit tests
- `src/steering_path_lib`: wrapper library for selecting planners and computing paths
- `src/python_bindings`: pybind11 module source and CMake integration
- `python`: Python package, demo, and Python tests
