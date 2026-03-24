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

Useful CMake options:

- `-DSTEERING_FUNCTIONS_ENABLE_EXTENDED_TESTS=ON` to include the larger upstream integration test executable in `ctest`

## Python

The repository now also exposes the complete high-level planning API as a Python package backed by the existing implementation.

Install it from the repository root:

```bash
python3 -m pip install .
```

Run the focused Python binding tests:

```bash
python3 -m unittest discover -s python/tests -v
```

Example usage:

```python
from steering_functions import PathType, State, SteeringPath

planner = SteeringPath(PathType.RS, 1.0, 1.0, 0.1)
start = State(0.0, 0.0, 0.0)
goal = State(4.0, 2.0, 1.0)
path = planner.computeShortestPath(start, goal)
```

## Layout

- `include/steering_state/state.h`: shared public `State` type used by both libraries
- `src/steering_functions`: steering function implementations and unit tests
- `src/steering_path_lib`: wrapper library for selecting planners and computing paths
