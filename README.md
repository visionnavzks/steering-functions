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

## Layout

- `include/steering_state/state.h`: shared public `State` type used by both libraries
- `src/steering_functions`: steering function implementations and unit tests
- `src/steering_path_lib`: wrapper library for selecting planners and computing paths
- `python`: pure-Python port of the public steering APIs in a separate package directory
