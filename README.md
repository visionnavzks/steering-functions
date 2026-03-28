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

## Steering Function Comparison

All steering functions target car-like robots with a bounded turning radius. They differ in **driving direction**, **path smoothness**, and whether curvature at the endpoints is constrained.

### Summary Table

| Steering Function | Driving Direction | Smoothness | Path Length | Comp. Time (mean) |
|---|---|---|---|---|
| Dubins | forwards only | G¹ | optimal | ~1.3 µs |
| CC±±-Dubins | forwards only | G² | suboptimal | ~7.1 µs |
| CC±0-Dubins | forwards only | G² | suboptimal | ~5.6 µs |
| CC0±-Dubins | forwards only | G² | suboptimal | ~5.6 µs |
| CC00-Dubins | forwards only | G² | suboptimal | ~4.9 µs |
| CC-Dubins | forwards only | G² | suboptimal | ~14.8 µs |
| Reeds-Shepp | forwards + backwards | G¹ | optimal | ~7.4 µs |
| HC±±-Reeds-Shepp | forwards + backwards | G² between cusps | suboptimal | ~58 µs |
| HC±0-Reeds-Shepp | forwards + backwards | G² between cusps | suboptimal | ~56 µs |
| HC0±-Reeds-Shepp | forwards + backwards | G² between cusps | suboptimal | ~51 µs |
| HC00-Reeds-Shepp | forwards + backwards | G² between cusps | suboptimal | ~56 µs |
| HC-Reeds-Shepp | forwards + backwards | G² between cusps | suboptimal | ~464 µs |
| CC00-Reeds-Shepp | forwards + backwards | G² | suboptimal | ~54 µs |

Computation times are measured on a single core of an Intel Xeon E5 @ 3.50 GHz over 10⁵ random steering queries.

### Key Dimensions

**Driving direction**
- *Dubins* / *CC-Dubins*: forwards only.
- *Reeds-Shepp* / *HC-Reeds-Shepp* / *CC-Reeds-Shepp*: forwards and backwards (path may include cusps / direction reversals).

**Path smoothness**
- G¹ (Dubins, Reeds-Shepp): orientation is continuous but curvature may jump instantaneously. Suitable for kinematic planning where steering-angle discontinuities are acceptable.
- G² between cusps (HC-Reeds-Shepp family): curvature is continuous within each driving segment; a curvature discontinuity is only permitted at a cusp (direction reversal). Achieved by inserting clothoid transition arcs.
- G² everywhere (CC-Dubins, CC-Reeds-Shepp): curvature is continuous along the full path including at cusps. Strictest smoothness guarantee; useful when actuator bandwidth limits how fast the steering angle can change.

**Endpoint curvature constraint** (HC-Reeds-Shepp and CC-Dubins variants)

The superscripts in the name encode what curvature is enforced at the start (first superscript) and goal (second superscript):

| Superscript | Meaning |
|---|---|
| `0` | Zero curvature is enforced regardless of the input value |
| `±` | The planner picks positive or negative maximum curvature; if the caller provides a non-zero value, the corresponding signed maximum curvature is used |
| *(none / bare name)* | Accepts an arbitrary user-specified curvature at both ends |

Choosing a constrained variant (e.g. HC00-RS) instead of the fully general one (HC-RS) reduces computation time by roughly 8×, making it preferable in tight motion-planning loops when the boundary conditions are known in advance.

**Optimality**
- Dubins and Reeds-Shepp are provably length-optimal within their respective model classes.
- All G² variants are suboptimal: the clothoid arcs added for curvature continuity increase path length slightly compared to the G¹ baseline. The histograms in `src/steering_functions/doc/images/statistics.png` show the empirical length overhead and the reduction in curvature discontinuities.

### How to Choose

1. **Need backwards driving?** → Reeds-Shepp family; otherwise Dubins family.
2. **Need smooth curvature (G²)?** → CC or HC variants; otherwise standard Dubins / Reeds-Shepp.
3. **Need G² at cusps too?** → CC-Reeds-Shepp; HC-Reeds-Shepp only guarantees G² *between* cusps.
4. **Boundary curvature known at planning time?** → Pick the matching HC superscript variant to save ~8× compute versus HC-RS.
5. **Arbitrary boundary curvature at runtime?** → HC-RS or CC-Dubins (accepts any start/goal curvature).