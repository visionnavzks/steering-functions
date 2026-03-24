# Python Port of steering-functions

Pure-Python implementation of the steering functions library.  No compiled
extensions or third-party dependencies are required — only the Python standard
library.

## Implemented planners

| Class | Description |
|---|---|
| `DubinsStateSpace` | Dubins curves (forward-only car) |
| `ReedsSheppStateSpace` | Reeds-Shepp curves (forwards + reverse car) |

## Installation

No installation is required.  Add the `python/` directory to your `PYTHONPATH`
and import directly:

```python
import sys
sys.path.insert(0, "path/to/steering-functions/python")

from steering_functions import State, DubinsStateSpace, ReedsSheppStateSpace
```

## Quick start

```python
from steering_functions import State, DubinsStateSpace, ReedsSheppStateSpace

# --- Dubins (forward-only) ---
ds = DubinsStateSpace(kappa=1.0)   # kappa = max curvature [1/m]
s1 = State(x=0, y=0, theta=0)
s2 = State(x=5, y=2, theta=1.0)

print(ds.get_distance(s1, s2))     # path length
controls = ds.get_controls(s1, s2) # list of Control objects
path = ds.get_path(s1, s2)         # list of State objects (discretised)

# --- Reeds-Shepp (forward + reverse) ---
rs = ReedsSheppStateSpace(kappa=1.0)
path = rs.get_path(s1, s2)
```

## Running the tests

From the repository root:

```bash
python3 -m unittest discover -s python/tests
```

Or run a single test file:

```bash
python3 python/tests/test_dubins.py
python3 python/tests/test_reeds_shepp.py
```

## API

### `State`

Dataclass representing the full kinematic state of the robot.

| Field | Description |
|---|---|
| `x`, `y` | Position (m) |
| `theta` | Heading angle (rad) |
| `kappa` | Curvature (1/m) |
| `sigma` | Sharpness (1/m²) |
| `d` | Driving direction {−1, +1} |
| `s` | Arc length along path (m) |

### `Control`

Dataclass describing one path segment.

| Field | Description |
|---|---|
| `delta_s` | Signed arc length (m) |
| `kappa` | Curvature at segment start (1/m) |
| `sigma` | Sharpness (1/m²) |

### `DubinsStateSpace(kappa, discretization=0.1, forwards=True)`

| Method | Returns |
|---|---|
| `get_distance(s1, s2)` | `float` — shortest path length |
| `get_controls(s1, s2)` | `List[Control]` |
| `get_path(s1, s2)` | `List[State]` (discretised) |
| `get_all_controls(s1, s2)` | `List[List[Control]]` (forward + reverse) |
| `integrate(state, controls)` | `List[State]` |
| `interpolate(state, controls, t)` | `State` at normalised arc-length *t* ∈ [0, 1] |

### `ReedsSheppStateSpace(kappa, discretization=0.1)`

Same methods as `DubinsStateSpace` plus:

| Method | Returns |
|---|---|
| `reeds_shepp(s1, s2)` | `ReedsSheppPath` |
| `get_all_rs_paths(s1, s2)` | `List[ReedsSheppPath]` |
| `extract_controls_from_path(path)` | `List[Control]` |
