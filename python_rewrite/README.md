# steering-functions (Python)

A pure-Python port of the C++ **steering_functions** library for computing
optimal and near-optimal paths for car-like robots.

## Supported steering functions

| Type | Continuity | Forwards only | Notes |
|------|-----------|---------------|-------|
| Dubins | G¹ | yes | optimal |
| Reeds-Shepp | G¹ | no | optimal |
| CC-Dubins (CC00, CC0pm, CCpm0, CCpmpm) | G² | yes | suboptimal, continuous curvature |
| HC-Reeds-Shepp (HC00, HC0pm, HCpm0, HCpmpm) | G² between cusps | no | suboptimal |
| CC-Reeds-Shepp (CC00) | G² | no | suboptimal |

## Quick start

```python
from steering_functions import SteeringPath, PathType, State

planner = SteeringPath(PathType.DUBINS, kappa_max=1.0, sigma_max=0.0, discretization=0.1)

start = State(x=0, y=0, theta=0)
goal  = State(x=5, y=2, theta=1.57)

path = planner.compute_shortest_path(start, goal)
for s in path:
    print(f"  x={s.x:.3f}  y={s.y:.3f}  θ={s.theta:.3f}")
```

## Running tests

```bash
cd python
python -m pytest tests/ -v
```
