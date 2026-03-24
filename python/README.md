# Python steering-functions

This directory provides a Python rewrite of the package interface, centered on
the `SteeringPath` API and associated `State`, `Control`, and `PathType` types.

## Run tests

From repository root:

```bash
python3 -m unittest discover -s python/tests
```

## Basic usage

```python
from steering_functions import PathType, State, SteeringPath

planner = SteeringPath(PathType.DUBINS, 1.0, 1.0, 0.1)
path = planner.compute_shortest_path(State(0, 0, 0), State(5, 0, 0))
```
