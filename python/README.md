# Python steering-functions port

This directory contains a pure-Python rewrite of the repository's public steering API.

## Scope

The Python port mirrors the shared data types and the high-level `SteeringPath` facade in a
new package:

- `steering_functions.State`
- `steering_functions.Control`
- `steering_functions.DubinsStateSpace`
- `steering_functions.ReedsSheppStateSpace`
- `steering_functions.PathType`
- `steering_functions.SteeringPath`

The current Python rewrite supports the `DUBINS` and `RS` planners. Continuous-curvature and
hybrid-curvature planner types remain available in the enum for API compatibility and raise
`NotImplementedError` in the Python facade.

## Run the Python tests

```bash
cd /home/runner/work/steering-functions/steering-functions
python3 -m unittest discover -s python/tests
```
