import math
from dataclasses import dataclass, field


@dataclass
class State:
    """Description of a kinematic car's state."""

    x: float = 0.0
    y: float = 0.0
    theta: float = 0.0
    kappa: float = 0.0
    sigma: float = 0.0
    d: float = 0.0
    s: float = 0.0
    vel: float = 0.0
    acc: float = 0.0
    time: float = 0.0
    fork_y: float = 0.0

    def __eq__(self, other):
        if not isinstance(other, State):
            return NotImplemented
        eps = 1e-6
        return (
            abs(self.x - other.x) < eps
            and abs(self.y - other.y) < eps
            and abs(self.theta - other.theta) < eps
            and abs(self.kappa - other.kappa) < eps
            and abs(self.sigma - other.sigma) < eps
            and abs(self.d - other.d) < eps
            and abs(self.s - other.s) < eps
            and abs(self.vel - other.vel) < eps
            and abs(self.acc - other.acc) < eps
            and abs(self.time - other.time) < eps
            and abs(self.fork_y - other.fork_y) < eps
        )



@dataclass
class Control:
    """Description of a path segment with its corresponding control inputs."""

    delta_s: float = 0.0
    kappa: float = 0.0
    sigma: float = 0.0
