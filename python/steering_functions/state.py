# Copyright (c) 2017 - for information on the respective copyright owner
# see the NOTICE file and/or the repository
#     https://github.com/hbanzhaf/steering_functions.git
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""State and Control data structures for the steering functions library."""

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

    def __eq__(self, other: object) -> bool:
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

    def __repr__(self) -> str:
        return (
            f"State(x={self.x}, y={self.y}, theta={self.theta}, "
            f"kappa={self.kappa}, sigma={self.sigma}, d={self.d}, "
            f"s={self.s}, vel={self.vel}, acc={self.acc}, "
            f"time={self.time}, fork_y={self.fork_y})"
        )


@dataclass
class Control:
    """Description of a path segment with its corresponding control inputs."""

    delta_s: float = 0.0
    """Signed arc length of a segment."""

    kappa: float = 0.0
    """Curvature at the beginning of a segment."""

    sigma: float = 0.0
    """Sharpness (derivative of curvature w.r.t. arc length) of a segment."""

    def __repr__(self) -> str:
        return (
            f"Control(delta_s={self.delta_s}, kappa={self.kappa}, sigma={self.sigma})"
        )
