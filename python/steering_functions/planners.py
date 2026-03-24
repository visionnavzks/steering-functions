from __future__ import annotations

from typing import List

from .base import BaseStateSpace
from .types import Control, State


class _SimplifiedPlanner(BaseStateSpace):
    def get_controls(self, state1: State, state2: State) -> List[Control]:
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        distance = (dx**2 + dy**2) ** 0.5
        return [Control(delta_s=distance, kappa=0.0, sigma=0.0)]

    def get_all_controls(self, state1: State, state2: State) -> List[List[Control]]:
        return [self.get_controls(state1, state2)]


class DubinsStateSpace(_SimplifiedPlanner):
    pass


class CCDubinsStateSpace(_SimplifiedPlanner):
    pass


class CC00DubinsStateSpace(_SimplifiedPlanner):
    pass


class CC0PMDubinsStateSpace(_SimplifiedPlanner):
    pass


class CCpm0DubinsStateSpace(_SimplifiedPlanner):
    pass


class CCpmpmDubinsStateSpace(_SimplifiedPlanner):
    pass


class ReedsSheppStateSpace(_SimplifiedPlanner):
    pass


class CC00ReedsSheppStateSpace(_SimplifiedPlanner):
    pass


class HCReedsSheppStateSpace(_SimplifiedPlanner):
    pass


class HC00ReedsSheppStateSpace(_SimplifiedPlanner):
    pass


class HC0PMReedsSheppStateSpace(_SimplifiedPlanner):
    pass


class HCpm0ReedsSheppStateSpace(_SimplifiedPlanner):
    pass


class HCpmpmReedsSheppStateSpace(_SimplifiedPlanner):
    pass
