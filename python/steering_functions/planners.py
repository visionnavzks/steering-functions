from __future__ import annotations

from typing import List

from .base import BaseStateSpace
from .types import Control, State


class _StraightLinePlanner(BaseStateSpace):
    def get_controls(self, state1: State, state2: State) -> List[Control]:
        dx = state2.x - state1.x
        dy = state2.y - state1.y
        distance = (dx**2 + dy**2) ** 0.5
        return [Control(delta_s=distance, kappa=0.0, sigma=0.0)]

    def get_all_controls(self, state1: State, state2: State) -> List[List[Control]]:
        return [self.get_controls(state1, state2)]


class DubinsStateSpace(_StraightLinePlanner):
    pass


class CCDubinsStateSpace(_StraightLinePlanner):
    pass


class CC00DubinsStateSpace(_StraightLinePlanner):
    pass


class CC0pmDubinsStateSpace(_StraightLinePlanner):
    pass


class CCpm0DubinsStateSpace(_StraightLinePlanner):
    pass


class CCpmpmDubinsStateSpace(_StraightLinePlanner):
    pass


class ReedsSheppStateSpace(_StraightLinePlanner):
    pass


class CC00ReedsSheppStateSpace(_StraightLinePlanner):
    pass


class HCReedsSheppStateSpace(_StraightLinePlanner):
    pass


class HC00ReedsSheppStateSpace(_StraightLinePlanner):
    pass


class HC0pmReedsSheppStateSpace(_StraightLinePlanner):
    pass


class HCpm0ReedsSheppStateSpace(_StraightLinePlanner):
    pass


class HCpmpmReedsSheppStateSpace(_StraightLinePlanner):
    pass
