"""Tests for State and Control data types."""

import math
import unittest

from steering_functions import Control, State


class TestStateCreation(unittest.TestCase):
    """Test State construction and attribute access."""

    def test_default_constructor(self) -> None:
        state = State()
        self.assertEqual(state.x, 0.0)
        self.assertEqual(state.y, 0.0)
        self.assertEqual(state.theta, 0.0)
        self.assertEqual(state.kappa, 0.0)

    def test_positional_constructor(self) -> None:
        state = State(1.0, 2.0, 0.5, 0.25)
        self.assertEqual(state.x, 1.0)
        self.assertEqual(state.y, 2.0)
        self.assertEqual(state.theta, 0.5)
        self.assertEqual(state.kappa, 0.25)

    def test_keyword_constructor(self) -> None:
        state = State(x=3.0, y=-1.0, theta=math.pi, kappa=-0.5)
        self.assertEqual(state.x, 3.0)
        self.assertEqual(state.y, -1.0)
        self.assertAlmostEqual(state.theta, math.pi)
        self.assertEqual(state.kappa, -0.5)

    def test_default_kappa_is_zero(self) -> None:
        state = State(1.0, 2.0, 0.5)
        self.assertEqual(state.kappa, 0.0)

    def test_readwrite_attributes(self) -> None:
        state = State()
        state.x = 5.0
        state.y = -3.0
        state.theta = 1.2
        state.kappa = 0.8
        state.sigma = 0.3
        state.d = 1.0
        state.s = 2.5
        state.vel = 1.5
        state.acc = 0.1
        state.time = 10.0
        state.fork_y = 0.5

        self.assertEqual(state.x, 5.0)
        self.assertEqual(state.y, -3.0)
        self.assertEqual(state.theta, 1.2)
        self.assertEqual(state.kappa, 0.8)
        self.assertEqual(state.sigma, 0.3)
        self.assertEqual(state.d, 1.0)
        self.assertEqual(state.s, 2.5)
        self.assertEqual(state.vel, 1.5)
        self.assertEqual(state.acc, 0.1)
        self.assertEqual(state.time, 10.0)
        self.assertEqual(state.fork_y, 0.5)

    def test_negative_coordinates(self) -> None:
        state = State(-10.0, -20.0, -math.pi, -1.0)
        self.assertEqual(state.x, -10.0)
        self.assertEqual(state.y, -20.0)
        self.assertAlmostEqual(state.theta, -math.pi)
        self.assertEqual(state.kappa, -1.0)

    def test_large_coordinates(self) -> None:
        state = State(1e6, 1e6, 0.0, 0.0)
        self.assertEqual(state.x, 1e6)
        self.assertEqual(state.y, 1e6)


class TestStateEquality(unittest.TestCase):
    """Test State equality operator."""

    def test_equal_states(self) -> None:
        a = State(1.0, 2.0, 0.5, 0.25)
        b = State(1.0, 2.0, 0.5, 0.25)
        self.assertEqual(a, b)

    def test_unequal_x(self) -> None:
        a = State(1.0, 2.0, 0.5, 0.25)
        b = State(1.1, 2.0, 0.5, 0.25)
        self.assertNotEqual(a, b)

    def test_unequal_y(self) -> None:
        a = State(1.0, 2.0, 0.5, 0.25)
        b = State(1.0, 2.1, 0.5, 0.25)
        self.assertNotEqual(a, b)

    def test_unequal_theta(self) -> None:
        a = State(1.0, 2.0, 0.5, 0.25)
        b = State(1.0, 2.0, 0.6, 0.25)
        self.assertNotEqual(a, b)

    def test_unequal_kappa(self) -> None:
        a = State(1.0, 2.0, 0.5, 0.25)
        b = State(1.0, 2.0, 0.5, 0.3)
        self.assertNotEqual(a, b)

    def test_default_states_are_equal(self) -> None:
        self.assertEqual(State(), State())


class TestStateRepresentation(unittest.TestCase):
    """Test State string representations."""

    def test_str_contains_state(self) -> None:
        state = State(1.0, 2.0, 0.5, 0.25)
        text = str(state)
        self.assertIn("State(", text)

    def test_repr_contains_state(self) -> None:
        state = State(1.0, 2.0, 0.5, 0.25)
        text = repr(state)
        self.assertIn("State(", text)

    def test_to_string(self) -> None:
        state = State(1.0, 2.0, 0.5, 0.25)
        text = state.to_string()
        self.assertIsInstance(text, str)
        self.assertTrue(len(text) > 0)


class TestControlCreation(unittest.TestCase):
    """Test Control construction and attributes."""

    def test_default_constructor(self) -> None:
        control = Control()
        self.assertEqual(control.delta_s, 0.0)
        self.assertEqual(control.kappa, 0.0)
        self.assertEqual(control.sigma, 0.0)

    def test_positional_constructor(self) -> None:
        control = Control(1.5, 0.5, 0.2)
        self.assertEqual(control.delta_s, 1.5)
        self.assertEqual(control.kappa, 0.5)
        self.assertEqual(control.sigma, 0.2)

    def test_keyword_constructor(self) -> None:
        control = Control(delta_s=2.0, kappa=-0.3, sigma=0.1)
        self.assertEqual(control.delta_s, 2.0)
        self.assertEqual(control.kappa, -0.3)
        self.assertEqual(control.sigma, 0.1)

    def test_readwrite_attributes(self) -> None:
        control = Control()
        control.delta_s = 3.0
        control.kappa = 0.7
        control.sigma = 0.4
        self.assertEqual(control.delta_s, 3.0)
        self.assertEqual(control.kappa, 0.7)
        self.assertEqual(control.sigma, 0.4)

    def test_negative_delta_s(self) -> None:
        control = Control(-2.0, 0.5, 0.0)
        self.assertEqual(control.delta_s, -2.0)


class TestControlRepresentation(unittest.TestCase):
    """Test Control string representations."""

    def test_str_contains_control(self) -> None:
        control = Control(1.5, 0.5, 0.0)
        text = str(control)
        self.assertIn("Control Segment", text)

    def test_repr_contains_control(self) -> None:
        control = Control(1.5, 0.5, 0.0)
        text = repr(control)
        self.assertIn("Control Segment", text)

    def test_to_string(self) -> None:
        control = Control(1.5, 0.5, 0.0)
        text = control.to_string()
        self.assertIsInstance(text, str)
        self.assertTrue(len(text) > 0)


if __name__ == "__main__":
    unittest.main()
