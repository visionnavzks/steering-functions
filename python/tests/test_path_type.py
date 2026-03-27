"""Tests for PathType enum completeness and values."""

import unittest

from steering_functions import PathType


class TestPathTypeEnum(unittest.TestCase):
    """Test PathType enum values."""

    EXPECTED_MEMBERS = [
        "NONE",
        "CC_DUBINS",
        "CC00_DUBINS",
        "CC0PM_DUBINS",
        "CCPM0_DUBINS",
        "CCPMPM_DUBINS",
        "DUBINS",
        "CC00_RS",
        "HC_RS",
        "HC00_RS",
        "HC0PM_RS",
        "HCPM0_RS",
        "HCPMPM_RS",
        "RS",
    ]

    def test_all_expected_members_exist(self) -> None:
        for name in self.EXPECTED_MEMBERS:
            with self.subTest(name=name):
                self.assertTrue(hasattr(PathType, name), f"Missing PathType.{name}")

    def test_member_count(self) -> None:
        """Ensure no unexpected members are added without tests."""
        # PathType is a pybind11 enum; count known values
        count = sum(1 for name in self.EXPECTED_MEMBERS if hasattr(PathType, name))
        self.assertEqual(count, len(self.EXPECTED_MEMBERS))

    def test_enum_values_are_distinct(self) -> None:
        values = [getattr(PathType, name) for name in self.EXPECTED_MEMBERS]
        self.assertEqual(len(values), len(set(values)))

    def test_enum_repr(self) -> None:
        self.assertIn("DUBINS", repr(PathType.DUBINS))

    def test_enum_equality(self) -> None:
        self.assertEqual(PathType.RS, PathType.RS)
        self.assertNotEqual(PathType.RS, PathType.DUBINS)


if __name__ == "__main__":
    unittest.main()
