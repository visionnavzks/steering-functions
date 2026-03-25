from pathlib import Path
import sys

from setuptools import setup


ROOT = Path(__file__).parent.resolve()
PYDEPS_DIR = ROOT / ".pydeps"

if PYDEPS_DIR.exists():
    sys.path.insert(0, str(PYDEPS_DIR))

from pybind11.setup_helpers import Pybind11Extension, build_ext


def glob_sources(pattern: str) -> list[str]:
    return sorted(str(path) for path in ROOT.glob(pattern))


extension = Pybind11Extension(
    "steering_functions._core",
    sorted(
        [
            str(ROOT / "src" / "python_bindings" / "module.cpp"),
            str(ROOT / "src" / "steering_path_lib" / "src" / "steering_path.cpp"),
            str(ROOT / "src" / "steering_functions" / "src" / "base_state_space.cpp"),
            str(ROOT / "src" / "steering_functions" / "src" / "dubins_state_space" / "dubins_state_space.cpp"),
            str(ROOT / "src" / "steering_functions" / "src" / "reeds_shepp_state_space" / "reeds_shepp_state_space.cpp"),
            str(ROOT / "src" / "steering_functions" / "src" / "utilities" / "utilities.cpp"),
            *glob_sources("src/steering_functions/src/hc_cc_core/*.cpp"),
            *glob_sources("src/steering_functions/src/cc_dubins_state_space/*.cpp"),
            *glob_sources("src/steering_functions/src/cc_reeds_shepp_state_space/*.cpp"),
            *glob_sources("src/steering_functions/src/hc_reeds_shepp_state_space/*.cpp"),
        ]
    ),
    include_dirs=[
        str(ROOT / "include"),
        str(ROOT / "src" / "steering_functions" / "include"),
        str(ROOT / "src" / "steering_path_lib" / "include"),
    ],
    cxx_std=17,
)


setup(
    name="steering-functions",
    version="0.1.0",
    description="Python package for steering functions path planning",
    package_dir={"": "python"},
    packages=["steering_functions"],
    ext_modules=[extension],
    cmdclass={"build_ext": build_ext},
    python_requires=">=3.9",
)
