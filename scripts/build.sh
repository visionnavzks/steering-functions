#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_PYTHON="${BUILD_PYTHON:-0}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
BUILD_TESTING="${BUILD_TESTING:-OFF}"
RUN_TESTS="${RUN_TESTS:-0}"
PYTHON_BIN="${PYTHON_BIN:-${PYTHON_EXECUTABLE:-$(command -v python3)}}"
RUN_PYTHON_TESTS="${RUN_PYTHON_TESTS:-0}"
RUN_PYTHON_DEMO="${RUN_PYTHON_DEMO:-0}"

if [ "${BUILD_PYTHON}" = "1" ]; then
    BUILD_DIR="${BUILD_DIR:-${ROOT_DIR}/build-python}"
else
    BUILD_DIR="${BUILD_DIR:-${ROOT_DIR}/build}"
fi

cmake_args=(
    -S "${ROOT_DIR}"
    -B "${BUILD_DIR}"
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
    -DBUILD_TESTING="${BUILD_TESTING}"
)

if [ "${BUILD_PYTHON}" = "1" ]; then
    cmake_args+=(
        -DBUILD_PYTHON_BINDINGS=ON
        -DPython3_EXECUTABLE="${PYTHON_BIN}"
    )
fi

if command -v ninja >/dev/null 2>&1; then
    cmake_args+=(-G Ninja)
fi

if [ "$#" -gt 0 ]; then
    cmake_args+=("$@")
fi

cmake "${cmake_args[@]}"

if [ "${BUILD_PYTHON}" = "1" ]; then
    cmake --build "${BUILD_DIR}" --parallel --target steering_functions_python
    export PYTHONPATH="${ROOT_DIR}/python${PYTHONPATH:+:${PYTHONPATH}}"
else
    cmake --build "${BUILD_DIR}" --parallel
fi

if [ "${RUN_TESTS}" = "1" ]; then
    ctest --test-dir "${BUILD_DIR}" --output-on-failure
fi

if [ "${RUN_PYTHON_TESTS}" = "1" ]; then
    "${PYTHON_BIN}" -m unittest discover -s "${ROOT_DIR}/python/tests" -v
fi

if [ "${RUN_PYTHON_DEMO}" = "1" ]; then
    "${PYTHON_BIN}" "${ROOT_DIR}/python/demo_interactive.py"
fi