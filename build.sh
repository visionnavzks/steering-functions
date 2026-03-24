#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-${ROOT_DIR}/build}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
BUILD_TESTING="${BUILD_TESTING:-OFF}"
RUN_TESTS="${RUN_TESTS:-0}"

cmake_args=(
    -S "${ROOT_DIR}"
    -B "${BUILD_DIR}"
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
    -DBUILD_TESTING="${BUILD_TESTING}"
)

if command -v ninja >/dev/null 2>&1; then
    cmake_args+=(-G Ninja)
fi

if [ "$#" -gt 0 ]; then
    cmake_args+=("$@")
fi

cmake "${cmake_args[@]}"
cmake --build "${BUILD_DIR}" --parallel

if [ "${RUN_TESTS}" = "1" ]; then
    ctest --test-dir "${BUILD_DIR}" --output-on-failure
fi