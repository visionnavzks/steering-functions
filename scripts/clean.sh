#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

dirs=(
    "${ROOT_DIR}/build"
    "${ROOT_DIR}/build-debug"
    "${ROOT_DIR}/build-python"
)

for dir in "${dirs[@]}"; do
    if [ -d "${dir}" ]; then
        echo "Removing ${dir}"
        rm -rf "${dir}"
    fi
done

echo "Clean complete."
