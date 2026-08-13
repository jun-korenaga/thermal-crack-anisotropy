#!/usr/bin/env bash
set -euo pipefail

release_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
if command -v cmake >/dev/null 2>&1; then
  cmake -S "$release_root" -B "$release_root/build" -DCMAKE_BUILD_TYPE=Release
  cmake --build "$release_root/build" --parallel
  ctest --test-dir "$release_root/build" --output-on-failure
else
  echo "CMake was not found; using the included GNU Make fallback."
  make -C "$release_root" all
  make -C "$release_root" test
fi
