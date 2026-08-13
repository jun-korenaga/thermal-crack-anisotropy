#!/usr/bin/env bash
set -euo pipefail

code_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
make -C "$code_dir" all
mkdir -p "$code_dir/results"
"$code_dir/build/benchmark_kanaun09" -P360 -Q64 -N800 -W0.4 \
  > "$code_dir/results/benchmark.csv"
echo "Wrote $code_dir/results/benchmark.csv"
