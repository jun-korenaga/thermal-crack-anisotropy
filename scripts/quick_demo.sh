#!/usr/bin/env bash
set -euo pipefail

release_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_dir="$release_root/results/quick_demo"
mkdir -p "$output_dir"

"$release_root/build/emt_withers_benchmark" --aspect 0.01 --ntheta 100 \
  > "$output_dir/withers.csv"
"$release_root/build/emt_guo_comparison" --aspect 0.01 --ntheta 100 \
  --angle-step 15 > "$output_dir/guo.csv"
"$release_root/build/emt_random_orientations" --samples 20 --ntheta 40 \
  --output-every 5 --seed 20260813 > "$output_dir/random.csv"
"$release_root/build/emt_effective_stiffness" --configuration hexagonal \
  --aspect 0.01 --porosity 0.001 --ntheta 100 \
  > "$output_dir/effective_stiffness.csv"

echo "Wrote quick demonstration outputs under $output_dir"
