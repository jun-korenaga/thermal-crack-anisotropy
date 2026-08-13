#!/usr/bin/env bash
set -euo pipefail

release_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
executable="$release_root/build/emt_withers_benchmark"
output_dir="$release_root/results"
mkdir -p "$output_dir"
output="$output_dir/appendix_a2_withers.csv"

if [[ ! -x "$executable" ]]; then
  echo "Build the project first with scripts/build.sh" >&2
  exit 1
fi

# These defaults reproduce the parameter grid in manuscript Figure A2 and can
# be computationally expensive. Override NTHETA_VALUES for a smaller check.
ntheta_values=${NTHETA_VALUES:-"100 200 400 800 1600 3200"}
weight_values=${WEIGHT_VALUES:-"0.0 0.2 0.4"}
aspect_values=${ASPECT_VALUES:-"1e5 1e4 1e3 1e2 1e1 1 1e-1 1e-2 1e-3 1e-4 1e-5"}

echo "aspect,ntheta,weight,integration_points,relative_rms_error,max_absolute_major_asymmetry,max_level_hits" > "$output"
for aspect in $aspect_values; do
  for ntheta in $ntheta_values; do
    for weight in $weight_values; do
      "$executable" --aspect "$aspect" --ntheta "$ntheta" --weight "$weight" |
        tail -n 1 >> "$output"
    done
  done
done
echo "Wrote $output"
