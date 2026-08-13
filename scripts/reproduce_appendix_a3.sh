#!/usr/bin/env bash
set -euo pipefail

release_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
executable="$release_root/build/emt_guo_comparison"
output_dir="$release_root/results"
mkdir -p "$output_dir"

if [[ ! -x "$executable" ]]; then
  echo "Build the project first with scripts/build.sh" >&2
  exit 1
fi

ntheta=${NTHETA:-800}
angle_step=${ANGLE_STEP:-1}
"$executable" --aspect 1e-2 --ntheta "$ntheta" --weight 0.4 \
  --angle-step "$angle_step" > "$output_dir/appendix_a3_guo_comparison.csv"
echo "Wrote $output_dir/appendix_a3_guo_comparison.csv"
