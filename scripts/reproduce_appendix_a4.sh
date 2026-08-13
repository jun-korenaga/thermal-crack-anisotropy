#!/usr/bin/env bash
set -euo pipefail

release_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
executable="$release_root/build/emt_random_orientations"
output_dir="$release_root/results"
mkdir -p "$output_dir"

if [[ ! -x "$executable" ]]; then
  echo "Build the project first with scripts/build.sh" >&2
  exit 1
fi

samples=${SAMPLES:-10000}
ntheta=${NTHETA:-400}
seed=${SEED:-20260813}
output_every=${OUTPUT_EVERY:-10}
"$executable" --samples "$samples" --output-every "$output_every" --seed "$seed" \
  --aspect 1e-2 --ntheta "$ntheta" --weight 0.4 \
  > "$output_dir/appendix_a4_random_orientations.csv"
echo "Wrote $output_dir/appendix_a4_random_orientations.csv"
