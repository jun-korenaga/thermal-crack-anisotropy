#!/usr/bin/env bash
set -euo pipefail

release_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
executable="$release_root/build/emt_effective_stiffness"
output_dir="$release_root/results"
mkdir -p "$output_dir"
output="$output_dir/appendix_a5_effective_stiffness.csv"

if [[ ! -x "$executable" ]]; then
  echo "Build the project first with scripts/build.sh" >&2
  exit 1
fi

ntheta=${NTHETA:-3200}
aspect_values=${ASPECT_VALUES:-"0.1 1e-2 1e-3"}
dip_values=${DIP_VALUES:-"$(seq 0 1 90)"}
: > "$output"
first=1
for aspect in $aspect_values; do
  for dip in $dip_values; do
    result=$("$executable" --configuration horizontal --aspect "$aspect" \
      --porosity 0.0209 --fluid-bulk 2.9 --dip "$dip" --ntheta "$ntheta")
    if [[ $first -eq 1 ]]; then
      echo "$result" > "$output"
      first=0
    else
      echo "$result" | tail -n 1 >> "$output"
    fi
  done
done
echo "Wrote $output"
