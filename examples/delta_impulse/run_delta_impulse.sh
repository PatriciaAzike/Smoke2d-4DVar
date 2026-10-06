#!/bin/bash
# Run the instantaneous delta-impulse verification entirely inside this directory.

set -euo pipefail

case_dir="$(cd "$(dirname "$0")" && pwd)"
default_exe="$case_dir/../../../smoke2d-4DVar-build/examples/delta_impulse/delta_impulse"
exe="${DELTA_IMPULSE_EXE:-$default_exe}"

if [ ! -x "$exe" ]; then
    echo "delta_impulse executable not found: $exe" >&2
    echo "build target delta_impulse first" >&2
    exit 1
fi

for input in adjoint_options.ini forward_options.ini model_options.ini gauges.data; do
    if [ ! -f "$case_dir/$input" ]; then
        echo "missing input file: $case_dir/$input" >&2
        exit 1
    fi
done

mkdir -p "$case_dir/adjoint" "$case_dir/forward" "$case_dir/model" \
         "$case_dir/panels"


for run_dir in adjoint forward model; do
    cp -p "$case_dir/gauges.data" "$case_dir/$run_dir/gauges.data"
done

cd "$case_dir"
exec "$exe" "$@"
