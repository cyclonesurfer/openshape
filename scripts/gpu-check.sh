#!/usr/bin/env bash
# Builds OpenShape, runs the test suite, the headless demo and the frame-time
# benchmark on this machine's GPU, and collects everything in one folder.
#
#   scripts/gpu-check.sh [out-dir]        (default: gpu-check-out)
#
# Prints PASS/FAIL per step and writes <out-dir>/report.txt (paste that back),
# the demo PNGs in <out-dir>/demo/, and <out-dir>.tar.gz (next to <out-dir>)
# with all of it.
# See docs/RUNNING.md for prerequisites.

set -uo pipefail
cd "$(dirname "$0")/.."

out="${1:-gpu-check-out}"
rm -rf "$out" "$out.tar.gz"
mkdir -p "$out"
report="$out/report.txt"
export OPENSHAPE_WGPU_BACKEND="${OPENSHAPE_WGPU_BACKEND:-vulkan}"
export RUST_LOG="${RUST_LOG:-openshape_render=info}"

say() { echo "$*" | tee -a "$report"; }

say "== OpenShape GPU check"
say "date:     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
say "commit:   $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
say "kernel:   $(uname -srm)"
say "os:       $(. /etc/os-release 2>/dev/null && echo "$PRETTY_NAME")"
say "rustc:    $(rustc --version 2>/dev/null || echo 'not found')"
say "backend:  OPENSHAPE_WGPU_BACKEND=$OPENSHAPE_WGPU_BACKEND"
say ""
say "== Vulkan devices (vulkaninfo --summary)"
if command -v vulkaninfo >/dev/null; then
    vulkaninfo --summary 2>/dev/null \
        | grep -E 'deviceName|deviceType|driverName|driverInfo|apiVersion' \
        | sed 's/^[[:space:]]*/  /' | tee -a "$report"
else
    say "  vulkaninfo not installed (package vulkan-tools)"
fi
say ""

failed=0
step() {
    local name=$1
    shift
    echo "-- $name: running"
    if "$@" >"$out/$name.log" 2>&1; then
        say "PASS  $name"
    else
        say "FAIL  $name   (log: $out/$name.log)"
        tail -n 30 "$out/$name.log" | sed 's/^/      /' >>"$report"
        failed=1
    fi
}

step build cargo build --release --workspace --all-targets
step tests cargo test --release --workspace
step demo cargo run --release -p openshape -- --headless-demo "$out/demo"
step bench cargo run --release -p openshape -- --bench "$out/bench" 600
step render-timing cargo test --release -p openshape-render --test offscreen \
    thousand_instances_timing -- --nocapture

say ""
say "== Headless demo"
grep -E '^(adapter:|\[PASS\]|\[FAIL\])' "$out/demo.log" 2>/dev/null | tee -a "$report"
say ""
say "== Benchmark"
cat "$out/bench/bench.txt" 2>/dev/null | tee -a "$report"
say ""
say "== Render timing (1000 instances)"
grep -E '^1000 x' "$out/render-timing.log" 2>/dev/null | tee -a "$report"
say ""
say "== Test totals"
grep -E '^test result:' "$out/tests.log" 2>/dev/null \
    | awk '{p+=$4; f+=$6} END {printf "passed %d, failed %d\n", p, f}' | tee -a "$report"

if grep -qiE '^adapter:.*(llvmpipe|lavapipe|softpipe)' "$out/demo.log" 2>/dev/null; then
    say ""
    say "WARNING: rendering ran on Mesa's software driver, not the GPU."
    say "         Check that RADV is installed and that /dev/dri/renderD* is accessible."
fi

say ""
if [ "$failed" = 0 ]; then
    say "RESULT: all steps passed"
else
    say "RESULT: some steps FAILED (see above)"
fi
# Last, so the report inside the archive has the RESULT line.
tar czf "$out.tar.gz" -C "$(dirname "$out")" "$(basename "$out")"
echo "Report: $report   Images: $out/demo/   Archive: $out.tar.gz"
exit "$failed"
