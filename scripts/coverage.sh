#!/usr/bin/env bash
#
# Measure and view test coverage for mphys_lib using LLVM source-based coverage.
#
# Usage:
#   ./scripts/coverage.sh            # build, run tests, print summary, open HTML report
#   ./scripts/coverage.sh --no-open  # same, but don't open the browser
#
# Only the project's own code (src/ and include/mphys/) is reported — external
# dependencies and the test sources themselves are excluded.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

LLVM_BIN="/opt/homebrew/opt/llvm/bin"
PROFDATA="$LLVM_BIN/llvm-profdata"
COV="$LLVM_BIN/llvm-cov"

BUILD_DIR="build/mac-coverage"
TEST_BIN="$BUILD_DIR/tests/mphys_tests"
COV_DIR="$BUILD_DIR/coverage"
RAW_PROFILE="$COV_DIR/mphys_tests.profraw"
MERGED_PROFILE="$COV_DIR/mphys_tests.profdata"
HTML_DIR="$COV_DIR/html"

# Report only on our own source. llvm-cov keeps files whose path matches a
# positive arg and drops anything matching -ignore-filename-regex.
SOURCES=("src" "include/mphys")
IGNORE_REGEX='(external/|/tests/|test_main)'

OPEN_REPORT=1
for arg in "$@"; do
    case "$arg" in
        --no-open) OPEN_REPORT=0 ;;
        *) echo "unknown arg: $arg" >&2; exit 1 ;;
    esac
done

echo "==> Configuring (coverage preset)"
cmake --preset mac-coverage >/dev/null

echo "==> Building tests"
cmake --build --preset mac-coverage

mkdir -p "$COV_DIR"

echo "==> Running tests"
LLVM_PROFILE_FILE="$RAW_PROFILE" "$TEST_BIN"

echo "==> Merging profile data"
"$PROFDATA" merge -sparse "$RAW_PROFILE" -o "$MERGED_PROFILE"

echo "==> Generating HTML report"
"$COV" show "$TEST_BIN" \
    -instr-profile="$MERGED_PROFILE" \
    -format=html \
    -output-dir="$HTML_DIR" \
    -ignore-filename-regex="$IGNORE_REGEX" \
    -show-line-counts-or-regions \
    "${SOURCES[@]}"

echo
echo "================= Coverage summary ================="
"$COV" report "$TEST_BIN" \
    -instr-profile="$MERGED_PROFILE" \
    -ignore-filename-regex="$IGNORE_REGEX" \
    "${SOURCES[@]}"
echo "===================================================="
echo
echo "HTML report: $ROOT/$HTML_DIR/index.html"

if [[ "$OPEN_REPORT" == "1" ]]; then
    open "$HTML_DIR/index.html"
fi
