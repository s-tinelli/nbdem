#!/bin/sh
# Clones Stim at the tested commit, applies the non-binary DEM patch and builds the command line tool.
# Usage: ./setup.sh [target_dir] [--tests]
#   target_dir  where Stim is cloned (default: ./Stim)
#   --tests     also build and run the patch's unit tests (needs GTest, see README)
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
STIM_COMMIT=ffe0d043598ed5df0851a3e6b4920359a23fafbd
DST="$HERE/Stim"
TESTS=0
for a in "$@"; do
  case "$a" in
    --tests) TESTS=1 ;;
    *) DST="$a" ;;
  esac
done
JOBS="${JOBS:-4}"

git clone https://github.com/quantumlib/Stim.git "$DST"
cd "$DST"
git checkout "$STIM_COMMIT"
git apply "$HERE/nb_dem.patch"
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make stim -j"$JOBS"
if [ "$TESTS" = 1 ]; then
  make stim_test -j"$JOBS"
  (cd .. && ./build/out/stim_test --gtest_filter='nb_dem.*')
fi
echo "OK: $DST/build/out/stim analyze_errors_nb --in circuit.stim"
