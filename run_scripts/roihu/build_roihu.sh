#!/bin/bash
# Builds on Roihu, inside a job (BUILD_TARGETS picks the programs, default dipole):
#   srun --account=lappi --partition=small --time=00:15:00 --cpus-per-task=8 ./run_scripts/roihu/build_roihu.sh

set -e
cd "$(dirname "$0")/../.."
source run_scripts/roihu/roihu_modules.sh
load_roihu_modules
module load spack/x86_64/v2026_03/Core/cmake/3.31.11

BUILD_TARGETS=${BUILD_TARGETS:-"dipole"}
mkdir -p build
cmake -S . -B build
cmake --build build -j"$(nproc)" --target $BUILD_TARGETS

echo "Build OK: $(for t in $BUILD_TARGETS; do printf 'build/bin/%s ' "$t"; done)"
