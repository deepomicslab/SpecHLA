#!/bin/bash
set -euo pipefail

# Ensure the dir is set relative to the actual path, links should be resolved
# so that this script can be linked somewhere else.
script_path=$(dirname $(realpath $0))
dir=$(cd $script_path; pwd)

# install spechap
mkdir -p $dir/bin/SpecHap/build
cd $dir/bin/SpecHap/build
cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
make

# install spechap
mkdir -p $dir/bin/extractHairs/build
cd $dir/bin/extractHairs/build
cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
make
