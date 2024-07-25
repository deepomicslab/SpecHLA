#!/bin/bash
set -euo pipefail

dir=$(cd `dirname $0`; pwd)

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
