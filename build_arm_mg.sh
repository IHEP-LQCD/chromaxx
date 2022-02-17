#!/bin/bash
set -e

module load myqcd/arm/chroma/double-mgproto-qdpxx
module load lqcd/arm/cmake/3.20.0

if [ -d build_arm_mg ]; then
    rm -rf build_arm_mg
fi
mkdir build_arm_mg

pushd build_arm_mg
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_INSTALL_PREFIX=../install/arm_mg \
      -DCMAKE_CXX_FLAGS="-O3 -g" ..
make -j 8 install
popd
