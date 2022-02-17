#!/bin/bash
set -e

module load lqcd/arm/chroma/double-qopqdp-qdpxx
module load lqcd/arm/cmake/3.20.0

if [ -d build_arm ]; then
    rm -rf build_arm
fi
mkdir build_arm

pushd build_arm
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_INSTALL_PREFIX=../install/arm \
      -DCMAKE_CXX_FLAGS="-O3 -g" ..
make -j 8 install
popd
