#!/bin/bash
set -e

module load lqcd/x86/chroma/double-qphix-qdpxx
module load lqcd/x86/cmake/3.20.0

if [ -d build_x86 ]; then
    rm -rf build_x86
fi
mkdir build_x86

pushd build_x86
cmake -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_C_COMPILER=mpiicc \
      -DCMAKE_INSTALL_PREFIX=../install/x86 \
      -DCMAKE_CXX_FLAGS="-O3 -g" ..
make -j 8 install
popd
