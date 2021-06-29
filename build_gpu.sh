#!/bin/bash
set -e

module load lqcd/gpu/chroma/single-cuda11-qdpjit-nompi
module load lqcd/gpu/cmake/3.20.0

if [ ! -d build_gpu ]; then
    mkdir build_gpu
fi

pushd build_gpu
cmake -DCMAKE_INSTALL_PREFIX=../install/gpu \
      -DCMAKE_CXX_FLAGS="-O3 -g -DBUILD_JIT_CONTRACTION_KERNELS" ..
make -j 8 install
popd
