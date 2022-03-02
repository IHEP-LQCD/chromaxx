#!/bin/bash
set -e

module load lqcd/gpu/chroma/double-cuda11-qdpjit
#module load swlqcd/gpu/chroma/double-cuda11-qdpjit-llvm13
module load lqcd/gpu/cmake/3.20.0

source_dir="$(pwd)"
build_dir="$(pwd)/build/gpu_double"
install_dir="$(pwd)/install/gpu_double"

if [ ! -d ${build_dir} ]; then
    mkdir -p ${build_dir}
fi

pushd ${build_dir}
cmake -DCMAKE_INSTALL_PREFIX=${install_dir} \
	  -DCMAKE_CXX_COMPILER=mpicxx \
      -DCMAKE_CXX_FLAGS="-O3 -g -DBUILD_JIT_CONTRACTION_KERNELS" ${source_dir}
make -j 8 install
popd
