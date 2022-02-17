#!/bin/bash
set -e

module load myqcd/x86_gpu/chroma/double-cuda11-qdpjit-llvm13
module load lqcd/gpu/cmake/3.20.0

build_dir="build_gpu_new"
install_dir="gpu_new"
if [ ! -d $build_dir ]; then
    mkdir $build_dir
fi

pushd $build_dir
cmake -DCMAKE_INSTALL_PREFIX=../install/$install_dir \
	  -DChroma_DIR=${CHROMA_HOME}/lib/cmake/Chroma \
	  -DQDPXX_DIR=/dg_hpc/LQCD/sunwei/software/lqcd202111/x86_gpu/chroma-gcc7.3.1-cuda11.0/install/qdpjit-double-llvm13/lib/cmake/QDPXX \
	  -DQMP_DIR=/dg_hpc/LQCD/sunwei/software/lqcd202111/x86_gpu/chroma-gcc7.3.1-cuda11.0/install/qmp-2.5.4-openmpi-gcc7/lib/cmake/QMP \
      -DCMAKE_CXX_FLAGS="-O3 -g -DBUILD_JIT_CONTRACTION_KERNELS" ..
make -j 8 install
popd
