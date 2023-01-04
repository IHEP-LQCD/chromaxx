#!/bin/bash
set -e

module load lqcd/x86/chroma/double-qphix-qdpxx
module load lqcd/x86/cmake/3.20.0

source_dir="$(pwd)/../"
build_dir="$(pwd)/build/x86"
install_dir="$(pwd)/install/x86"

if [ -d ${build_dir} ]; then
    rm -rf ${build_dir}
fi
mkdir ${build_dir}

pushd ${build_dir}
cmake -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_C_COMPILER=mpiicc \
      -DCMAKE_INSTALL_PREFIX=${install_dir} \
      -DCMAKE_CXX_FLAGS="-O3 -g" ${source_dir}
make -j 8 install
popd
