#!/bin/bash
set -e

module load lqcd/arm/chroma/double-qopqdp-qdpxx
module load lqcd/arm/cmake/3.20.0

source_dir="$(pwd)"
build_dir="$(pwd)/build/arm_double"
install_dir="$(pwd)/install/arm_double"

if [ -d ${build_dir} ]; then
    rm -rf ${build_dir}
fi
mkdir ${build_dir}

pushd ${build_dir}
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_INSTALL_PREFIX=${install_dir} \
      -DCMAKE_CXX_FLAGS="-O3 -g" ${source_dir}
make -j 8 install
popd
