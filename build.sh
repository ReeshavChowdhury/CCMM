#!/bin/bash
set -e

ROOT_DIR=$(pwd)
INSTALL_DIR=$ROOT_DIR/install

echo "========================================"
echo "1. Building and installing OpenFHE..."
echo "========================================"
cd openfhe-development
mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DWITH_OPENMP=ON -DBUILD_EXAMPLES=OFF -DBUILD_UNITTESTS=OFF -DBUILD_BENCHMARKS=OFF
make -j$(nproc) install

echo "========================================"
echo "2. Building Benchmarks..."
echo "========================================"
cd $ROOT_DIR/benchmarks
mkdir -p build && cd build
cmake ..
make -j$(nproc)

echo "========================================"
echo "Build Complete!"
echo "Binaries are located in: benchmarks/build/"
echo "Run them with:"
echo "  export LD_LIBRARY_PATH=$INSTALL_DIR/lib:\$LD_LIBRARY_PATH"
echo "  export OMP_NUM_THREADS=20"
echo "  ./benchmarks/build/levelA_B_bench"
echo "  ./benchmarks/build/rect_bench"
echo "========================================"
