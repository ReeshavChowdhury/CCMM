# Artifact: Automorphism-Fused Ciphertext Matrix Transpose (AFCMT)

This repository contains the C++ artifact evaluating the **AFCMT** architecture and Rectangular **Pre-CMT Accumulation** for the CKKS scheme.

## Artifact Structure

1. **`openfhe-development/`**: The modified OpenFHE C++ cryptographic library. We integrated our AFCMT Transpose algorithm and memory layout optimizations natively into the core OpenFHE framework.
   - **Modified Files**: 
     - `src/pke/lib/scheme/ckksrns/ckksrns-cmt.cpp`
     - `src/pke/include/scheme/ckksrns/ckksrns-cmt.h`
2. **`benchmarks/`**: Contains the reproducible evaluation scripts, isolated from the OpenFHE source tree with their own CMake configuration.
   - `levelA_B_bench.cpp`: Reproduces the Square CCMM isolated (Level A) and end-to-end (Level B) speedup tables.
   - `rect_bench.cpp`: Reproduces the Rectangular CCMM speedup (depth $B \in \{1, 2, 4, 8\}$) driven by $O(1)$ Deferred Pre-CMT accumulation.

---

## 1. Compilation Instructions

We provide an automated script to cleanly build and install OpenFHE to a local directory (`install/`), and then automatically configure and build the external benchmarks against it.

```bash
# Run the automated build script
./build.sh
```

Alternatively, to compile manually:
```bash
# 1. Compile and install OpenFHE locally
cd openfhe-development
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../install -DWITH_OPENMP=ON -DBUILD_EXAMPLES=OFF
make -j$(nproc) install
cd ../..

# 2. Compile the Benchmarks
cd benchmarks
mkdir build && cd build
cmake ..
make -j$(nproc)
```

---

## 2. Evaluation & Testing

Once compiled, you can directly execute the benchmark binaries from the `benchmarks/build/` directory. Be sure to point your library path to the local install directory so it uses the modified OpenFHE library and not a system-installed version.

```bash
export LD_LIBRARY_PATH=$PWD/install/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=20
```

### Run Level A and Level B Benchmarks (Square CCMM)
This benchmark mathematically locks the polynomial multiplication kernel (`FusedPPMM4`) across all evaluated variants to ensure the $1.25\times-2.7\times$ speedup represents the true, isolated performance gain generated strictly by the AFCMT transpose module.

```bash
./benchmarks/build/levelA_B_bench
```

### Run Rectangular CCMM Benchmark (Pre-CMT Accumulation)
This validates the Deferred Transpose theorem for rectangular CCMM of depth $B$. It strictly compares the Standard approach ($3B$ transposes) against our Deferred Pre-CMT approach ($B+2$ transposes).

```bash
./benchmarks/build/rect_bench
```
