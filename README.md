# Artifact: Fast Encrypted Matrix Multiplication (AFCMT)

This folder contains the code to run the experiments from our paper on making CKKS encrypted matrix multiplication (CCMM) faster.

We modified the standard OpenFHE library to include our faster "AFCMT" method.

## Code Structure

1. **`openfhe-development/`**: The modified OpenFHE library with our AFCMT code added inside.
2. **`benchmarks/`**: The scripts to run the experiments and print the results shown in the paper's tables.
   - `precision_test.cpp`: Checks that our method gets the exact same answers as the original method.
   - `levelA_B_bench.cpp`: Measures the speed of the isolated transpose step and the full matrix multiplication.
   - `rect_bench.cpp`: Measures the speed when multiplying rectangular matrices using our "Deferred Accumulation" trick.

---

## 1. How to Compile

Run the automated script to build OpenFHE and compile the benchmarks:

```bash
# Run the automated build script
./build.sh
```

---

## 2. How to Run the Experiments

After compiling, you can run the programs in the `benchmarks/build/` folder. First, set your library path so the system uses our modified OpenFHE:

```bash
export LD_LIBRARY_PATH=$PWD/install/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=20
```

### Reproduce Table 4: Pointwise Verification
This runs the precision test to prove that our fast method produces the exact same numbers as the original method with 0.00 difference.
```bash
./benchmarks/build/precision_test
```

### Reproduce Table 1 (Isolated CMT) and Table 2 (End-to-End CCMM)
This runs the main speed comparisons for square matrices.
- **Table 1** (Isolated CMT Benchmark): Compares the speed of just the transpose step.
- **Table 2** (End-to-End CCMM Benchmark): Compares the total time of the full encrypted matrix multiplication.
```bash
./benchmarks/build/levelA_B_bench
```

### Reproduce Table 5 (Pre-CMT Accumulation)
This runs the benchmark for rectangular matrices (where depth B is 1, 2, 4, or 8). It compares the standard approach to our new "Deferred" approach.
```bash
./benchmarks/build/rect_bench
```
