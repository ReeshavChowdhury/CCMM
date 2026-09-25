#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-cmt.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <iomanip>

using namespace lbcrypto;
using Clock = std::chrono::high_resolution_clock;

double MsSince(Clock::time_point t) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
}

int main(int argc, char* argv[]) {
    uint32_t N = 8192;
    uint32_t d = 32; 
    std::vector<uint32_t> Bs = {1, 2, 4, 8};

    std::cout << "B\tStandard_ms\tDeferred_ms\tSpeedup\n";

    CMTContext ctx = MakeCMTContext(N, d);
    std::vector<std::vector<int64_t>> M1(d, std::vector<int64_t>(d, 1));
    auto encA = EncryptMatrixColumns(ctx, M1);
    auto encB = EncryptMatrixColumns(ctx, M1);

    for (uint32_t B : Bs) {
        ctx.afcmtMode = AFCMTMode::AFCMT_FUSED_RADIX2;
        
        auto run_std = [&]() {
            double total_ms = 0;
            int reps = 3;
            for(int w=0; w<1; ++w) {
                CCMMTimings tm;
                for(uint32_t b=0; b<B; ++b) RunCCMM(ctx, encA, encB, CMTVariant::AFCMT, &tm, false, true);
            }
            for(int r=0; r<reps; ++r) {
                CCMMTimings tm;
                auto t0 = Clock::now();
                for(uint32_t b=0; b<B; ++b) RunCCMM(ctx, encA, encB, CMTVariant::AFCMT, &tm, false, true);
                total_ms += MsSince(t0);
            }
            return total_ms / reps;
        };

        auto run_def = [&]() {
            double total_ms = 0;
            int reps = 3;
            for(int w=0; w<1; ++w) {
                CCMMTimings tm;
                RunRectangularCCMM(ctx, encA, encB, B, CMTVariant::AFCMT, &tm);
            }
            for(int r=0; r<reps; ++r) {
                CCMMTimings tm;
                auto t0 = Clock::now();
                RunRectangularCCMM(ctx, encA, encB, B, CMTVariant::AFCMT, &tm);
                total_ms += MsSince(t0);
            }
            return total_ms / reps;
        };

        double t_std = run_std();
        double t_def = run_def();
        double speedup = t_std / t_def;

        std::cout << B << "\t" 
                  << std::fixed << std::setprecision(2) << t_std << "\t\t" 
                  << t_def << "\t\t" 
                  << speedup << "x\n" << std::flush;
    }
    return 0;
}
