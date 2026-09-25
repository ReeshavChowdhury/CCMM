#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-cmt.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace lbcrypto;
using Clock = std::chrono::high_resolution_clock;

double MsSince(Clock::time_point t) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
}

int main(int argc, char* argv[]) {
    uint32_t N = 8192;
    std::vector<uint32_t> ds = {8, 16, 32, 64};

    std::cout << "=========================================================\n";
    std::cout << "Level A — CMT-to-CMT comparison (Isolated Transpose)\n";
    std::cout << "=========================================================\n";
    std::cout << "N\td\tT_Ref(ms)\tT_FCMT(ms)\tT_AFCMT(ms)\tS_fusion\tS_kernel\n";

    for (uint32_t d : ds) {
        CMTContext ctx = MakeCMTContext(N, d);
        std::vector<Ciphertext<DCRTPoly>> A(d);
        for(uint32_t i=0; i<d; ++i){
            std::vector<double> v(N/2, 1.0);
            A[i] = ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(v));
        }

        auto run_cmt = [&](CMTVariant variant, AFCMTMode mode) {
            ctx.afcmtMode = mode;
            double total_ms = 0;
            int reps = 3;
            for(int w=0; w<2; ++w) {
                CCMMTimings dummy;
                RunCMT(ctx, A, variant, &dummy.cmt[0]);
            }
            for(int r=0; r<reps; ++r) {
                CCMMTimings tm;
                auto t0 = Clock::now();
                RunCMT(ctx, A, variant, &tm.cmt[0]);
                total_ms += MsSince(t0);
            }
            return total_ms / reps;
        };

        double t_ref = run_cmt(CMTVariant::CKL_REFERENCE, AFCMTMode::AFCMT_GENERIC);
        double t_fcmt = run_cmt(CMTVariant::FCMT, AFCMTMode::AFCMT_GENERIC);
        double t_afcmt = run_cmt(CMTVariant::AFCMT, AFCMTMode::AFCMT_FUSED_RADIX2);

        double s_fusion = t_ref / t_fcmt;
        double s_kernel = t_fcmt / t_afcmt;

        std::cout << N << "\t" << d << "\t" 
                  << std::fixed << std::setprecision(2) << t_ref << "\t\t" 
                  << t_fcmt << "\t\t" << t_afcmt << "\t\t"
                  << s_fusion << "x\t\t" << s_kernel << "x\n" << std::flush;
    }

    std::cout << "\n=========================================================\n";
    std::cout << "Level B — complete CKL CCMM comparison (End-to-End)\n";
    std::cout << "=========================================================\n";
    std::cout << "N\td\tCCMM_Ref(ms)\tCCMM_FCMT(ms)\tCCMM_AFCMT(ms)\tOverall_Speedup\n";

    for (uint32_t d : ds) {
        CMTContext ctx = MakeCMTContext(N, d);
        std::vector<std::vector<int64_t>> M1(d, std::vector<int64_t>(d, 1));
        std::vector<std::vector<int64_t>> M2(d, std::vector<int64_t>(d, 1));
        auto encA = EncryptMatrixColumns(ctx, M1);
        auto encB = EncryptMatrixColumns(ctx, M2);

        auto run_ccmm = [&](CMTVariant variant, AFCMTMode mode) {
            ctx.afcmtMode = mode;
            double total_ms = 0;
            int reps = 3;
            for(int w=0; w<2; ++w) {
                CCMMTimings dummy;
                RunCCMM(ctx, encA, encB, variant, &dummy, false, true);
            }
            for(int r=0; r<reps; ++r) {
                CCMMTimings tm;
                auto t0 = Clock::now();
                RunCCMM(ctx, encA, encB, variant, &tm, false, true);
                total_ms += MsSince(t0);
            }
            return total_ms / reps;
        };

        double t_ref = run_ccmm(CMTVariant::CKL_REFERENCE, AFCMTMode::AFCMT_GENERIC);
        double t_fcmt = run_ccmm(CMTVariant::FCMT, AFCMTMode::AFCMT_GENERIC);
        double t_afcmt = run_ccmm(CMTVariant::AFCMT, AFCMTMode::AFCMT_FUSED_RADIX2);

        double speedup = t_ref / t_afcmt;

        std::cout << N << "\t" << d << "\t" 
                  << std::fixed << std::setprecision(2) << t_ref << "\t\t" 
                  << t_fcmt << "\t\t" << t_afcmt << "\t\t"
                  << speedup << "x\n" << std::flush;
    }

    return 0;
}
