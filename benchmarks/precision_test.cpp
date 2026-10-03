#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-cmt.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace lbcrypto;

int main(int argc, char* argv[]) {
    std::vector<uint32_t> Ns = {8192, 16384, 32768};
    
    std::cout << "=========================================================\n";
    std::cout << "Precision Validation: AFCMT vs CKL_REFERENCE\n";
    std::cout << "=========================================================\n";
    std::cout << "N\td\tDiff_vs_Plain\tDiff_vs_CKL\tMatch_CKL?\n";

    for (uint32_t N : Ns) {
        for (uint32_t d = 8; d * d <= N / 2; d *= 2) {
            
            CMTContext ctx = MakeCMTContext(N, d, 2, 50);
            
            std::vector<std::vector<int64_t>> M(d, std::vector<int64_t>(d, 0));
            for(uint32_t r = 0; r < d; ++r) {
                for(uint32_t c = 0; c < d; ++c) {
                    M[r][c] = (rand() % 100) - 50;
                }
            }
            
            auto A = EncryptMatrixColumns(ctx, M);

            ctx.afcmtMode = AFCMTMode::AFCMT_GENERIC;
            CCMMTimings tm1;
            auto out_ref = RunCMT(ctx, A, CMTVariant::CKL_REFERENCE, &tm1.cmt[0]);
            
            ctx.afcmtMode = AFCMTMode::AFCMT_FUSED_RADIX2;
            CCMMTimings tm2;
            auto out_afcmt = RunCMT(ctx, A, CMTVariant::AFCMT, &tm2.cmt[0]);
            
            auto M_ref = DecryptMatrixColumns(ctx, out_ref);
            auto M_afcmt = DecryptMatrixColumns(ctx, out_afcmt);
            
            double max_diff = 0.0; double max_diff_pt = 0.0;
            
            for (uint32_t r = 0; r < d; ++r) {
                for (uint32_t c = 0; c < d; ++c) {
                    double diff = std::abs(M_ref[r][c] - M_afcmt[r][c]);
                    double diff_pt = std::abs(M[c][r] - M_afcmt[r][c]);
                    if (diff_pt > max_diff_pt) max_diff_pt = diff_pt;
                    if (diff > max_diff) {
                        max_diff = diff;
                    }
                }
            }
            
            bool match = max_diff < 1e-3;
            std::cout << N << "\t" << d << "\t" 
                      << std::scientific << std::setprecision(2) << max_diff_pt << "\t" << max_diff << "\t" 
                      << (match ? "YES" : "NO") << "\n" << std::flush;
        }
    }
    return 0;
}
