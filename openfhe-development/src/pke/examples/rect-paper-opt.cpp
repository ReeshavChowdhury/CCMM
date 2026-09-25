#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-cmt.h"
#include <iostream>
#include <chrono>

using namespace lbcrypto;
using Clock = std::chrono::high_resolution_clock;

double MsSince(Clock::time_point t) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
}

int main() {
    uint32_t N = 8192;
    
    // Experiment 1: 128 x 64 x 32 (d = 64)
    {
        uint32_t d = 64;
        std::cout << "\n=== EXPERIMENT 1: 128 x 64 x 32 (d=64, 2 blocks) ===" << std::endl;
        CMTContext ctx = MakeCMTContext(N, d);
        
        std::vector<Ciphertext<DCRTPoly>> A1(d), A2(d), B1(d);
        for(uint32_t i=0; i<d; ++i){
            std::vector<double> v(N/2, 1.0);
            Plaintext pt = ctx.cc->MakeCKKSPackedPlaintext(v);
            A1[i] = ctx.cc->Encrypt(ctx.pk, pt);
            A2[i] = ctx.cc->Encrypt(ctx.pk, pt);
            B1[i] = ctx.cc->Encrypt(ctx.pk, pt);
        }
        
        // Proper Warmup
        CCMMTimings dummy;
        for(int w=0; w<5; ++w) {
            RunCCMM(ctx, A1, B1, CMTVariant::AFCMT, &dummy, false);
        }
        
        double total_ms = 0;
        int reps = 10;
        for(int r=0; r<reps; ++r) {
            auto t0 = Clock::now();
            CCMMTimings tm1, tm2;
            auto B1_CMT = RunCMT(ctx, B1, CMTVariant::AFCMT, &tm1.cmt[0]);
            RunCCMM(ctx, A1, B1_CMT, CMTVariant::AFCMT, &tm1, true);
            RunCCMM(ctx, A2, B1_CMT, CMTVariant::AFCMT, &tm2, true);
            total_ms += MsSince(t0);
        }
        std::cout << "Optimized AFCMT: " << (total_ms / reps) / 1000.0 << " seconds." << std::endl;
    }
    
    // Experiment 2: 256 x 32 x 16 (d = 32)
    {
        uint32_t d = 32;
        std::cout << "\n=== EXPERIMENT 2: 256 x 32 x 16 (d=32, 8 blocks) ===" << std::endl;
        CMTContext ctx = MakeCMTContext(N, d);
        
        std::vector<std::vector<Ciphertext<DCRTPoly>>> A(8, std::vector<Ciphertext<DCRTPoly>>(d));
        std::vector<Ciphertext<DCRTPoly>> B1(d);
        for(uint32_t i=0; i<d; ++i){
            std::vector<double> v(N/2, 1.0);
            Plaintext pt = ctx.cc->MakeCKKSPackedPlaintext(v);
            B1[i] = ctx.cc->Encrypt(ctx.pk, pt);
            for(int b=0; b<8; ++b) {
                A[b][i] = ctx.cc->Encrypt(ctx.pk, pt);
            }
        }
        
        // Proper Warmup
        CCMMTimings dummy;
        for(int w=0; w<5; ++w) {
            RunCCMM(ctx, A[0], B1, CMTVariant::AFCMT, &dummy, false);
        }
        
        double total_ms = 0;
        int reps = 10;
        for(int r=0; r<reps; ++r) {
            auto t0 = Clock::now();
            CCMMTimings tm1, tm_rest;
            auto B1_CMT = RunCMT(ctx, B1, CMTVariant::AFCMT, &tm1.cmt[0]);
            for(int b=0; b<8; ++b) {
                RunCCMM(ctx, A[b], B1_CMT, CMTVariant::AFCMT, &tm_rest, true);
            }
            total_ms += MsSince(t0);
        }
        std::cout << "Optimized AFCMT: " << (total_ms / reps) / 1000.0 << " seconds." << std::endl;
    }
    
    return 0;
}
