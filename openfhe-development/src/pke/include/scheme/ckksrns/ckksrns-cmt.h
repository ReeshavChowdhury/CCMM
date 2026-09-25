//==================================================================================
// CKL CMT, FCMT, and affine-orbit AFCMT on the OpenFHE CKKS path.
// Key-switch core, evaluation keys, and parameters are unchanged.
//==================================================================================

#ifndef LBCRYPTO_CRYPTO_CKKSRNS_CMT_H
#define LBCRYPTO_CRYPTO_CKKSRNS_CMT_H

#include "openfhe.h"

#include <cstdint>
#include <string>
#include <vector>

namespace lbcrypto {

enum class CMTVariant {
    CKL_REFERENCE   = 0,
    OPENFHE_HOISTED = 1,
    SPECTRAL_TWEAK  = 2,
    FCMT            = 3,
    AFCMT           = 4
};

enum class AFCMTMode {
    AFCMT_GENERIC      = 0,
    AFCMT_RESIDUE      = 1,
    AFCMT_FUSED_RADIX2 = 2,
    AFCMT_BATCHED      = 3  // class-batched radix-2 (d/2 twiddle classes)
};

enum class StorageMode { COLUMN_STORAGE = 0, ROW_STORAGE = 1, DUAL_STORAGE = 2 };

inline const char* CMTVariantName(CMTVariant v) {
    switch (v) {
        case CMTVariant::CKL_REFERENCE:
            return "CKL_REFERENCE";
        case CMTVariant::OPENFHE_HOISTED:
            return "OPENFHE_HOISTED";
        case CMTVariant::SPECTRAL_TWEAK:
            return "SPECTRAL_TWEAK";
        case CMTVariant::FCMT:
            return "FCMT";
        case CMTVariant::AFCMT:
            return "AFCMT";
        default:
            return "UNKNOWN";
    }
}

inline const char* AFCMTModeName(AFCMTMode m) {
    switch (m) {
        case AFCMTMode::AFCMT_GENERIC:
            return "AFCMT_GENERIC";
        case AFCMTMode::AFCMT_RESIDUE:
            return "AFCMT_RESIDUE";
        case AFCMTMode::AFCMT_FUSED_RADIX2:
            return "AFCMT_FUSED_RADIX2";
        case AFCMTMode::AFCMT_BATCHED:
            return "AFCMT_BATCHED";
        default:
            return "UNKNOWN";
    }
}

inline bool CKLConditionSatisfied(uint32_t N, uint32_t d) {
    if (d == 0 || N == 0 || (N % d) != 0)
        return false;
    const uint64_t d2   = static_cast<uint64_t>(d) * static_cast<uint64_t>(d);
    const uint64_t twoN = 2ull * static_cast<uint64_t>(N);
    return (twoN % d2) == 0;
}

// Exact π_t(ν) = ν + k t (2ν+1)  (mod N), natural (non-bit-reversed) index.
uint32_t PiNatural(uint32_t nu, uint32_t t, uint32_t N, uint32_t k);

// General affine: b' = ((1+2kt)b + t(2a+1)) mod d.
uint32_t LaneAffineGeneral(uint32_t a, uint32_t b, uint32_t t, uint32_t k, uint32_t d);
// CKL-simplified (requires d | 2k): b' = b + t(2a+1) mod d.
uint32_t LaneAffineSimplified(uint32_t a, uint32_t b, uint32_t t, uint32_t d);
// Uses simplified formula; asserts CKL admissibility.
uint32_t LaneAffine(uint32_t a, uint32_t b, uint32_t t, uint32_t k, uint32_t d);
uint32_t PiAffine(uint32_t nu, uint32_t t, uint32_t N, uint32_t k, uint32_t d);

struct ClassCountReport {
    uint32_t unique_gamma{0};
    uint32_t expected_gamma{0};
    uint32_t min_gamma_count{0};
    uint32_t max_gamma_count{0};
    uint32_t expected_gamma_count{0};  // 2k over N frequencies
    uint32_t unique_step{0};
    uint32_t expected_step{0};
    uint32_t min_step_count{0};
    uint32_t max_step_count{0};
    uint32_t expected_step_count{0};  // 2k/d over a in [0,k)
    bool gamma_ok{false};
    bool step_ok{false};
};

ClassCountReport CountTransformClasses(uint32_t N, uint32_t d);

struct IndexKernelStats {
    double ns_per_index{0};
    double indices_per_sec{0};
    uint64_t checksum{0};
};

// Experiment A: no CKKS arithmetic.
IndexKernelStats BenchIndexDirect(uint32_t N, uint32_t d, uint32_t iters);
IndexKernelStats BenchIndexLookup(uint32_t N, uint32_t d, uint32_t iters);
IndexKernelStats BenchIndexRecurrence(uint32_t N, uint32_t d, uint32_t iters);

struct CMTIndexTable {
    uint32_t N{0};
    uint32_t d{0};
    uint32_t k{0};
    uint32_t logN{0};
    std::vector<uint32_t> h;
    std::vector<uint32_t> tstar;
    std::vector<uint32_t> hstar;
    std::vector<uint32_t> bitrev;              // natural ν -> OpenFHE storage
    std::vector<uint32_t> alpha;               // α_t mod d
    std::vector<std::vector<uint32_t>> beta;   // β[a][t]
    std::vector<std::vector<std::vector<uint32_t>>> lane_map;  // [a][t][b]
    std::vector<std::vector<uint32_t>> pi_nat; // [t][ν] natural π
    std::vector<std::vector<uint32_t>> pi_sto; // [t][ν_nat] storage source index
    std::vector<uint32_t> dft_rev;             // bit-reverse of 0..d-1
    std::vector<uint32_t> h_for;               // stream ts -> Galois h
};

struct CMTTraffic {
    uint64_t bytes_written{0};
    uint64_t bytes_read{0};
    uint64_t gather_bytes{0};
    uint64_t auto_materialize_bytes{0};
    uint64_t residue_layout_bytes{0};
    uint64_t automorphism_output_bytes_written{0};
    uint64_t automorphism_output_bytes_read{0};
};

struct CMTTimings {
    double t_forward_tweak_ms{0};
    double t_hoist_ms{0};
    double t_keyprod_ms{0};  // EvalFastKeySwitchCore only
    double t_automorphism_ms{0};
    double t_spectral_tweak_ms{0};
    double t_gather_ms{0};
    double t_affine_index_ms{0};
    double t_layout_ms{0};
    double t_lane_permutation_ms{0};
    double t_dft_ms{0};
    double t_final_phase_ms{0};
    double t_twiddle_load_ms{0};
    double t_tweak_fwd_ms{0};
    double t_tweak_inv_ms{0};
    double t_monomial_ms{0};
    double t_cmt_ms{0};
    CMTTraffic traffic;
};

struct CMTContext {
    CryptoContext<DCRTPoly> cc;
    PublicKey<DCRTPoly> pk;
    PrivateKey<DCRTPoly> sk;
    CMTIndexTable idx;
    AFCMTMode afcmtMode{AFCMTMode::AFCMT_FUSED_RADIX2};
    int64_t encodeScale{1LL << 40};
    // Parameter-only monomials, built once at context setup (not on the timed path).
    std::vector<DCRTPoly> evalXi;
    DCRTPoly evalX2k;   // X^{+2k}  (forward TWEAK)
    DCRTPoly monoM2k;   // X^{-2k}  (fused consume γ)
    DCRTPoly monoM1;
    bool haveMonos{false};
    // Setup-only DFT stage roots: [tow][sto * logd + stage] = omega_sto^{d/len}.
    // Used by FastEvalTweak / SpectralConsume so the timed path does no ModExp.
    uint32_t logd{0};
    std::vector<std::vector<NativeInteger>> wlenFwd;
    std::vector<std::vector<NativeInteger>> wlenInv;
    bool haveWlen{false};
};

CMTIndexTable BuildCMTIndexTable(uint32_t N, uint32_t d);
std::vector<uint32_t> CMTGaloisIndices(const CMTIndexTable& idx);
void GenerateCMTEvalKeys(const CMTContext& ctx);
CMTContext MakeCMTContext(uint32_t N, uint32_t d, uint32_t multDepth = 2, uint32_t scaleModSize = 50);

// Named CKKS sets from batchmm_sec (S12 / S13b / S13H / F16) or HE standard 128-bit.
struct CKKSParamSet {
    const char* name{""};
    uint32_t logN{0};           // 0 = let OpenFHE choose (STD128)
    uint32_t scaleBits{50};     // log2 Δ
    uint32_t firstModBits{60};  // first Q limb
    uint32_t depth{2};          // extra scale limbs (multiplicative depth)
    uint32_t dnum{0};           // HYBRID digits; 0 = OpenFHE default
    bool std128{false};         // HEStd_128_classic
    bool sparseHwt{false};      // SPARSE_TERNARY
};
CMTContext MakeCMTContextForSet(uint32_t d, const CKKSParamSet& ps);

uint64_t EvalFormCiphertextBytes(const Ciphertext<DCRTPoly>& ct);

std::vector<Ciphertext<DCRTPoly>> EncryptMatrixColumns(const CMTContext& ctx,
                                                       const std::vector<std::vector<int64_t>>& M);
std::vector<Ciphertext<DCRTPoly>> EncryptMatrixRows(const CMTContext& ctx,
                                                    const std::vector<std::vector<int64_t>>& M);
std::vector<std::vector<double>> DecryptMatrixColumns(const CMTContext& ctx,
                                                      const std::vector<Ciphertext<DCRTPoly>>& cts);
std::vector<std::vector<int64_t>> CleartextCMT(const std::vector<std::vector<int64_t>>& M);

std::vector<Ciphertext<DCRTPoly>> RunCMT(const CMTContext& ctx, const std::vector<Ciphertext<DCRTPoly>>& cols,
                                         CMTVariant variant, CMTTimings* timings = nullptr);

struct CCMMTimings {
    double t_ccmm_ms{0};  // exclusive wall-clock of one online CCMM
    double t_cmt_ms[3]{};
    double t_ppmm_ms[4]{};
    double t_layout_ms{0};      // ModuleTranspose + PackCipher
    double t_relin_ms{0};
    double t_accumulate_ms{0};
    double t_cmt_sum_ms{0};     // exclusive sum of the 3 CMT scopes in this rep
    double t_ppmm_sum_ms{0};    // exclusive sum of the 4 PPMM scopes in this rep
    uint32_t n_cmt{0};
    uint32_t n_ppmm{0};
    uint32_t n_relin{0};
    CMTTimings cmt[3];
    // Legacy aliases used by older harnesses.
    double t_ppmm_ms_legacy{0};
    double t_cmt_total_ms{0};
};

// Square CKL CCMM: 3 CMT + 4 PPMM + relin + accumulate.
// skipFirstCMT: right operand is already row-oriented. CMT_1 timed as 0.
// PPMM / module-transpose / relin are the same OpenMP fused kernels for every
// variant. Only the CMT backend (CKL vs AFCMT) changes.
// fastKernels is accepted for ABI compatibility and is ignored (always on).
std::vector<Ciphertext<DCRTPoly>> RunCCMM(const CMTContext& ctx, const std::vector<Ciphertext<DCRTPoly>>& leftCols,
                                          const std::vector<Ciphertext<DCRTPoly>>& rightCols, CMTVariant variant,
                                          CCMMTimings* timings = nullptr, bool skipFirstCMT = false,
                                          bool fastKernels = true);

std::vector<Ciphertext<DCRTPoly>> RunRectangularCCMM(
    const CMTContext& ctx,
    const std::vector<Ciphertext<DCRTPoly>>& leftCols,
    const std::vector<Ciphertext<DCRTPoly>>& rightCols,
    uint32_t B,
    CMTVariant variant,
    CCMMTimings* tm);


// Rectangular d x 2d  times  2d x d, as two square blocks.
std::vector<Ciphertext<DCRTPoly>> RunCCMMRectangular(const CMTContext& ctx,
                                                     const std::vector<std::vector<int64_t>>& A,  // d x 2d
                                                     const std::vector<std::vector<int64_t>>& B,  // 2d x d
                                                     CMTVariant variant, CCMMTimings* timings = nullptr);

struct MathTestReport {
    bool pass{true};
    std::string name;
    std::string detail;
};

std::vector<MathTestReport> RunMathUnitTests(uint32_t N, uint32_t d, CMTContext* ctxOrNull);

}  // namespace lbcrypto

#endif
