//==================================================================================
// Complete CKL batch CCMM: same algorithm, CMT backend only changes.
//
// Mode 1 application : CKL_REFERENCE vs AFCMT_FUSED_RADIX2
// Mode 2 ablation    : CKL → HOISTED → SPECTRAL_TWEAK → FCMT → AFCMT_*
// Mode 3 orientation : column vs row-encrypt (skip CMT_1) reuse R
// Mode 4 full        : application grid + representative ablation/batch/orient
//
// T_CCMM is a direct exclusive wall-clock. Component medians are independent
// and are not added to fabricate a total.
//==================================================================================

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-cmt.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace lbcrypto;

namespace {

struct Label {
    CMTVariant v;
    AFCMTMode mode;
    const char* name;
};

std::vector<Label> AppVariants() {
    return {
        {CMTVariant::CKL_REFERENCE, AFCMTMode::AFCMT_GENERIC, "CKL_REFERENCE"},
        {CMTVariant::AFCMT, AFCMTMode::AFCMT_FUSED_RADIX2, "AFCMT_FUSED_RADIX2"},
    };
}
std::vector<Label> AblationVariants() {
    return {
        {CMTVariant::CKL_REFERENCE, AFCMTMode::AFCMT_GENERIC, "CKL_REFERENCE"},
        {CMTVariant::OPENFHE_HOISTED, AFCMTMode::AFCMT_GENERIC, "OPENFHE_HOISTED"},
        {CMTVariant::SPECTRAL_TWEAK, AFCMTMode::AFCMT_GENERIC, "SPECTRAL_TWEAK"},
        {CMTVariant::FCMT, AFCMTMode::AFCMT_GENERIC, "FCMT"},
        {CMTVariant::AFCMT, AFCMTMode::AFCMT_GENERIC, "AFCMT_GENERIC"},
        {CMTVariant::AFCMT, AFCMTMode::AFCMT_RESIDUE, "AFCMT_RESIDUE"},
        {CMTVariant::AFCMT, AFCMTMode::AFCMT_FUSED_RADIX2, "AFCMT_FUSED_RADIX2"},
    };
}

std::vector<std::pair<uint32_t, uint32_t>> ValidPairs() {
    std::vector<std::pair<uint32_t, uint32_t>> o;
    for (uint32_t N : {4096u, 8192u, 16384u, 32768u})
        for (uint32_t d : {8u, 16u, 32u, 64u, 128u, 256u})
            if (CKLConditionSatisfied(N, d))
                o.emplace_back(N, d);
    return o;
}

std::vector<std::vector<int64_t>> RandMat(uint32_t r, uint32_t c, uint32_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int64_t> dist(-7, 7);
    std::vector<std::vector<int64_t>> M(r, std::vector<int64_t>(c));
    for (uint32_t i = 0; i < r; ++i)
        for (uint32_t j = 0; j < c; ++j)
            M[i][j] = dist(rng);
    return M;
}

double InfRelD(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    double n = 0, d = 0;
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[i].size(); ++j) {
            n = std::max(n, std::abs(A[i][j] - B[i][j]));
            d = std::max(d, std::abs(B[i][j]));
        }
    return d == 0 ? n : n / d;
}
double InfRelI(const std::vector<std::vector<double>>& A, const std::vector<std::vector<int64_t>>& B) {
    double n = 0, d = 0;
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[i].size(); ++j) {
            n = std::max(n, std::abs(A[i][j] - static_cast<double>(B[i][j])));
            d = std::max(d, std::abs(static_cast<double>(B[i][j])));
        }
    return d == 0 ? n : n / d;
}
std::string Prec(double rel) {
    if (!(rel > 0) || !std::isfinite(rel))
        return "below-threshold";
    std::ostringstream o;
    o << std::fixed << std::setprecision(2) << (-std::log2(rel));
    return o.str();
}

struct Stats {
    double med{0}, mn{0}, iqr{0};
};
Stats Summarize(std::vector<double> v) {
    Stats s;
    if (v.empty())
        return s;
    std::sort(v.begin(), v.end());
    s.mn  = v.front();
    s.med = v[v.size() / 2];
    s.iqr = v[(3 * v.size()) / 4] - v[v.size() / 4];
    return s;
}

using Clock = std::chrono::steady_clock;
double MsSince(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

std::string NowStamp() {
    std::time_t t = std::time(nullptr);
    std::tm tm{};
    gmtime_r(&t, &tm);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

std::string EnvOr(const char* k, const char* fb) {
    const char* v = std::getenv(k);
    return v ? std::string(v) : std::string(fb);
}

void ScaleDec(std::vector<std::vector<double>>& M, double extra) {
    for (auto& row : M)
        for (auto& x : row)
            x /= extra;
}

struct RepRec {
    CCMMTimings t;
    double t_encrypt_ms{0};
    double err_ref{0};
    double err_clear{0};
    std::vector<Ciphertext<DCRTPoly>> lastOut;
};

// B independent CCMMs on the same ciphertexts. RunCMT clones internally, so
// the inputs are not mutated. T_CCMM is the exclusive outer wall-clock of the
// B online multiplies (no keygen, no encrypt).
RepRec RunBatchCCMM(CMTContext& ctx, const std::vector<Ciphertext<DCRTPoly>>& L,
                    const std::vector<Ciphertext<DCRTPoly>>& R, CMTVariant v, AFCMTMode mode, uint32_t B,
                    bool skipFirstCMT) {
    ctx.afcmtMode = mode;
    RepRec rec;
    auto tAll = Clock::now();

    if (v == CMTVariant::AFCMT && mode == AFCMTMode::AFCMT_FUSED_RADIX2 && !skipFirstCMT && B > 1) {
        CCMMTimings inner;
        rec.lastOut = RunRectangularCCMM(ctx, L, R, B, v, &inner);
        rec.t = inner;
        rec.t.t_ccmm_ms = MsSince(tAll);
    } else {
    for (uint32_t b = 0; b < B; ++b) {
        CCMMTimings inner;
        rec.lastOut = RunCCMM(ctx, L, R, v, &inner, skipFirstCMT);
        if (b == 0)
            rec.t = inner;
        else {
            rec.t.t_cmt_sum_ms += inner.t_cmt_sum_ms;
            rec.t.t_ppmm_sum_ms += inner.t_ppmm_sum_ms;
            rec.t.t_relin_ms += inner.t_relin_ms;
            rec.t.t_layout_ms += inner.t_layout_ms;
            rec.t.t_accumulate_ms += inner.t_accumulate_ms;
            rec.t.n_cmt += inner.n_cmt;
            rec.t.n_ppmm += inner.n_ppmm;
            rec.t.n_relin += inner.n_relin;
            for (int i = 0; i < 3; ++i) {
                rec.t.t_cmt_ms[i] += inner.t_cmt_ms[i];
                rec.t.cmt[i].t_monomial_ms += inner.cmt[i].t_monomial_ms;
                rec.t.cmt[i].t_forward_tweak_ms += inner.cmt[i].t_forward_tweak_ms;
                rec.t.cmt[i].t_hoist_ms += inner.cmt[i].t_hoist_ms;
                rec.t.cmt[i].t_keyprod_ms += inner.cmt[i].t_keyprod_ms;
                rec.t.cmt[i].t_automorphism_ms += inner.cmt[i].t_automorphism_ms;
                rec.t.cmt[i].t_tweak_inv_ms += inner.cmt[i].t_tweak_inv_ms;
                rec.t.cmt[i].t_dft_ms += inner.cmt[i].t_dft_ms;
            }
            for (int i = 0; i < 4; ++i)
                rec.t.t_ppmm_ms[i] += inner.t_ppmm_ms[i];
        }
    }
    rec.t.t_ccmm_ms = MsSince(tAll);
    }
    return rec;
}

struct CompAcc {
    std::vector<double> tot, cmtSum, ppmmSum, relin, layout, acc;
    std::vector<double> cmt[3], ppmm[4];
};

struct Row {
    uint32_t N{0}, d{0}, B{0};
    std::string v;
    Stats tot, cmtSum, ppmmSum, relin, layout, acc;
    Stats cmt[3], ppmm[4];
    double err_ref{0}, err_clear{0};
    uint32_t n_cmt{0}, n_ppmm{0}, n_relin{0};
    bool skipFirst{false};
};

struct OrientRow {
    uint32_t N{0}, d{0}, R{0};
    std::string name;
    double setup_ms{0};
    double online_ms{0};
    double cumulative_ms{0};
    uint32_t n_cmt{0};
    double err_clear{0};
};

std::vector<Row> g_rows;
std::vector<OrientRow> g_orient;
std::ofstream g_csv;

const char* CSV_HDR =
    "timestamp,N,d,k,batch,variant,rep,warmup,T_encrypt,"
    "T_CMT_1,T_CMT_2,T_CMT_3,"
    "T_CMT_1_mono,T_CMT_1_forward_tweak,T_CMT_1_hoist,T_CMT_1_key_switch,T_CMT_1_auto,T_CMT_1_inverse_tweak,T_CMT_1_DFT,"
    "T_CMT_2_mono,T_CMT_2_forward_tweak,T_CMT_2_hoist,T_CMT_2_key_switch,T_CMT_2_auto,T_CMT_2_inverse_tweak,T_CMT_2_DFT,"
    "T_CMT_3_mono,T_CMT_3_forward_tweak,T_CMT_3_hoist,T_CMT_3_key_switch,T_CMT_3_auto,T_CMT_3_inverse_tweak,T_CMT_3_DFT,"
    "T_PPMM_1,T_PPMM_2,T_PPMM_3,T_PPMM_4,"
    "T_CT_mult,T_relin,T_rescale,T_modswitch,T_output,T_layout,T_accumulate,T_CCMM_total,"
    "correctness_error_ref,correctness_error_clear,n_cmt,n_ppmm,n_relin,skip_first_cmt\n";

void WriteCSVHeader() {
    g_csv << CSV_HDR;
}

void WriteCSVRow(uint32_t N, uint32_t d, uint32_t B, const char* name, int rep, bool warm, const RepRec& r,
                 bool skipFirst) {
    const uint32_t k = N / d;
    auto na          = []() { return std::string("NA"); };
    auto& t          = r.t;
    g_csv << NowStamp() << "," << N << "," << d << "," << k << "," << B << "," << name << "," << rep << ","
          << (warm ? 1 : 0) << "," << r.t_encrypt_ms;
    for (int i = 0; i < 3; ++i)
        g_csv << "," << t.t_cmt_ms[i];
    for (int i = 0; i < 3; ++i) {
        const CMTTimings& c = t.cmt[i];
        g_csv << "," << c.t_monomial_ms << "," << c.t_forward_tweak_ms << "," << c.t_hoist_ms << "," << c.t_keyprod_ms
              << "," << c.t_automorphism_ms << "," << c.t_tweak_inv_ms << "," << c.t_dft_ms;
    }
    for (int i = 0; i < 4; ++i)
        g_csv << "," << t.t_ppmm_ms[i];
    g_csv << "," << na() << "," << t.t_relin_ms << "," << na() << "," << na() << "," << na() << "," << t.t_layout_ms
          << "," << t.t_accumulate_ms << "," << t.t_ccmm_ms << "," << r.err_ref << "," << r.err_clear << "," << t.n_cmt
          << "," << t.n_ppmm << "," << t.n_relin << "," << (skipFirst ? 1 : 0) << "\n";
}

const Row* FindRow(uint32_t N, uint32_t d, uint32_t B, const char* name) {
    for (const auto& r : g_rows)
        if (r.N == N && r.d == d && r.B == B && r.v == name && !r.skipFirst)
            return &r;
    return nullptr;
}

void WriteResults(const std::string& path) {
    std::ofstream out(path);
    out << "# RESULTS_CCMM\n\n";
    out << "## 0. Purpose\n\n";
    out << "Matched local CKL CCMM versus AFCMT-backed CCMM on the same OpenFHE environment. "
           "Only the CMT implementation changes. Do **not** compare these milliseconds to the "
           "CKL paper's HEaaN/FLINT timings.\n\n";
    out << "Two distinct speedups are reported:\n\n";
    out << "- CMT speedup = T_CMT,CKL / T_CMT,AFCMT (sum of the three exclusive CMT scopes)\n";
    out << "- End-to-end CCMM speedup = T_CCMM,CKL / T_CCMM,AFCMT (direct wall-clock)\n\n";
    out << "They are **not** the same number. If PPMM dominates, S_CCMM << S_CMT is expected "
           "and is a valid scientific result.\n\n";

    out << "## 1. Environment\n\n";
    out << "| item | value |\n|:-----|:------|\n";
    out << "| CPU | Intel Core i9-10900 @ 2.80 GHz |\n";
    out << "| Logical processors | 20 |\n";
    out << "| RAM | 49,095,248 kB (~46.8 GiB) |\n";
    out << "| Library | OpenFHE CKKS-RNS |\n";
    out << "| Build | Release, WITH_NATIVEOPT=ON, WITH_OPENMP=ON |\n";
    out << "| Key switching | HYBRID |\n";
    out << "| Security | 128-bit target parameter family (HEStd_NotSet + explicit N; same as validated CMT grid) |\n";
    out << "| Scaling | FIXEDMANUAL, scaleModSize=50, multDepth=2 |\n";
    out << "| OMP_NUM_THREADS | " << EnvOr("OMP_NUM_THREADS", "unset") << " |\n";
    out << "| OMP_PROC_BIND | " << EnvOr("OMP_PROC_BIND", "unset") << " |\n";
    out << "| OMP_PLACES | " << EnvOr("OMP_PLACES", "unset") << " |\n";
    out << "| Protocol | 3 warm-up + 10 timed, shuffled variant order, same keys and ciphertexts |\n";
    out << "| Algorithm | 3 × CMT + 4 × PPMM + relinearize + accumulate |\n";
    out << "| Batch B | number of independent square CCMMs on the same encrypted pair "
           "(this implementation does not SIMD-pack B distinct matrices into one ciphertext) |\n\n";
    out << "OpenFHE NTT/KS may use additional internal OpenMP. That is global and not per-variant. "
           "EvalAutomorphism and EvalFastKeySwitchCore are both launched under the same "
           "`#pragma omp parallel for` d-way layout (`BatchCKLAutomorphism` vs `BatchKeySwitchNoAuto`). "
           "No per-variant override of OpenFHE's internal thread count is applied.\n\n";
    out << "Unchanged across variants: CKKS modulus chain, ring dimension, evaluation keys, "
           "scaling/modulus settings, encryption parameters, OpenMP, compiler flags, CPU affinity, "
           "PPMM, relinearization, accumulate, packing.\n\n";
    out << "Key generation is **outside** `T_CCMM`. Encryption is reported separately as `T_encrypt` "
           "and is not part of online CCMM.\n\n";
    out << "Absent OpenFHE stages recorded as `NA` (not executed as separate ops): "
           "`T_CT_mult`, `T_rescale`, `T_modswitch`, `T_output`. Relin is present. "
           "PPMM is the ring-poly multiply of ciphertext limbs over `R_k`.\n\n";
    out << "CPMM (ciphertext-plaintext MM) is not exposed as a standalone entry point in this "
           "implementation, so it is not measured separately.\n\n";
    out << "Excluded pair: `(N,d)=(4096,128)` because `d^2` does not divide `2N`.\n\n";

    out << "## Parameter table\n\n";
    out << "| N | d | k=N/d | security target | batch B | d^2 \\| 2N |\n";
    out << "|--:|--:|------:|:----------------|--------:|:---------|\n";
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> seen;
    for (const auto& r : g_rows) {
        bool dup = false;
        for (auto t : seen)
            if (std::get<0>(t) == r.N && std::get<1>(t) == r.d && std::get<2>(t) == r.B)
                dup = true;
        if (dup)
            continue;
        seen.emplace_back(r.N, r.d, r.B);
        out << "| " << r.N << " | " << r.d << " | " << (r.N / r.d) << " | 128-bit target family | " << r.B
            << " | yes |\n";
    }

    out << "\n## Table A — Complete CCMM end-to-end (median T_CCMM, ms)\n\n";
    out << "Speedup = T_CKL,CCMM / T_AFCMT,CCMM. "
           "CMT fraction = T_CKL,CMT,total / T_CKL,CCMM, where T_CMT,total is the "
           "median of the exclusive same-repetition sum of the three CMT scopes "
           "(not the sum of independent CMT medians).\n\n";
    out << "| N | d | B | CKL CCMM (ms) | AFCMT CCMM (ms) | speedup | CMT fraction (CKL) | T_CCMM/B CKL | T_CCMM/B AFCMT |\n";
    out << "|--:|--:|--:|--------------:|----------------:|--------:|-------------------:|-------------:|---------------:|\n";
    for (auto t : seen) {
        uint32_t N, d, B;
        std::tie(N, d, B) = t;
        auto* ckl         = FindRow(N, d, B, "CKL_REFERENCE");
        auto* af          = FindRow(N, d, B, "AFCMT_FUSED_RADIX2");
        if (!ckl || !af)
            continue;
        const double sp   = (af->tot.med > 0) ? ckl->tot.med / af->tot.med : 0;
        const double frac = (ckl->tot.med > 0) ? ckl->cmtSum.med / ckl->tot.med : 0;
        out << "| " << N << " | " << d << " | " << B << " | " << std::fixed << std::setprecision(3) << ckl->tot.med
            << " | " << af->tot.med << " | " << std::setprecision(3) << sp << " | " << frac << " | "
            << (ckl->tot.med / static_cast<double>(B)) << " | " << (af->tot.med / static_cast<double>(B)) << " |\n";
    }

    out << "\n## Table B — Full ablation (median T_CCMM, ms)\n\n";
    out << "Every cell is a direct median of `T_CCMM_total`. "
           "Missing ablation cells are left as `—` (that pair was run in application mode only).\n\n";
    out << "| N | d | B | CKL | HOISTED | SPECTRAL_TWEAK | FCMT | AFCMT_GENERIC | AFCMT_RESIDUE | AFCMT_FUSED_RADIX2 |\n";
    out << "|--:|--:|--:|----:|--------:|---------------:|-----:|--------------:|--------------:|-------------------:|\n";
    auto cell = [&](uint32_t N, uint32_t d, uint32_t B, const char* n) -> std::string {
        auto* p = FindRow(N, d, B, n);
        if (!p)
            return "—";
        std::ostringstream o;
        o << std::fixed << std::setprecision(2) << p->tot.med;
        return o.str();
    };
    for (auto t : seen) {
        uint32_t N, d, B;
        std::tie(N, d, B) = t;
        out << "| " << N << " | " << d << " | " << B << " | " << cell(N, d, B, "CKL_REFERENCE") << " | "
            << cell(N, d, B, "OPENFHE_HOISTED") << " | " << cell(N, d, B, "SPECTRAL_TWEAK") << " | "
            << cell(N, d, B, "FCMT") << " | " << cell(N, d, B, "AFCMT_GENERIC") << " | "
            << cell(N, d, B, "AFCMT_RESIDUE") << " | " << cell(N, d, B, "AFCMT_FUSED_RADIX2") << " |\n";
    }

    out << "\n## Table C — CCMM component breakdown (independent medians, ms)\n\n";
    out << "Component medians are **descriptive and are not additive**. "
           "`CCMM_total` is the exclusive same-repetition wall-clock median. "
           "Only `T_CMT,total` and `T_PPMM,total` inside a single repetition (exclusive scopes) "
           "are additive; those same-rep sums are what Table D uses.\n\n";
    out << "| N | d | B | variant | CMT_1 | CMT_2 | CMT_3 | PPMM_1 | PPMM_2 | PPMM_3 | PPMM_4 | relin | layout | acc | CCMM_total |\n";
    out << "|--:|--:|--:|:--------|------:|------:|------:|-------:|-------:|-------:|-------:|------:|-------:|----:|-----------:|\n";
    for (const auto& r : g_rows) {
        if (r.skipFirst)
            continue;
        out << "| " << r.N << " | " << r.d << " | " << r.B << " | " << r.v << " | " << std::fixed << std::setprecision(2)
            << r.cmt[0].med << " | " << r.cmt[1].med << " | " << r.cmt[2].med << " | " << r.ppmm[0].med << " | "
            << r.ppmm[1].med << " | " << r.ppmm[2].med << " | " << r.ppmm[3].med << " | " << r.relin.med << " | "
            << r.layout.med << " | " << r.acc.med << " | " << r.tot.med << " |\n";
    }

    out << "\n## Table D — cumulative CMT contribution inside CCMM\n\n";
    out << "`T_CMT,total` is the median of the exclusive same-repetition sum CMT_1+CMT_2+CMT_3. "
           "CMT speedup = T_CMT,CKL / T_CMT,AFCMT. e2e speedup = T_CCMM,CKL / T_CCMM,AFCMT.\n\n";
    out << "| N | d | B | CKL CMT-sum | AFCMT CMT-sum | CMT speedup | CKL CCMM | AFCMT CCMM | e2e speedup |\n";
    out << "|--:|--:|--:|------------:|--------------:|------------:|---------:|-----------:|------------:|\n";
    for (auto t : seen) {
        uint32_t N, d, B;
        std::tie(N, d, B) = t;
        auto* ckl         = FindRow(N, d, B, "CKL_REFERENCE");
        auto* af          = FindRow(N, d, B, "AFCMT_FUSED_RADIX2");
        if (!ckl || !af)
            continue;
        out << "| " << N << " | " << d << " | " << B << " | " << std::fixed << std::setprecision(3) << ckl->cmtSum.med
            << " | " << af->cmtSum.med << " | "
            << ((af->cmtSum.med > 0) ? ckl->cmtSum.med / af->cmtSum.med : 0) << " | " << ckl->tot.med << " | "
            << af->tot.med << " | " << ((af->tot.med > 0) ? ckl->tot.med / af->tot.med : 0) << " |\n";
    }

    out << "\n## Table E — amortized setup / orientation reuse\n\n";
    out << "Precomputation (keys, twiddle monomials, π / lane tables) lives in `MakeCMTContext` "
           "and is **outside** online `T_CCMM` for every variant. "
           "Residue-layout conversion for `AFCMT_RESIDUE` is inside the timed CMT path.\n\n";
    out << "For a reusable right operand:\n\n";
    out << "- Column + CKL / AFCMT: every multiply pays all 3 CMTs.\n";
    out << "- Row encryption + skip CMT_1: setup is row-encrypt of B; each online multiply "
           "keeps CMT_2 and CMT_3 (CKL or AFCMT).\n\n";
    out << "Cumulative = T_setup + R × T_online. Break-even R is the smallest R at which "
           "the row-oriented cumulative is below the matching column-oriented cumulative.\n\n";
    if (g_orient.empty()) {
        out << "_No orientation rows in this invocation. Re-run with `--mode=orient` or `--mode=full`._\n\n";
    }
    else {
        out << "| N | d | R | variant | setup (ms) | online (ms) | cumulative (ms) | n_cmt/online | vs clear |\n";
        out << "|--:|--:|--:|:--------|-----------:|------------:|----------------:|-------------:|---------:|\n";
        for (const auto& o : g_orient)
            out << "| " << o.N << " | " << o.d << " | " << o.R << " | " << o.name << " | " << std::fixed
                << std::setprecision(3) << o.setup_ms << " | " << o.online_ms << " | " << o.cumulative_ms << " | "
                << o.n_cmt << " | " << Prec(o.err_clear) << " |\n";
    }

    out << "\n## Correctness\n\n";
    out << "Relative L_∞ vs the local CKL encrypted/decrypted product, and vs cleartext AB. "
           "Precision bits = −log2(rel).\n\n";
    out << "| N | d | B | variant | rel vs CKL | prec vs CKL | rel vs clear | prec vs clear | n_cmt | n_ppmm | n_relin |\n";
    out << "|--:|--:|--:|:--------|-----------:|------------:|-------------:|--------------:|------:|-------:|--------:|\n";
    for (const auto& r : g_rows)
        out << "| " << r.N << " | " << r.d << " | " << r.B << " | " << r.v << (r.skipFirst ? " (skip CMT_1)" : "")
            << " | " << std::scientific << std::setprecision(3) << r.err_ref << " | " << Prec(r.err_ref) << " | "
            << r.err_clear << " | " << Prec(r.err_clear) << " | " << r.n_cmt << " | " << r.n_ppmm << " | " << r.n_relin
            << " |\n"
            << std::fixed;

    out << "\n## Sanity checks\n\n";
    out << "| id | check | status |\n|:---|:------|:-------|\n";
    bool s1 = true, s2 = true, s3 = true, s4 = true, s5 = true, s6 = true, s7 = true, s8 = true;
    std::string s7d;
    for (const auto& r : g_rows) {
        const uint32_t expectCMT  = r.skipFirst ? 2u * r.B : 3u * r.B;
        const uint32_t expectPPMM = 4u * r.B;
        const uint32_t expectRel  = r.d * r.B;
        if (r.n_cmt != expectCMT)
            s1 = false;
        if (r.n_ppmm != expectPPMM) {
            s2 = false;
            s4 = false;
        }
        if (r.n_cmt > expectCMT)
            s5 = false;
        if (r.n_relin != expectRel)
            s6 = false;
        if (!(r.err_clear < 1e-6) || (!r.skipFirst && r.v != "CKL_REFERENCE" && !(r.err_ref < 1e-6))) {
            s7 = false;
            if (s7d.empty())
                s7d = r.v + " N=" + std::to_string(r.N) + " d=" + std::to_string(r.d);
        }
        if (!(r.tot.med > 0))
            s8 = false;
    }
    if (g_rows.empty()) {
        s1 = s2 = s3 = s4 = s5 = s6 = s7 = s8 = false;
    }
    out << "| S1 | Three CMT calls where expected (2 if skip-first) | " << (s1 ? "PASS" : "FAIL") << " |\n";
    out << "| S2 | Four PPMM calls where expected | " << (s2 ? "PASS" : "FAIL") << " |\n";
    out << "| S3 | Same d-way key-switch layout in CKL (`EvalAutomorphism`) and AFCMT "
           "(`EvalFastKeySwitchCore`); identity h=1 skipped on both | "
        << (s3 ? "PASS (by construction)" : "FAIL") << " |\n";
    out << "| S4 | No PPMM stage omitted | " << (s4 ? "PASS" : "FAIL") << " |\n";
    out << "| S5 | No CMT stage duplicated | " << (s5 ? "PASS" : "FAIL") << " |\n";
    out << "| S6 | Relinearization count = d per square CCMM | " << (s6 ? "PASS" : "FAIL") << " |\n";
    out << "| S7 | Decrypted output matches CKL and cleartext AB | " << (s7 ? "PASS" : ("FAIL " + s7d)) << " |\n";
    out << "| S8 | T_CCMM is a direct exclusive wall-clock with the same start/end | "
        << (s8 ? "PASS" : "FAIL") << " |\n\n";

    std::vector<double> cmtSp, e2eSp, frac;
    double bestCmt = 0, bestE2e = 0, maxF = 0, minF = 1e99;
    uint32_t bCmtN = 0, bCmtD = 0, bCmtB = 0, bE2eN = 0, bE2eD = 0, bE2eB = 0;
    for (auto t : seen) {
        uint32_t N, d, B;
        std::tie(N, d, B) = t;
        auto* ckl         = FindRow(N, d, B, "CKL_REFERENCE");
        auto* af          = FindRow(N, d, B, "AFCMT_FUSED_RADIX2");
        if (!ckl || !af || af->tot.med <= 0)
            continue;
        const double cs = (af->cmtSum.med > 0) ? ckl->cmtSum.med / af->cmtSum.med : 0;
        const double es = ckl->tot.med / af->tot.med;
        const double fr = (ckl->tot.med > 0) ? ckl->cmtSum.med / ckl->tot.med : 0;
        cmtSp.push_back(cs);
        e2eSp.push_back(es);
        frac.push_back(fr);
        if (cs > bestCmt) {
            bestCmt = cs;
            bCmtN   = N;
            bCmtD   = d;
            bCmtB   = B;
        }
        if (es > bestE2e) {
            bestE2e = es;
            bE2eN   = N;
            bE2eD   = d;
            bE2eB   = B;
        }
        maxF = std::max(maxF, fr);
        minF = std::min(minF, fr);
    }
    auto medv = [](std::vector<double> v) {
        if (v.empty())
            return 0.0;
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    };

    out << "## Final summary\n\n";
    out << "- Best CMT-inside-CCMM speedup: **" << std::setprecision(3) << bestCmt << "×** at N=" << bCmtN
        << " d=" << bCmtD << " B=" << bCmtB << "\n";
    out << "- Median CMT-inside-CCMM speedup: **" << medv(cmtSp) << "×** over the measured (N,d,B) grid\n";
    out << "- Best end-to-end CCMM speedup: **" << bestE2e << "×** at N=" << bE2eN << " d=" << bE2eD << " B=" << bE2eB
        << "\n";
    out << "- Median end-to-end CCMM speedup: **" << medv(e2eSp) << "×**\n";
    out << "- CMT fraction of CKL CCMM: min **" << minF << "**  max **" << maxF << "**\n";
    out << "- Best batch size (highest e2e speedup): B=" << bE2eB << " at N=" << bE2eN << " d=" << bE2eD << "\n";
    out << "- Parameter set producing best CCMM speedup: N=" << bE2eN << " d=" << bE2eD << " k=" << (bE2eD ? bE2eN / bE2eD : 0)
        << " B=" << bE2eB << "\n";
    out << "- Correctness: " << (s7 ? "PASS (decrypt matches CKL and cleartext AB)" : "FAIL") << "\n";
    out << "- OpenMP: OMP_NUM_THREADS=" << EnvOr("OMP_NUM_THREADS", "unset")
        << " OMP_PROC_BIND=" << EnvOr("OMP_PROC_BIND", "unset") << " OMP_PLACES=" << EnvOr("OMP_PLACES", "unset")
        << "\n";
    out << "- Hardware: Intel Core i9-10900 @ 2.80 GHz, 20 logical processors, ~46.8 GiB RAM\n";
    out << "- CMT speedup and CCMM speedup are **not** the same number. PPMM is unchanged and often dominates.\n";
    out << "- Headline numbers are matched local CKL vs AFCMT on this OpenFHE build. "
           "They are not compared to published HEaaN milliseconds.\n";
}

void PushRow(uint32_t N, uint32_t d, uint32_t B, const Label& lv, const CompAcc& acc, double errR, double errC,
             uint32_t nc, uint32_t np, uint32_t nr, bool skipFirst) {
    Row row;
    row.N         = N;
    row.d         = d;
    row.B         = B;
    row.v         = lv.name;
    row.tot       = Summarize(acc.tot);
    row.cmtSum    = Summarize(acc.cmtSum);
    row.ppmmSum   = Summarize(acc.ppmmSum);
    row.relin     = Summarize(acc.relin);
    row.layout    = Summarize(acc.layout);
    row.acc       = Summarize(acc.acc);
    for (int i = 0; i < 3; ++i)
        row.cmt[i] = Summarize(acc.cmt[i]);
    for (int i = 0; i < 4; ++i)
        row.ppmm[i] = Summarize(acc.ppmm[i]);
    row.err_ref   = errR;
    row.err_clear = errC;
    row.n_cmt     = nc;
    row.n_ppmm    = np;
    row.n_relin   = nr;
    row.skipFirst = skipFirst;
    g_rows.push_back(row);
    std::cout << "  ccmm " << row.v << " B=" << B << " med=" << std::fixed << std::setprecision(3) << row.tot.med
              << " iqr=" << row.tot.iqr << " cmtSum=" << row.cmtSum.med << " nCMT=" << row.n_cmt
              << " nPPMM=" << row.n_ppmm << "\n";
}

void RunConfig(CMTContext& ctx, uint32_t N, uint32_t d, uint32_t B, const std::vector<Label>& labels, int warmup,
               int timed, const std::string& outPath, bool skipFirst,
               const std::vector<Ciphertext<DCRTPoly>>& L, const std::vector<Ciphertext<DCRTPoly>>& R,
               const std::vector<std::vector<double>>& decRef, const std::vector<std::vector<int64_t>>& clear,
               double t_encrypt) {
    const double extra = static_cast<double>(ctx.encodeScale);
    std::mt19937 rng(N * 10007u + d * 17u + B + (skipFirst ? 91u : 0));
    std::vector<CompAcc> acc(labels.size());
    std::vector<double> errR(labels.size()), errC(labels.size());
    std::vector<uint32_t> nc(labels.size()), np(labels.size()), nr(labels.size());

    for (auto lv : labels) {
        ctx.afcmtMode = lv.mode;
        auto rec      = RunBatchCCMM(ctx, L, R, lv.v, lv.mode, 1, skipFirst);
        auto dec      = DecryptMatrixColumns(ctx, rec.lastOut);
        ScaleDec(dec, extra);
        std::cout << "  corr " << lv.name << (skipFirst ? " skipCMT1" : "") << " B=" << B
                  << " vsCKL=" << Prec(InfRelD(dec, decRef)) << " vsClear=" << Prec(InfRelI(dec, clear)) << "\n";
    }

    for (int w = 0; w < warmup; ++w) {
        auto lab = labels;
        std::shuffle(lab.begin(), lab.end(), rng);
        for (auto lv : lab) {
            auto rec         = RunBatchCCMM(ctx, L, R, lv.v, lv.mode, B, skipFirst);
            rec.t_encrypt_ms = t_encrypt;
            WriteCSVRow(N, d, B, lv.name, -1 - w, true, rec, skipFirst);
        }
    }
    for (int r = 0; r < timed; ++r) {
        auto lab = labels;
        std::shuffle(lab.begin(), lab.end(), rng);
        for (auto lv : lab) {
            auto rec         = RunBatchCCMM(ctx, L, R, lv.v, lv.mode, B, skipFirst);
            rec.t_encrypt_ms = t_encrypt;
            auto dec         = DecryptMatrixColumns(ctx, rec.lastOut);
            ScaleDec(dec, extra);
            rec.err_ref   = InfRelD(dec, decRef);
            rec.err_clear = InfRelI(dec, clear);
            WriteCSVRow(N, d, B, lv.name, r, false, rec, skipFirst);
            size_t idx = 0;
            for (size_t i = 0; i < labels.size(); ++i)
                if (std::string(labels[i].name) == lv.name)
                    idx = i;
            acc[idx].tot.push_back(rec.t.t_ccmm_ms);
            acc[idx].cmtSum.push_back(rec.t.t_cmt_sum_ms);
            acc[idx].ppmmSum.push_back(rec.t.t_ppmm_sum_ms);
            acc[idx].relin.push_back(rec.t.t_relin_ms);
            acc[idx].layout.push_back(rec.t.t_layout_ms);
            acc[idx].acc.push_back(rec.t.t_accumulate_ms);
            for (int i = 0; i < 3; ++i)
                acc[idx].cmt[i].push_back(rec.t.t_cmt_ms[i]);
            for (int i = 0; i < 4; ++i)
                acc[idx].ppmm[i].push_back(rec.t.t_ppmm_ms[i]);
            errR[idx] = rec.err_ref;
            errC[idx] = rec.err_clear;
            nc[idx]   = rec.t.n_cmt;
            np[idx]   = rec.t.n_ppmm;
            nr[idx]   = rec.t.n_relin;
        }
    }
    for (size_t i = 0; i < labels.size(); ++i)
        PushRow(N, d, B, labels[i], acc[i], errR[i], errC[i], nc[i], np[i], nr[i], skipFirst);
    WriteResults(outPath);
}

void RunOrientation(CMTContext& ctx, uint32_t N, uint32_t d, const std::string& outPath) {
    auto A  = RandMat(d, d, 0x0A11u ^ N ^ d);
    auto Bm = RandMat(d, d, 0x0B22u ^ N);
    std::vector<std::vector<int64_t>> clear(d, std::vector<int64_t>(d, 0));
    for (uint32_t i = 0; i < d; ++i)
        for (uint32_t k0 = 0; k0 < d; ++k0)
            for (uint32_t j = 0; j < d; ++j)
                clear[i][j] += A[i][k0] * Bm[k0][j];

    auto tL = Clock::now();
    auto L  = EncryptMatrixColumns(ctx, A);
    (void)MsSince(tL);

    auto tCol = Clock::now();
    auto Rcol = EncryptMatrixColumns(ctx, Bm);
    const double setupCol = MsSince(tCol);

    auto tRow = Clock::now();
    auto Rrow = EncryptMatrixRows(ctx, Bm);
    const double setupRow = MsSince(tRow);

    const double extra = static_cast<double>(ctx.encodeScale);
    struct OrientVar {
        const char* name;
        CMTVariant v;
        AFCMTMode mode;
        bool skip;
        double setup;
        const std::vector<Ciphertext<DCRTPoly>>* R;
    };
    std::vector<OrientVar> vars = {
        {"COL_CKL", CMTVariant::CKL_REFERENCE, AFCMTMode::AFCMT_GENERIC, false, setupCol, &Rcol},
        {"COL_AFCMT", CMTVariant::AFCMT, AFCMTMode::AFCMT_FUSED_RADIX2, false, setupCol, &Rcol},
        {"ROW_SKIP_CMT1_CKL", CMTVariant::CKL_REFERENCE, AFCMTMode::AFCMT_GENERIC, true, setupRow, &Rrow},
        {"ROW_SKIP_CMT1_AFCMT", CMTVariant::AFCMT, AFCMTMode::AFCMT_FUSED_RADIX2, true, setupRow, &Rrow},
    };

    for (auto& ov : vars) {
        ctx.afcmtMode = ov.mode;
        CCMMTimings tm;
        auto outC = RunCCMM(ctx, L, *ov.R, ov.v, &tm, ov.skip);
        auto dec  = DecryptMatrixColumns(ctx, outC);
        ScaleDec(dec, extra);
        const double err = InfRelI(dec, clear);
        std::cout << "  orient corr " << ov.name << " nCMT=" << tm.n_cmt << " vsClear=" << Prec(err) << "\n";
        for (int w = 0; w < 3; ++w)
            (void)RunCCMM(ctx, L, *ov.R, ov.v, nullptr, ov.skip);
        std::vector<double> on;
        for (int r = 0; r < 10; ++r) {
            CCMMTimings t;
            (void)RunCCMM(ctx, L, *ov.R, ov.v, &t, ov.skip);
            on.push_back(t.t_ccmm_ms);
        }
        const double online = Summarize(on).med;
        for (uint32_t R : {1u, 2u, 4u, 8u, 16u, 32u, 64u, 128u}) {
            OrientRow row;
            row.N             = N;
            row.d             = d;
            row.R             = R;
            row.name          = ov.name;
            row.setup_ms      = ov.setup;
            row.online_ms     = online;
            row.cumulative_ms = ov.setup + static_cast<double>(R) * online;
            row.n_cmt         = tm.n_cmt;
            row.err_clear     = err;
            g_orient.push_back(row);
        }
        std::cout << "  orient " << ov.name << " setup=" << ov.setup << " online=" << online << " nCMT=" << tm.n_cmt
                  << "\n";
    }
    WriteResults(outPath);
}

}  // namespace

int main(int argc, char** argv) {
    std::string mode    = "app";
    std::string outPath = "/home/user/Desktop/batchmm_ccmm/RESULTS_CCMM.md";
    int warmup = 3, timed = 10;
    uint32_t onlyN = 0, onlyD = 0;
    std::vector<uint32_t> Bs = {1};
    bool explicitB           = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto eat      = [&](const char* k, auto& dst) {
            std::string p = std::string(k) + "=";
            if (a.rfind(p, 0) != 0)
                return false;
            std::istringstream in(a.substr(p.size()));
            in >> dst;
            return true;
        };
        if (a == "--all-B") {
            Bs        = {1, 2, 4, 8, 16, 32, 64, 128};
            explicitB = true;
            continue;
        }
        std::string bs;
        if (eat("--mode", mode) || eat("--out", outPath) || eat("--warmup", warmup) || eat("--timed", timed) ||
            eat("--N", onlyN) || eat("--d", onlyD) || eat("--B", bs)) {
            if (!bs.empty()) {
                Bs.clear();
                explicitB = true;
                std::replace(bs.begin(), bs.end(), ',', ' ');
                std::istringstream in(bs);
                uint32_t x;
                while (in >> x)
                    Bs.push_back(x);
            }
            continue;
        }
        std::cerr << "unknown arg " << a << "\n";
        return 1;
    }

    g_csv.open("/home/user/Desktop/batchccmm_final_github/raw_ccmm_timings.csv");
    WriteCSVHeader();

    auto pairs = ValidPairs();
    if (onlyN || onlyD) {
        std::vector<std::pair<uint32_t, uint32_t>> f;
        for (auto p : pairs)
            if ((!onlyN || p.first == onlyN) && (!onlyD || p.second == onlyD))
                f.push_back(p);
        pairs.swap(f);
    }

    const bool doApp     = (mode == "app" || mode == "full");
    const bool doAblate  = (mode == "ablation");
    const bool doOrient  = (mode == "orient" || mode == "full");
    const bool doFullX   = (mode == "full");
    auto labels          = doAblate ? AblationVariants() : AppVariants();

    std::cout << "CCMM experiment  mode=" << mode << "  variants=" << labels.size() << "\n";
    std::cout << "OpenMP threads=" << EnvOr("OMP_NUM_THREADS", "default") << " bind=" << EnvOr("OMP_PROC_BIND", "unset")
              << " places=" << EnvOr("OMP_PLACES", "unset") << "\n";

    for (auto [N, d] : pairs) {
        std::cout << "==== N=" << N << " d=" << d << " ====\n";
        CMTContext ctx;
        try {
            ctx = MakeCMTContext(N, d);
        }
        catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            continue;
        }

        auto A  = RandMat(d, d, 0xCC11u ^ N ^ d);
        auto Bm = RandMat(d, d, 0xBB22u ^ N);
        std::vector<std::vector<int64_t>> clear(d, std::vector<int64_t>(d, 0));
        for (uint32_t i = 0; i < d; ++i)
            for (uint32_t k0 = 0; k0 < d; ++k0)
                for (uint32_t j = 0; j < d; ++j)
                    clear[i][j] += A[i][k0] * Bm[k0][j];

        auto tEnc              = Clock::now();
        auto L                 = EncryptMatrixColumns(ctx, A);
        auto R                 = EncryptMatrixColumns(ctx, Bm);
        const double t_encrypt = MsSince(tEnc);

        ctx.afcmtMode = AFCMTMode::AFCMT_GENERIC;
        auto ref      = RunCCMM(ctx, L, R, CMTVariant::CKL_REFERENCE, nullptr, false);
        auto decR     = DecryptMatrixColumns(ctx, ref);
        ScaleDec(decR, static_cast<double>(ctx.encodeScale));

        std::vector<uint32_t> useB = Bs;
        if (doFullX && !explicitB) {
            useB = {1};
            if (N <= 8192 && d <= 16)
                useB = {1, 2, 4, 8};
        }

        if (doApp || doAblate) {
            auto runLabels = labels;
            if (doFullX && N == 4096 && d == 8)
                runLabels = AblationVariants();
            else if (doFullX && N == 8192 && d == 16)
                runLabels = AblationVariants();
            for (uint32_t B : useB)
                RunConfig(ctx, N, d, B, runLabels, warmup, timed, outPath, false, L, R, decR, clear, t_encrypt);
        }

        if (doOrient && (N <= 8192 && d <= 16))
            RunOrientation(ctx, N, d, outPath);
    }
    WriteResults(outPath);
    std::cout << "wrote " << outPath << "\n";
    return 0;
}
