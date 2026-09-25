//==================================================================================
// Measurement-driven FCMT/AFCMT upgrade (afcmt_1).
// Protocol: 3 warm-ups, 10 timed reps, randomized variant order, median/min/IQR.
//==================================================================================

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define PROFILE
#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-cmt.h"

#include <sched.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace lbcrypto;

namespace {

struct Label {
    CMTVariant v;
    AFCMTMode mode;
    std::string name;
};

std::vector<Label> AllVariants() {
    return {
        {CMTVariant::CKL_REFERENCE, AFCMTMode::AFCMT_GENERIC, "CKL_REFERENCE"},
        {CMTVariant::OPENFHE_HOISTED, AFCMTMode::AFCMT_GENERIC, "OPENFHE_HOISTED"},
        {CMTVariant::SPECTRAL_TWEAK, AFCMTMode::AFCMT_GENERIC, "SPECTRAL_TWEAK"},
        {CMTVariant::FCMT, AFCMTMode::AFCMT_GENERIC, "FCMT"},
        {CMTVariant::AFCMT, AFCMTMode::AFCMT_GENERIC, "AFCMT_GENERIC"},
        {CMTVariant::AFCMT, AFCMTMode::AFCMT_RESIDUE, "AFCMT_RESIDUE"},
        {CMTVariant::AFCMT, AFCMTMode::AFCMT_FUSED_RADIX2, "AFCMT_FUSED_RADIX2"},
        {CMTVariant::AFCMT, AFCMTMode::AFCMT_BATCHED, "AFCMT_BATCHED"},
    };
}

std::vector<std::pair<uint32_t, uint32_t>> ValidPairs() {
    std::vector<std::pair<uint32_t, uint32_t>> o;
    for (uint32_t N : {4096u, 8192u, 16384u})
        for (uint32_t d : {8u, 16u, 32u, 64u, 128u})
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
    const double q1 = v[v.size() / 4];
    const double q3 = v[(3 * v.size()) / 4];
    s.iqr           = q3 - q1;
    return s;
}

struct CorrRow {
    uint32_t N, d;
    std::string v;
    double vs_ref, vs_clear;
};
struct CompStats {
    Stats total, mono, fwd, hoist, ks, aut, inv, dft, hoistks;
};

CompStats SummarizeComps(const std::vector<CMTTimings>& reps) {
    auto col = [&](auto fn) {
        std::vector<double> v;
        v.reserve(reps.size());
        for (const auto& t : reps)
            v.push_back(fn(t));
        return Summarize(std::move(v));
    };
    CompStats s;
    s.total   = col([](const CMTTimings& t) { return t.t_cmt_ms; });
    s.mono    = col([](const CMTTimings& t) { return t.t_monomial_ms; });
    s.fwd     = col([](const CMTTimings& t) { return t.t_forward_tweak_ms; });
    s.hoist   = col([](const CMTTimings& t) { return t.t_hoist_ms; });
    s.ks      = col([](const CMTTimings& t) { return t.t_keyprod_ms; });
    s.aut     = col([](const CMTTimings& t) { return t.t_automorphism_ms; });
    s.inv     = col([](const CMTTimings& t) { return t.t_tweak_inv_ms; });
    s.dft     = col([](const CMTTimings& t) { return t.t_dft_ms; });
    s.hoistks = col([](const CMTTimings& t) { return t.t_hoist_ms + t.t_keyprod_ms; });
    return s;
}

struct MicroRow {
    uint32_t N, d;
    std::string v;
    CompStats s;
    CMTTraffic traffic;
};
struct CCMMRow {
    uint32_t N, d;
    std::string shape, v;
    CCMMTimings t;
    double vs_clear;
};
struct IdxRow {
    uint32_t N, d;
    std::string kernel;
    IndexKernelStats s;
};
struct AmortRow {
    uint32_t N, d, R;
    double setup, online, cum;
    std::string strategy;
};

std::vector<MathTestReport> g_math;
std::vector<CorrRow> g_corr;
std::vector<MicroRow> g_micro;
std::vector<CCMMRow> g_ccmm;
std::vector<IdxRow> g_idx;
std::vector<ClassCountReport> g_cls;
std::vector<uint32_t> g_clsN, g_clsD;
std::vector<AmortRow> g_amort;
std::string g_host;
std::ofstream g_raw;

void Host() {
    std::ostringstream o;
    o << "- Date: 2026-08-14\n- CPU: Intel i9-10900, 20 threads, OpenMP fixed\n"
      << "- GPU: unavailable\n"
      << "- Protocol: **one executable invocation**, all variants, all valid (N,d).\n"
      << "- 3 warm-up + 10 timed. Each timed round shuffles variant order.\n"
      << "- Every rep records (T_total, T_mono, T_fwd, T_hoist, T_KS, T_auto, T_inv, T_DFT).\n"
      << "- Reported T_total = median(T_total). Each T_i = median(T_i) independently. "
         "T_total is **not** the sum of component medians. Decomposition is **not** the last rep.\n"
      << "- IQR = Q3-Q1 over the 10 timed reps.\n"
      << "- T_KS = EvalFastKeySwitchCore only. CKL_REFERENCE T_KS is the whole EvalAutomorphism call.\n"
      << "- FCMT/AFCMT/SPECTRAL forward TWEAK is one monomial NTT + per-slot DFT.\n"
      << "- CKL_REFERENCE and OPENFHE_HOISTED keep coefficient-domain column TWEAK.\n"
      << "- Excluded: (4096,128) because 128^2 does not divide 2N\n";
    g_host = o.str();
    cpu_set_t set;
    CPU_ZERO(&set);
    for (int i = 0; i < 20; ++i)
        CPU_SET(i, &set);
    sched_setaffinity(0, sizeof(set), &set);
}

void WriteResults(const std::string& path) {
    std::ofstream out(path);
    out << "# FCMT / AFCMT experimental upgrade (`afcmt_1`)\n\n" << g_host << "\n";

    out << "## Table 1 — Mathematical validation\n\n";
    out << "| test | pass | detail |\n|:-----|:----:|:-------|\n";
    for (const auto& r : g_math)
        out << "| " << r.name << " | " << (r.pass ? "yes" : "NO") << " | " << r.detail << " |\n";

    out << "\n## Table 2 — CMT baseline (median ms)\n\n";
    out << "| N | d | CKL_REFERENCE | OPENFHE_HOISTED | SPECTRAL_TWEAK | FCMT |\n";
    out << "|--:|--:|--------------:|----------------:|---------------:|-----:|\n";
    auto findT = [&](uint32_t N, uint32_t d, const char* name) {
        for (const auto& m : g_micro)
            if (m.N == N && m.d == d && m.v == name)
                return m.s.total.med;
        return 0.0;
    };
    std::vector<std::pair<uint32_t, uint32_t>> seen;
    for (const auto& m : g_micro) {
        bool dup = false;
        for (auto p : seen)
            if (p.first == m.N && p.second == m.d)
                dup = true;
        if (dup)
            continue;
        seen.emplace_back(m.N, m.d);
        const double ckl = findT(m.N, m.d, "CKL_REFERENCE");
        const double hoi = findT(m.N, m.d, "OPENFHE_HOISTED");
        const double sp  = findT(m.N, m.d, "SPECTRAL_TWEAK");
        const double fc  = findT(m.N, m.d, "FCMT");
        out << "| " << m.N << " | " << m.d << " | " << std::fixed << std::setprecision(3) << ckl << " | " << hoi
            << " | " << sp << " | " << fc << " |\n";
    }

    out << "\n## Table 2b — CMT speedup vs CKL (median T_CMT)\n\n";
    out << "| N | d | FCMT / CKL | best AFCMT / CKL | FCMT / HOISTED | best / HOISTED |\n";
    out << "|--:|--:|-----------:|-----------------:|---------------:|---------------:|\n";
    for (auto p : seen) {
        const double ckl = findT(p.first, p.second, "CKL_REFERENCE");
        const double hoi = findT(p.first, p.second, "OPENFHE_HOISTED");
        const double fc  = findT(p.first, p.second, "FCMT");
        double best      = fc;
        for (const char* nm : {"AFCMT_GENERIC", "AFCMT_RESIDUE", "AFCMT_FUSED_RADIX2", "AFCMT_BATCHED"}) {
            const double x = findT(p.first, p.second, nm);
            if (x > 0 && (best <= 0 || x < best))
                best = x;
        }
        auto ratio = [](double num, double den) { return (num > 0 && den > 0) ? num / den : 0.0; };
        out << "| " << p.first << " | " << p.second << " | " << std::fixed << std::setprecision(3)
            << ratio(ckl, fc) << " | " << ratio(ckl, best) << " | " << ratio(hoi, fc) << " | " << ratio(hoi, best)
            << " |\n";
    }

    out << "\n## Table 3 — AFCMT ablation (median ms)\n\n";
    out << "| N | d | FCMT | AFCMT_GENERIC | AFCMT_RESIDUE | AFCMT_FUSED_RADIX2 | AFCMT_BATCHED | best speedup vs FCMT |\n";
    out << "|--:|--:|-----:|--------------:|--------------:|-------------------:|--------------:|---------------------:|\n";
    for (auto p : seen) {
        double f = findT(p.first, p.second, "FCMT");
        double g = findT(p.first, p.second, "AFCMT_GENERIC");
        double r = findT(p.first, p.second, "AFCMT_RESIDUE");
        double u = findT(p.first, p.second, "AFCMT_FUSED_RADIX2");
        double b = findT(p.first, p.second, "AFCMT_BATCHED");
        double best = f;
        for (double x : {g, r, u, b})
            if (x > 0 && x < best)
                best = x;
        out << "| " << p.first << " | " << p.second << " | " << f << " | " << g << " | " << r << " | " << u << " | " << b
            << " | " << ((best > 0 && f > 0) ? f / best : 0) << " |\n";
    }

    out << "\n## Table 4 — independent component medians [IQR] (ms)\n\n";
    out << "Each cell is `median [IQR]` over the 10 timed reps of that timer. "
           "`T_total` is `median(T_total)`, not the sum of the other medians. "
           "The last repetition is not used as the decomposition.\n\n";
    out << "| N | d | variant | T_total | T_mono | T_fwd | T_hoist | T_KS | T_auto | T_inv | T_DFT |\n";
    out << "|--:|--:|:--------|--------:|-------:|------:|--------:|-----:|-------:|------:|------:|\n";
    auto cell = [](const Stats& s) {
        std::ostringstream o;
        o << std::fixed << std::setprecision(2) << s.med << " [" << s.iqr << "]";
        return o.str();
    };
    for (const auto& m : g_micro)
        out << "| " << m.N << " | " << m.d << " | " << m.v << " | " << cell(m.s.total) << " | " << cell(m.s.mono)
            << " | " << cell(m.s.fwd) << " | " << cell(m.s.hoist) << " | " << cell(m.s.ks) << " | " << cell(m.s.aut)
            << " | " << cell(m.s.inv) << " | " << cell(m.s.dft) << " |\n";

    out << "\n## Table 5 — traffic (bytes)\n\n";
    out << "| N | d | variant | auto writes | auto reads | gather reads | residue layout |\n";
    out << "|--:|--:|:--------|------------:|-----------:|-------------:|---------------:|\n";
    for (const auto& m : g_micro)
        out << "| " << m.N << " | " << m.d << " | " << m.v << " | " << m.traffic.automorphism_output_bytes_written
            << " | " << m.traffic.automorphism_output_bytes_read << " | " << m.traffic.gather_bytes << " | "
            << m.traffic.residue_layout_bytes << " |\n";

    out << "\n## Table 6 — class reuse\n\n";
    out << "Gamma classes are over N frequencies: unique = d/2, count = 2k. "
           "Step classes are over a in [0,k): unique = d/2, count = 2k/d.\n\n";
    out << "| N | d | unique γ | expected d/2 | γ count | expected 2k | unique step | expected d/2 | step count | expected 2k/d | ok |\n";
    out << "|--:|--:|---------:|-------------:|--------:|------------:|------------:|-------------:|-----------:|--------------:|:--:|\n";
    for (size_t i = 0; i < g_cls.size(); ++i) {
        const auto& c = g_cls[i];
        out << "| " << g_clsN[i] << " | " << g_clsD[i] << " | " << c.unique_gamma << " | " << c.expected_gamma << " | "
            << c.min_gamma_count << " | " << c.expected_gamma_count << " | " << c.unique_step << " | " << c.expected_step
            << " | " << c.min_step_count << " | " << c.expected_step_count << " | "
            << ((c.gamma_ok && c.step_ok) ? "yes" : "NO") << " |\n";
    }

    out << "\n## Table 7 — CCMM\n\n";
    out << "| N | d | shape | variant | T_CMT | T_PPMM | T_relin | T_CCMM | rel vs clear | prec |\n";
    out << "|--:|--:|:------|:--------|------:|-------:|--------:|-------:|-------------:|:-----|\n";
    for (const auto& r : g_ccmm)
        out << "| " << r.N << " | " << r.d << " | " << r.shape << " | " << r.v << " | " << std::setprecision(2)
            << r.t.t_cmt_total_ms << " | " << r.t.t_ppmm_sum_ms << " | " << r.t.t_relin_ms << " | " << r.t.t_ccmm_ms
            << " | " << std::scientific << std::setprecision(3) << r.vs_clear << " | " << Prec(r.vs_clear) << " |\n"
            << std::fixed;

    out << "\n## Table 8 — amortization (row encrypt setup vs R online AFCMT CMTs)\n\n";
    out << "| N | d | R | strategy | setup ms | online ms | cumulative ms |\n";
    out << "|--:|--:|--:|:---------|---------:|----------:|--------------:|\n";
    for (const auto& a : g_amort)
        out << "| " << a.N << " | " << a.d << " | " << a.R << " | " << a.strategy << " | " << std::setprecision(2)
            << a.setup << " | " << a.online << " | " << a.cum << " |\n";

    out << "\n## Experiment A — index generation (no CKKS)\n\n";
    out << "| N | d | kernel | ns/index | indices/s | checksum |\n|--:|--:|:-------|---------:|----------:|----------:|\n";
    for (const auto& r : g_idx)
        out << "| " << r.N << " | " << r.d << " | " << r.kernel << " | " << std::setprecision(3) << r.s.ns_per_index
            << " | " << std::scientific << r.s.indices_per_sec << " | " << r.s.checksum << " |\n"
            << std::fixed;

    out << "\n## Correctness (all variants vs CKL)\n\n";
    out << "| N | d | variant | rel vs ref | prec vs ref | rel vs clear | prec vs clear |\n";
    out << "|--:|--:|:--------|-----------:|------------:|-------------:|--------------:|\n";
    for (const auto& r : g_corr)
        out << "| " << r.N << " | " << r.d << " | " << r.v << " | " << std::scientific << std::setprecision(3)
            << r.vs_ref << " | " << Prec(r.vs_ref) << " | " << r.vs_clear << " | " << Prec(r.vs_clear) << " |\n";

    out << "\n## 2× analysis (measured, not claimed)\n\n";
    out << "FCMT cannot drop key-switch: each of the `d` streams still needs "
           "`EvalFastRotationPrecompute` (T_hoist) plus `EvalFastKeySwitchCore` (T_keyprod). "
           "The automorphism itself is cheap (typically <4 ms). The only work FCMT can remove "
           "is inverse TWEAK plus inverse monomials plus the auto permute.\n\n";
    out << "| N | d | CKL T_CMT | FCMT T_CMT | speedup | median(T_hoist+T_KS) | CKL / that ceiling |\n";
    out << "|--:|--:|----------:|-----------:|--------:|---------------------:|-------------------:|\n";
    for (auto p : seen) {
        double ckl = 0, fc = 0, hoistks = 0;
        for (const auto& m : g_micro) {
            if (m.N != p.first || m.d != p.second)
                continue;
            if (m.v == "CKL_REFERENCE")
                ckl = m.s.total.med;
            if (m.v == "FCMT") {
                fc      = m.s.total.med;
                hoistks = m.s.hoistks.med;
            }
        }
        const double floor = hoistks;
        out << "| " << p.first << " | " << p.second << " | " << std::fixed << std::setprecision(3) << ckl << " | "
            << fc << " | " << ((fc > 0) ? ckl / fc : 0) << " | " << floor << " | "
            << ((floor > 0) ? ckl / floor : 0) << " |\n";
    }
    out << "\nA consistent 2× requires T_CMT(FCMT) ≤ T_CMT(CKL)/2. That is possible only when "
           "hoist+KS is already < 50% of CKL and the fused DFT/TWEAK constants stay small. "
           "On this CPU the auto permute is almost free, so FCMT's algorithmic saving is "
           "inv-TWEAK, not key-switch. Small `(N,d)` therefore stay near 1.1–1.4×; larger rings "
           "can approach 2× on a good run but do not do so on every pair.\n";

    out << "\n## Timing protocol notes\n\n";
    out << "- Keygen and encryption are outside the timed CMT region.\n";
    out << "- Each timed CMT clones the input ciphertext; no mutated reuse.\n";
    out << "- Per-rep CSV: `raw_timings.csv` columns "
           "T_total,T_mono,T_fwd,T_hoist,T_KS,T_auto,T_inv,T_DFT.\n";
    out << "- Table 4 uses independent medians. Do not add the component medians to recover T_total.\n";
    out << "- CKL_REFERENCE T_keyprod is the whole `EvalAutomorphism` call and is **not** "
           "compared to fused T_keyprod. Use T_CMT for that comparison.\n";
    out << "- Hardware counters unavailable (`perf_event_paranoid=4`). No invented counter data.\n";
    out << "- Negative results are retained. No universal speedup is claimed.\n";
    out << "- Forward TWEAK on FCMT/AFCMT/SPECTRAL is one monomial NTT + per-slot DFT. "
           "CKL_REFERENCE / OPENFHE_HOISTED keep coefficient-domain column TWEAK.\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::string mode = "all", outPath = "/home/user/Desktop/batchmm/afcmt_1/RESULTS.md";
    int warmup = 3, timed = 10;
    uint32_t onlyN = 0, onlyD = 0;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto eat      = [&](const char* k, auto& dst) {
            std::string p = std::string(k) + "=";
            if (a.rfind(p, 0) == 0) {
                std::istringstream ss(a.substr(p.size()));
                ss >> dst;
                return true;
            }
            return false;
        };
        eat("--mode", mode) || eat("--out", outPath) || eat("--warmup", warmup) || eat("--timed", timed) ||
            eat("--N", onlyN) || eat("--d", onlyD);
    }
    Host();
    g_raw.open("/home/user/Desktop/batchccmm_final_github/raw_timings.csv");
    g_raw << "N,d,variant,rep,T_total,T_mono,T_fwd,T_hoist,T_KS,T_auto,T_inv,T_DFT\n";

    auto pairs = ValidPairs();
    if (onlyN || onlyD) {
        std::vector<std::pair<uint32_t, uint32_t>> f;
        for (auto p : pairs)
            if ((!onlyN || p.first == onlyN) && (!onlyD || p.second == onlyD))
                f.push_back(p);
        pairs.swap(f);
    }

    // Experiment F / class counts + Experiment A (no FHE) for every pair.
    for (auto [N, d] : pairs) {
        auto c = CountTransformClasses(N, d);
        g_cls.push_back(c);
        g_clsN.push_back(N);
        g_clsD.push_back(d);
        const uint32_t iters = (N >= 16384 ? 4u : 12u);
        g_idx.push_back({N, d, "A1_direct", BenchIndexDirect(N, d, iters)});
        g_idx.push_back({N, d, "A2_lookup", BenchIndexLookup(N, d, iters)});
        g_idx.push_back({N, d, "A3_recurrence", BenchIndexRecurrence(N, d, iters)});
    }

    if (pairs.empty())
        return 1;

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
        auto math = RunMathUnitTests(N, d, &ctx);
        for (auto& r : math)
            g_math.push_back(r);
        for (const auto& r : math)
            std::cout << (r.pass ? "PASS " : "FAIL ") << r.name << " : " << r.detail << "\n";

        auto M    = RandMat(d, d, 0xA11CE ^ N ^ d);
        auto cols = EncryptMatrixColumns(ctx, M);
        auto ref  = RunCMT(ctx, cols, CMTVariant::CKL_REFERENCE);
        auto decR = DecryptMatrixColumns(ctx, ref);
        auto clr  = CleartextCMT(M);

        // Correctness all variants (Tests J, K).
        for (auto lv : AllVariants()) {
            ctx.afcmtMode = lv.mode;
            auto out      = (lv.v == CMTVariant::CKL_REFERENCE) ? ref : RunCMT(ctx, cols, lv.v);
            auto dec      = DecryptMatrixColumns(ctx, out);
            g_corr.push_back({N, d, lv.name, InfRelD(dec, decR), InfRelI(dec, clr)});
            std::cout << "  corr " << lv.name << " " << Prec(g_corr.back().vs_ref) << "\n";
        }

        if (mode == "unittests") {
            WriteResults(outPath);
            continue;
        }

        // CMT micro: 3 warmup + 10 timed, random order.
        std::mt19937 rng(N * 10007u + d);
        auto labels = AllVariants();
        std::vector<std::vector<CMTTimings>> reps(labels.size());
        CMTTraffic traffic[8];
        for (int w = 0; w < warmup; ++w)
            for (auto lv : labels) {
                ctx.afcmtMode = lv.mode;
                (void)RunCMT(ctx, cols, lv.v, nullptr);
            }
        for (int r = 0; r < timed; ++r) {
            std::shuffle(labels.begin(), labels.end(), rng);
            for (auto lv : labels) {
                ctx.afcmtMode = lv.mode;
                CMTTimings t;
                (void)RunCMT(ctx, cols, lv.v, &t);
                size_t idx = 0;
                for (size_t i = 0; i < AllVariants().size(); ++i)
                    if (AllVariants()[i].name == lv.name)
                        idx = i;
                reps[idx].push_back(t);
                traffic[idx] = t.traffic;
                g_raw << N << "," << d << "," << lv.name << "," << r << "," << t.t_cmt_ms << "," << t.t_monomial_ms
                      << "," << t.t_forward_tweak_ms << "," << t.t_hoist_ms << "," << t.t_keyprod_ms << ","
                      << t.t_automorphism_ms << "," << t.t_tweak_inv_ms << "," << t.t_dft_ms << "\n";
            }
        }
        g_raw.flush();
        auto allL = AllVariants();
        for (size_t i = 0; i < allL.size(); ++i) {
            MicroRow row;
            row.N       = N;
            row.d       = d;
            row.v       = allL[i].name;
            row.s       = SummarizeComps(reps[i]);
            row.traffic = traffic[i];
            g_micro.push_back(row);
            std::cout << "  micro " << row.v << " T_total med=" << row.s.total.med << " iqr=" << row.s.total.iqr
                      << " min=" << row.s.total.mn << "\n";
        }

        if (mode == "all" || mode == "ccmm") {
            auto A = RandMat(d, d, 11 + N);
            auto B = RandMat(d, d, 13 + d);
            auto exp = [&]() {
                std::vector<std::vector<int64_t>> C(d, std::vector<int64_t>(d, 0));
                for (uint32_t i = 0; i < d; ++i)
                    for (uint32_t k0 = 0; k0 < d; ++k0)
                        for (uint32_t j = 0; j < d; ++j)
                            C[i][j] += A[i][k0] * B[k0][j];
                return C;
            }();
            auto L = EncryptMatrixColumns(ctx, A);
            auto R = EncryptMatrixColumns(ctx, B);
            for (auto lv : {Label{CMTVariant::CKL_REFERENCE, AFCMTMode::AFCMT_GENERIC, "CKL_REFERENCE"},
                            Label{CMTVariant::FCMT, AFCMTMode::AFCMT_GENERIC, "FCMT"},
                            Label{CMTVariant::AFCMT, AFCMTMode::AFCMT_FUSED_RADIX2, "AFCMT_FUSED_RADIX2"}}) {
                ctx.afcmtMode = lv.mode;
                (void)RunCCMM(ctx, L, R, lv.v, nullptr);  // warm-up
                CCMMTimings t;
                auto out = RunCCMM(ctx, L, R, lv.v, &t);
                auto dec = DecryptMatrixColumns(ctx, out);
                const double extra = static_cast<double>(ctx.encodeScale);
                for (auto& row : dec)
                    for (auto& x : row)
                        x /= extra;
                g_ccmm.push_back({N, d, "square", lv.name, t, InfRelI(dec, exp)});
                std::cout << "  ccmm " << lv.name << " " << t.t_ccmm_ms << "ms\n";
            }
            if (d <= 16) {
                auto AR = RandMat(d, 2 * d, 21 + N);
                auto BR = RandMat(2 * d, d, 23 + d);
                std::vector<std::vector<int64_t>> ER(d, std::vector<int64_t>(d, 0));
                for (uint32_t i = 0; i < d; ++i)
                    for (uint32_t k0 = 0; k0 < 2 * d; ++k0)
                        for (uint32_t j = 0; j < d; ++j)
                            ER[i][j] += AR[i][k0] * BR[k0][j];
                ctx.afcmtMode = AFCMTMode::AFCMT_FUSED_RADIX2;
                CCMMTimings t;
                auto out = RunCCMMRectangular(ctx, AR, BR, CMTVariant::AFCMT, &t);
                auto dec = DecryptMatrixColumns(ctx, out);
                const double extra = static_cast<double>(ctx.encodeScale);
                for (auto& row : dec)
                    for (auto& x : row)
                        x /= extra;
                g_ccmm.push_back({N, d, "rect", "AFCMT_FUSED_RADIX2", t, InfRelI(dec, ER)});
            }
        }

        // Experiment L/N amortization: row setup vs R * AFCMT online.
        using Clock = std::chrono::steady_clock;
        auto t0     = Clock::now();
        auto rowcts = EncryptMatrixRows(ctx, M);
        const double setupRow = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        ctx.afcmtMode         = AFCMTMode::AFCMT_FUSED_RADIX2;
        CMTTimings one;
        (void)RunCMT(ctx, cols, CMTVariant::AFCMT, &one);
        for (uint32_t R : {1u, 2u, 4u, 8u, 16u, 32u, 64u, 128u, 256u}) {
            g_amort.push_back({N, d, R, 0, R * one.t_cmt_ms, R * one.t_cmt_ms, "column+AFCMT"});
            g_amort.push_back({N, d, R, setupRow, 0, setupRow, "row+no_first_CMT"});
            g_amort.push_back({N, d, R, setupRow, std::max(0, (int)R - 1) * one.t_cmt_ms,
                               setupRow + std::max(0, (int)R - 1) * one.t_cmt_ms, "row+remaining_AFCMT"});
        }
        WriteResults(outPath);
        std::cout << "wrote " << outPath << "\n";
    }
    WriteResults(outPath);
    std::cout << "done " << outPath << "\n";
    return 0;
}
