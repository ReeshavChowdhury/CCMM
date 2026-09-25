//==================================================================================
// Proposed CCMM system vs CKL-reference CCMM.
// CKL keeps the original sequential 4x PPMM.
// Proposed = AFCMT CMT + fused/shared-NTT PPMM + OpenMP transpose/relin.
// Same 3 CMT + 4 PPMM products + relin + accumulate arithmetic.
//==================================================================================

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-cmt.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace lbcrypto;

namespace {

std::vector<std::vector<int64_t>> RandMat(uint32_t r, uint32_t c, uint32_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int64_t> dist(-7, 7);
    std::vector<std::vector<int64_t>> M(r, std::vector<int64_t>(c));
    for (uint32_t i = 0; i < r; ++i)
        for (uint32_t j = 0; j < c; ++j)
            M[i][j] = dist(rng);
    return M;
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
double InfRelD(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    double n = 0, d = 0;
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < A[i].size(); ++j) {
            n = std::max(n, std::abs(A[i][j] - B[i][j]));
            d = std::max(d, std::abs(B[i][j]));
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
    double med{0}, iqr{0};
};
Stats Summarize(std::vector<double> v) {
    Stats s;
    if (v.empty())
        return s;
    std::sort(v.begin(), v.end());
    s.med = v[v.size() / 2];
    s.iqr = v[(3 * v.size()) / 4] - v[v.size() / 4];
    return s;
}

void ScaleDec(std::vector<std::vector<double>>& M, double extra) {
    for (auto& row : M)
        for (auto& x : row)
            x /= extra;
}

struct Label {
    const char* name;
    CMTVariant v;
    AFCMTMode mode;
    bool fast;
};

struct Row {
    uint32_t N{0}, d{0};
    std::string name;
    Stats tot, cmt, ppmm, relin, layout;
    double err_ref{0}, err_clear{0};
    uint32_t n_cmt{0}, n_ppmm{0};
};

std::vector<Row> g_rows;

bool CKLOk(uint32_t N, uint32_t d) {
    return CKLConditionSatisfied(N, d);
}

}  // namespace

int main(int argc, char** argv) {
    uint32_t onlyN = 0, onlyD = 0;
    int warmup = 3, timed = 10;
    std::string outPath = "/home/user/Desktop/batchmm_ccmm/RESULTS_CCMM_FAIR.md";
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
        if (!(eat("--N", onlyN) || eat("--d", onlyD) || eat("--warmup", warmup) || eat("--timed", timed) ||
              eat("--out", outPath))) {
            std::cerr << "unknown arg " << a << "\n";
            return 1;
        }
    }

    // Paper Alg. 4 CCMM vs ours. Same 3 CMT + 4 R_k PPMM + relin + accumulate.
    // Same fused OpenMP PPMM. Only the CMT backend changes.
    const std::vector<Label> labels = {
        {"PAPER_CCMM", CMTVariant::CKL_REFERENCE, AFCMTMode::AFCMT_GENERIC, true},
        {"OURS_CCMM", CMTVariant::AFCMT, AFCMTMode::AFCMT_FUSED_RADIX2, true},
    };

    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    for (uint32_t N : {4096u, 8192u, 16384u})
        for (uint32_t d : {8u, 16u, 32u, 64u, 128u})
            if (CKLOk(N, d) && (!onlyN || N == onlyN) && (!onlyD || d == onlyD))
                pairs.emplace_back(N, d);

    std::cout << "CCMM-OPT  variants=" << labels.size() << " pairs=" << pairs.size()
              << " OMP=" << (std::getenv("OMP_NUM_THREADS") ? std::getenv("OMP_NUM_THREADS") : "unset") << "\n";

    auto write = [&]() {
        std::ofstream out(outPath);
        out << "# RESULTS_CCMM_OPT\n\n";
        out << "Paper Algorithm 4 (Cheon–Kang–Lee `batcmm.pdf`) vs our CCMM.\n";
        out << "Same pipeline: 3 CMT + 4 PPMM over `R_k` + relin of the `sk^2` term + accumulate.\n";
        out << "Same fused OpenMP PPMM / transpose / relin. Only CMT backend changes "
               "(paper `EvalAutomorphism` TWEAK vs AFCMT fused consume).\n";
        out << "Decrypt must match paper CCMM and cleartext AB.\n\n";
        out << "| CPU | Intel Core i9-10900 @ 2.80 GHz, 20 logical processors |\n";
        out << "| OpenMP | OMP_NUM_THREADS="
            << (std::getenv("OMP_NUM_THREADS") ? std::getenv("OMP_NUM_THREADS") : "unset") << " |\n\n";

        out << "## End-to-end (median T_CCMM, ms)\n\n";
        out << "| N | d | PAPER_CCMM | OURS_CCMM | speedup | CMT-sum paper | CMT-sum ours | CMT speedup |\n";
        out << "|--:|--:|-----------:|----------:|--------:|--------------:|-------------:|------------:|\n";
        auto find = [&](uint32_t N, uint32_t d, const char* n) -> const Row* {
            for (const auto& r : g_rows)
                if (r.N == N && r.d == d && r.name == n)
                    return &r;
            return nullptr;
        };
        std::vector<std::pair<uint32_t, uint32_t>> seen;
        for (const auto& r : g_rows) {
            bool dup = false;
            for (auto p : seen)
                if (p.first == r.N && p.second == r.d)
                    dup = true;
            if (!dup)
                seen.emplace_back(r.N, r.d);
        }
        double best = 0, medacc = 0;
        uint32_t bN = 0, bD = 0;
        std::vector<double> sps;
        for (auto p : seen) {
            auto* ckl = find(p.first, p.second, "PAPER_CCMM");
            auto* pr  = find(p.first, p.second, "OURS_CCMM");
            if (!ckl || !pr)
                continue;
            const double s  = (pr->tot.med > 0) ? ckl->tot.med / pr->tot.med : 0;
            const double cs = (pr->cmt.med > 0) ? ckl->cmt.med / pr->cmt.med : 0;
            out << "| " << p.first << " | " << p.second << " | " << std::fixed << std::setprecision(2) << ckl->tot.med
                << " | " << pr->tot.med << " | " << std::setprecision(3) << s << " | " << ckl->cmt.med << " | "
                << pr->cmt.med << " | " << cs << " |\n";
            sps.push_back(s);
            if (s > best) {
                best = s;
                bN   = p.first;
                bD   = p.second;
            }
        }
        if (!sps.empty()) {
            std::sort(sps.begin(), sps.end());
            medacc = sps[sps.size() / 2];
        }

        out << "\n## Component medians (ms; not additive except exclusive same-rep sums)\n\n";
        out << "| N | d | variant | CMT-sum | PPMM-sum | relin | layout | T_CCMM |\n";
        out << "|--:|--:|:--------|--------:|---------:|------:|-------:|-------:|\n";
        for (const auto& r : g_rows)
            out << "| " << r.N << " | " << r.d << " | " << r.name << " | " << std::fixed << std::setprecision(2)
                << r.cmt.med << " | " << r.ppmm.med << " | " << r.relin.med << " | " << r.layout.med << " | "
                << r.tot.med << " |\n";

        out << "\n## Correctness\n\n";
        out << "| N | d | variant | vs CKL | vs clear | n_cmt | n_ppmm |\n";
        out << "|--:|--:|:--------|-------:|---------:|------:|-------:|\n";
        bool ok = !g_rows.empty();
        for (const auto& r : g_rows) {
            out << "| " << r.N << " | " << r.d << " | " << r.name << " | " << Prec(r.err_ref) << " | "
                << Prec(r.err_clear) << " | " << r.n_cmt << " | " << r.n_ppmm << " |\n";
            if (r.n_cmt != 3 || r.n_ppmm != 4 || r.err_clear > 1e-6)
                ok = false;
        }

        out << "\n## Summary\n\n";
        out << "- Best e2e OURS vs PAPER (same PPMM): **" << std::setprecision(3) << best << "×** at N=" << bN
            << " d=" << bD << "\n";
        out << "- Median e2e OURS vs PAPER: **" << medacc << "×**\n";
        out << "- Reached 2× e2e: **" << (best >= 2.0 ? "YES" : "NO") << "**\n";
        out << "- Correctness: " << (ok ? "PASS" : "FAIL") << "\n";
        out << "- Fair: PPMM / layout / relin / OpenMP / keys / params identical. Only CMT backend changes.\n";
    };

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
        auto A  = RandMat(d, d, 0xA11u ^ N ^ d);
        auto Bm = RandMat(d, d, 0xB22u ^ N);
        std::vector<std::vector<int64_t>> clear(d, std::vector<int64_t>(d, 0));
        for (uint32_t i = 0; i < d; ++i)
            for (uint32_t k0 = 0; k0 < d; ++k0)
                for (uint32_t j = 0; j < d; ++j)
                    clear[i][j] += A[i][k0] * Bm[k0][j];
        auto L = EncryptMatrixColumns(ctx, A);
        auto R = EncryptMatrixColumns(ctx, Bm);

        ctx.afcmtMode = AFCMTMode::AFCMT_GENERIC;
        auto ref      = RunCCMM(ctx, L, R, CMTVariant::CKL_REFERENCE, nullptr, false, false);
        auto decR     = DecryptMatrixColumns(ctx, ref);
        ScaleDec(decR, static_cast<double>(ctx.encodeScale));

        std::mt19937 rng(N * 13u + d);
        std::vector<std::vector<double>> tot(labels.size()), cmt(labels.size()), pp(labels.size()), rel(labels.size()),
            lay(labels.size());
        std::vector<double> errR(labels.size()), errC(labels.size());
        std::vector<uint32_t> nc(labels.size()), np(labels.size());

        for (auto lv : labels) {
            ctx.afcmtMode = lv.mode;
            CCMMTimings tm;
            auto out = RunCCMM(ctx, L, R, lv.v, &tm, false, lv.fast);
            auto dec = DecryptMatrixColumns(ctx, out);
            ScaleDec(dec, static_cast<double>(ctx.encodeScale));
            std::cout << "  corr " << lv.name << " vsCKL=" << Prec(InfRelD(dec, decR))
                      << " vsClear=" << Prec(InfRelI(dec, clear)) << " nCMT=" << tm.n_cmt << " nPPMM=" << tm.n_ppmm
                      << "\n";
        }

        for (int w = 0; w < warmup; ++w) {
            auto lab = labels;
            std::shuffle(lab.begin(), lab.end(), rng);
            for (auto lv : lab) {
                ctx.afcmtMode = lv.mode;
                (void)RunCCMM(ctx, L, R, lv.v, nullptr, false, lv.fast);
            }
        }
        for (int r = 0; r < timed; ++r) {
            auto lab = labels;
            std::shuffle(lab.begin(), lab.end(), rng);
            for (auto lv : lab) {
                ctx.afcmtMode = lv.mode;
                CCMMTimings tm;
                auto out = RunCCMM(ctx, L, R, lv.v, &tm, false, lv.fast);
                auto dec = DecryptMatrixColumns(ctx, out);
                ScaleDec(dec, static_cast<double>(ctx.encodeScale));
                size_t idx = 0;
                for (size_t i = 0; i < labels.size(); ++i)
                    if (std::string(labels[i].name) == lv.name)
                        idx = i;
                tot[idx].push_back(tm.t_ccmm_ms);
                cmt[idx].push_back(tm.t_cmt_sum_ms);
                pp[idx].push_back(tm.t_ppmm_sum_ms);
                rel[idx].push_back(tm.t_relin_ms);
                lay[idx].push_back(tm.t_layout_ms);
                errR[idx] = InfRelD(dec, decR);
                errC[idx] = InfRelI(dec, clear);
                nc[idx]   = tm.n_cmt;
                np[idx]   = tm.n_ppmm;
            }
        }
        for (size_t i = 0; i < labels.size(); ++i) {
            Row row;
            row.N         = N;
            row.d         = d;
            row.name      = labels[i].name;
            row.tot       = Summarize(tot[i]);
            row.cmt       = Summarize(cmt[i]);
            row.ppmm      = Summarize(pp[i]);
            row.relin     = Summarize(rel[i]);
            row.layout    = Summarize(lay[i]);
            row.err_ref   = errR[i];
            row.err_clear = errC[i];
            row.n_cmt     = nc[i];
            row.n_ppmm    = np[i];
            g_rows.push_back(row);
            std::cout << "  " << row.name << " med=" << row.tot.med << " iqr=" << row.tot.iqr
                      << " cmt=" << row.cmt.med << " ppmm=" << row.ppmm.med << "\n";
        }
        write();
    }
    write();
    std::cout << "wrote " << outPath << "\n";
    return 0;
}
