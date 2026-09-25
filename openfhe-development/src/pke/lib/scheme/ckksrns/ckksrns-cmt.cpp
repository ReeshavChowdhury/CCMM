//==================================================================================
// CKL CMT, FCMT, and AFCMT. EvalFastKeySwitchCore is unchanged.
//==================================================================================

#include "scheme/ckksrns/ckksrns-cmt.h"

#include "math/nbtheory.h"
#include "utils/utilities.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lbcrypto {

namespace {

using Clock = std::chrono::steady_clock;

inline double MsSince(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

NativeInteger ModExp(NativeInteger base, uint32_t e, const NativeInteger& q) {
    NativeInteger r(1);
    NativeInteger b = base;
    while (e) {
        if (e & 1)
            r.ModMulEq(b, q);
        b.ModMulEq(b, q);
        e >>= 1;
    }
    return r;
}

uint32_t ModInverseU32(uint32_t a, uint32_t m) {
    int64_t t = 0, newt = 1;
    int64_t r = static_cast<int64_t>(m), newr = static_cast<int64_t>(a);
    while (newr != 0) {
        const int64_t q = r / newr;
        int64_t tmp     = newt;
        newt            = t - q * newt;
        t               = tmp;
        tmp             = newr;
        newr            = r - q * newr;
        r               = tmp;
    }
    if (r > 1)
        OPENFHE_THROW("ModInverseU32: not invertible");
    if (t < 0)
        t += m;
    return static_cast<uint32_t>(t);
}

void SetAllFormat(std::vector<Ciphertext<DCRTPoly>>& cts, Format fmt) {
    for (auto& ct : cts)
        for (auto& el : ct->GetElements())
            el.SetFormat(fmt);
}

Ciphertext<DCRTPoly> AddCT(const Ciphertext<DCRTPoly>& a, const Ciphertext<DCRTPoly>& b) {
    auto r         = a->Clone();
    auto& re       = r->GetElements();
    const auto& be = b->GetElements();
    for (size_t i = 0; i < re.size(); ++i)
        re[i] += be[i];
    return r;
}

Ciphertext<DCRTPoly> SubCT(const Ciphertext<DCRTPoly>& a, const Ciphertext<DCRTPoly>& b) {
    auto r         = a->Clone();
    auto& re       = r->GetElements();
    const auto& be = b->GetElements();
    for (size_t i = 0; i < re.size(); ++i)
        re[i] -= be[i];
    return r;
}

void NegacyclicMulXInPlace(DCRTPoly& p, int64_t power) {
    p.SetFormat(Format::COEFFICIENT);
    const uint32_t N   = p.GetRingDimension();
    const int64_t twoN = 2 * static_cast<int64_t>(N);
    int64_t e          = power % twoN;
    if (e < 0)
        e += twoN;
    int sign = 1;
    if (e >= static_cast<int64_t>(N)) {
        sign = -1;
        e -= N;
    }
    if (e == 0 && sign == 1)
        return;
    for (auto& tw : p.GetAllElements()) {
        const uint32_t n      = tw.GetLength();
        const NativeInteger q = tw.GetModulus();
        std::vector<NativeInteger> tmp(n, NativeInteger(0));
        for (uint32_t i = 0; i < n; ++i) {
            const uint32_t dst = (i + static_cast<uint32_t>(e)) % n;
            const bool wrap    = (i + static_cast<uint32_t>(e)) >= n;
            NativeInteger v    = tw[i];
            int s              = sign;
            if (wrap)
                s = -s;
            if (s < 0)
                v = q.ModSub(v, q);
            tmp[dst] = v;
        }
        for (uint32_t i = 0; i < n; ++i)
            tw[i] = tmp[i];
    }
}

DCRTPoly MakeEvalMonomial(const DCRTPoly& like, int64_t power) {
    DCRTPoly p(like.GetParams(), Format::COEFFICIENT, true);
    const uint32_t N   = p.GetRingDimension();
    const int64_t twoN = 2 * static_cast<int64_t>(N);
    int64_t e          = power % twoN;
    if (e < 0)
        e += twoN;
    int sign = 1;
    if (e >= static_cast<int64_t>(N)) {
        sign = -1;
        e -= N;
    }
    for (auto& tw : p.GetAllElements()) {
        if (sign > 0)
            tw[static_cast<uint32_t>(e)] = 1;
        else
            tw[static_cast<uint32_t>(e)] = tw.GetModulus() - 1;
    }
    p.SetFormat(Format::EVALUATION);
    return p;
}

void MultByMonomialInPlace(Ciphertext<DCRTPoly>& ct, int64_t power, bool spectral) {
    auto& els = ct->GetElements();
    if (spectral) {
        if (els.empty())
            return;
        DCRTPoly mono = MakeEvalMonomial(els[0], power);
        for (auto& el : els) {
            el.SetFormat(Format::EVALUATION);
            el *= mono;
        }
    }
    else {
        for (auto& el : els)
            NegacyclicMulXInPlace(el, power);
    }
}

void ScaleByDInvInPlace(Ciphertext<DCRTPoly>& ct, uint32_t d) {
    for (auto& el : ct->GetElements()) {
        el.SetFormat(Format::EVALUATION);
        for (auto& tw : el.GetAllElements()) {
            const NativeInteger q   = tw.GetModulus();
            const NativeInteger inv = NativeInteger(d).ModInverse(q);
            const uint32_t n        = tw.GetLength();
            for (uint32_t i = 0; i < n; ++i)
                tw[i].ModMulEq(inv, q);
        }
    }
}

// In-place DIT radix-2, identical to the even/odd TWEAK recurrence.
// bit-reverse then butterflies with wLen = omega^{d/len}.
void TweakValsInplace(NativeInteger* a, uint32_t d, NativeInteger omega, const NativeInteger& q,
                      const uint32_t* rev = nullptr, const NativeInteger* wLenCached = nullptr) {
    if (d <= 1)
        return;
    if (rev) {
        NativeInteger tmp[4096];
        if (d > 4096)
            OPENFHE_THROW("d > 4096 is not supported in TweakValsInplace stack allocation");
        for (uint32_t i = 0; i < d; ++i)
            tmp[i] = a[rev[i]];
        for (uint32_t i = 0; i < d; ++i)
            a[i] = tmp[i];
    }
    else {
        for (uint32_t i = 1, j = 0; i < d; ++i) {
            uint32_t bit = d >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                std::swap(a[i], a[j]);
        }
    }
    uint32_t stage = 0;
    for (uint32_t len = 2; len <= d; len <<= 1, ++stage) {
        const uint32_t half      = len >> 1;
        const NativeInteger wLen = wLenCached ? wLenCached[stage] : ModExp(omega, d / len, q);
        for (uint32_t start = 0; start < d; start += len) {
            NativeInteger w(1);
            for (uint32_t j = 0; j < half; ++j) {
                const NativeInteger t = w.ModMul(a[start + j + half], q);
                const NativeInteger u = a[start + j];
                a[start + j]          = u.ModAdd(t, q);
                a[start + j + half]   = u.ModSub(t, q);
                w.ModMulEq(wLen, q);
            }
        }
    }
}

void TweakRecursive(std::vector<Ciphertext<DCRTPoly>>& cts, uint32_t k, int sgn, bool spectral) {
    const uint32_t d = static_cast<uint32_t>(cts.size());
    if (d <= 1)
        return;
    std::vector<Ciphertext<DCRTPoly>> even, odd;
    even.reserve(d / 2);
    odd.reserve(d / 2);
    for (uint32_t i = 0; i < d; i += 2) {
        even.push_back(cts[i]);
        odd.push_back(cts[i + 1]);
    }
    TweakRecursive(even, 2 * k, sgn, spectral);
    TweakRecursive(odd, 2 * k, sgn, spectral);
    for (uint32_t j = 0; j < d / 2; ++j) {
        auto twiddled           = odd[j]->Clone();
        const int64_t power     = static_cast<int64_t>(sgn) * 2 * static_cast<int64_t>(k) * static_cast<int64_t>(j);
        MultByMonomialInPlace(twiddled, power, spectral);
        cts[j]         = AddCT(even[j], twiddled);
        cts[j + d / 2] = SubCT(even[j], twiddled);
    }
}

// Eval-domain TWEAK: one monomial NTT, then a length-d DFT at each slot.
// Same math as TweakRecursive(..., spectral=true), without an NTT per butterfly.
// If scaleD != 0, each output slot is also multiplied by scaleD^{-1} (fused d^{-1}).
void FastEvalTweak(std::vector<Ciphertext<DCRTPoly>>& cts, uint32_t k, int sgn, uint32_t scaleD,
                   const DCRTPoly* cachedOmega, const uint32_t* dftRev,
                   const std::vector<std::vector<NativeInteger>>* wlenByTow = nullptr, uint32_t logd = 0) {
    const uint32_t d = static_cast<uint32_t>(cts.size());
    if (d <= 1)
        return;
    SetAllFormat(cts, Format::EVALUATION);
    DCRTPoly monoLocal;
    const DCRTPoly* mono = cachedOmega;
    if (!mono) {
        monoLocal = MakeEvalMonomial(cts[0]->GetElements()[0],
                                     static_cast<int64_t>(sgn) * 2 * static_cast<int64_t>(k));
        mono = &monoLocal;
    }
    const uint32_t nComp = static_cast<uint32_t>(cts[0]->GetElements().size());
    const uint32_t nTow  = cts[0]->GetElements()[0].GetNumOfElements();
    const uint32_t N     = cts[0]->GetElements()[0].GetRingDimension();
    const uint32_t nPl   = nComp * nTow;

    std::vector<NativePoly*> slot(static_cast<size_t>(nPl) * d);
    std::vector<NativeInteger> qv(nPl), dinv(nPl, NativeInteger(1));
    std::vector<const NativePoly*> wv(nPl);
    for (uint32_t comp = 0; comp < nComp; ++comp) {
        for (uint32_t tow = 0; tow < nTow; ++tow) {
            const uint32_t pl = comp * nTow + tow;
            qv[pl]            = cts[0]->GetElements()[comp].GetElementAtIndex(tow).GetModulus();
            wv[pl]            = &mono->GetElementAtIndex(tow);
            if (scaleD > 1)
                dinv[pl] = NativeInteger(scaleD).ModInverse(qv[pl]);
            for (uint32_t t = 0; t < d; ++t)
                slot[static_cast<size_t>(pl) * d + t] =
                    &cts[t]->GetElements()[comp].GetAllElements()[tow];
        }
    }

#pragma omp parallel
    {
        std::vector<NativeInteger> vals(d);
#pragma omp for schedule(static) collapse(2)
        for (int64_t pl = 0; pl < static_cast<int64_t>(nPl); ++pl) {
            for (int64_t s = 0; s < static_cast<int64_t>(N); ++s) {
                const uint32_t sto = static_cast<uint32_t>(s);
                const NativeInteger q = qv[static_cast<uint32_t>(pl)];
                NativePoly** sl       = &slot[static_cast<size_t>(pl) * d];
                for (uint32_t t = 0; t < d; ++t)
                    vals[t] = (*sl[t])[sto];
                const NativeInteger* wcached = nullptr;
                if (wlenByTow && !wlenByTow->empty()) {
                    const uint32_t nTowLoc = static_cast<uint32_t>((*wlenByTow).size());
                    const uint32_t tow     = (nTowLoc > 0) ? (static_cast<uint32_t>(pl) % nTowLoc) : 0;
                    wcached                = (*wlenByTow)[tow].data() + static_cast<size_t>(sto) * logd;
                }
                TweakValsInplace(vals.data(), d, (*wv[static_cast<uint32_t>(pl)])[sto], q, dftRev, wcached);
                if (scaleD > 1) {
                    const NativeInteger inv = dinv[static_cast<uint32_t>(pl)];
                    for (uint32_t t = 0; t < d; ++t)
                        (*sl[t])[sto] = vals[t].ModMul(inv, q);
                }
                else {
                    for (uint32_t t = 0; t < d; ++t)
                        (*sl[t])[sto] = vals[t];
                }
            }
        }
    }
}

void ApplyTWEAK(std::vector<Ciphertext<DCRTPoly>>& cts, uint32_t /*N*/, uint32_t k, int sgn, bool spectral,
                uint32_t scaleD = 0, const DCRTPoly* cachedOmega = nullptr, const uint32_t* dftRev = nullptr,
                const std::vector<std::vector<NativeInteger>>* wlenByTow = nullptr, uint32_t logd = 0) {
    if (spectral) {
        FastEvalTweak(cts, k, sgn, scaleD, cachedOmega, dftRev, wlenByTow, logd);
        return;
    }
    SetAllFormat(cts, Format::COEFFICIENT);
    TweakRecursive(cts, k, sgn, false);
    SetAllFormat(cts, Format::EVALUATION);
}

// Independent streams after TWEAK: hoist all, then KS-core all. Wall-clock of
// each parallel phase is T_hoist / T_keyprod (not the sum of per-stream times).
std::vector<Ciphertext<DCRTPoly>> BatchKeySwitchNoAuto(const CryptoContext<DCRTPoly>& cc,
                                                       const std::vector<Ciphertext<DCRTPoly>>& cts,
                                                       const std::vector<uint32_t>& hFor, CMTTimings* tm) {
    const uint32_t d = static_cast<uint32_t>(cts.size());
    std::vector<Ciphertext<DCRTPoly>> out(d);
    std::vector<std::shared_ptr<std::vector<DCRTPoly>>> digits(d);
    auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(cts[0]->GetKeyTag());

    auto t0 = Clock::now();
#pragma omp parallel for schedule(static)
    for (int64_t ts = 0; ts < static_cast<int64_t>(d); ++ts) {
        if (hFor[static_cast<uint32_t>(ts)] == 1)
            continue;
        auto c = cts[static_cast<uint32_t>(ts)];
        digits[static_cast<uint32_t>(ts)] = cc->EvalFastRotationPrecompute(c);
    }
    if (tm)
        tm->t_hoist_ms += MsSince(t0);

    t0 = Clock::now();
#pragma omp parallel for schedule(static)
    for (int64_t ts64 = 0; ts64 < static_cast<int64_t>(d); ++ts64) {
        const uint32_t ts = static_cast<uint32_t>(ts64);
        const uint32_t h  = hFor[ts];
        if (h == 1) {
            out[ts] = cts[ts]->Clone();
            continue;
        }
        auto it = evalKeyMap.find(h);
        if (it == evalKeyMap.end())
            OPENFHE_THROW("EvalKey for index [" + std::to_string(h) + "] is not found.");
        const auto& cv = cts[ts]->GetElements();
        auto baPtr     = cc->GetScheme()->EvalFastKeySwitchCore(digits[ts], it->second, cv[0].GetParams());
        auto ba        = *baPtr;
        ba[0] += cv[0];
        auto result = cts[ts]->CloneEmpty();
        result->SetElements(std::move(ba));
        out[ts] = result;
    }
    if (tm)
        tm->t_keyprod_ms += MsSince(t0);
    return out;
}

// Same d-way OpenMP as BatchKeySwitchNoAuto. CKL paper loop, just concurrent
// independent EvalAutomorphism calls. Wall-clock of the region is T_KS/T_auto.
std::vector<Ciphertext<DCRTPoly>> BatchCKLAutomorphism(const CryptoContext<DCRTPoly>& cc,
                                                       const std::vector<Ciphertext<DCRTPoly>>& aux,
                                                       const CMTIndexTable& idx, CMTTimings* tm) {
    const uint32_t d = idx.d;
    std::vector<Ciphertext<DCRTPoly>> out(d);
    auto& keymap = cc->GetEvalAutomorphismKeyMap(aux[0]->GetKeyTag());
    auto t0      = Clock::now();
#pragma omp parallel for schedule(static)
    for (int64_t t64 = 0; t64 < static_cast<int64_t>(d); ++t64) {
        const uint32_t t  = static_cast<uint32_t>(t64);
        const uint32_t ts = idx.tstar[t];
        const uint32_t h  = idx.h[t];
        if (h == 1)
            out[t] = aux[ts]->Clone();
        else
            out[t] = cc->EvalAutomorphism(aux[ts], h, keymap);
    }
    if (tm) {
        const double ms = MsSince(t0);
        tm->t_keyprod_ms += ms;
        tm->t_automorphism_ms += ms;
    }
    return out;
}

// Hoisted KS (parallel) then automorphism permute (parallel). Same thread
// layout as FCMT's BatchKeySwitchNoAuto.
std::vector<Ciphertext<DCRTPoly>> BatchHoistedAutomorphism(const CryptoContext<DCRTPoly>& cc,
                                                           const std::vector<Ciphertext<DCRTPoly>>& aux,
                                                           const CMTIndexTable& idx, CMTTimings* tm) {
    const uint32_t d = idx.d;
    auto u           = BatchKeySwitchNoAuto(cc, aux, idx.h_for, tm);
    std::vector<Ciphertext<DCRTPoly>> out(d);
    const uint32_t N = u[0]->GetElements()[0].GetRingDimension();
    auto t0          = Clock::now();
#pragma omp parallel for schedule(static)
    for (int64_t t64 = 0; t64 < static_cast<int64_t>(d); ++t64) {
        const uint32_t t  = static_cast<uint32_t>(t64);
        const uint32_t ts = idx.tstar[t];
        const uint32_t h  = idx.h[t];
        if (h == 1) {
            out[t] = u[ts]->Clone();
            continue;
        }
        std::vector<uint32_t> vec(N);
        PrecomputeAutoMap(N, h, &vec);
        auto result = u[ts]->Clone();
        for (auto& el : result->GetElements()) {
            el.SetFormat(Format::EVALUATION);
            el = el.AutomorphismTransform(h, vec);
        }
        out[t] = result;
    }
    if (tm)
        tm->t_automorphism_ms += MsSince(t0);
    return out;
}

// Shared fused consumer: for each natural ν, v_t = u[t*][storage(π_t(ν))], DFT, phase.
// Fast path gathers through the precomputed π table. Residue packing is only the
// AFCMT_RESIDUE ablation (an extra tensor copy; not used by FCMT / FUSED / BATCHED).
std::vector<Ciphertext<DCRTPoly>> SpectralConsume(const CMTContext& ctx,
                                                  const std::vector<Ciphertext<DCRTPoly>>& u, CMTTimings* tm,
                                                  AFCMTMode mode) {
    const auto& idx  = ctx.idx;
    const uint32_t d = idx.d;
    const uint32_t N = idx.N;
    const uint32_t k = idx.k;
    const int64_t m2k = -2 * static_cast<int64_t>(k);

    auto proto = u[0]->CloneEmpty();
    {
        const auto& src = u[0]->GetElements();
        std::vector<DCRTPoly> els;
        els.reserve(src.size());
        for (const auto& el : src)
            els.emplace_back(el.GetParams(), Format::EVALUATION, true);
        proto->SetElements(std::move(els));
    }
    const uint32_t nComp = static_cast<uint32_t>(u[0]->GetElements().size());
    const uint32_t nTow  = u[0]->GetElements()[0].GetNumOfElements();

    DCRTPoly monoM2kLocal, monoM1Local;
    if (!ctx.haveMonos) {
        monoM2kLocal = MakeEvalMonomial(u[0]->GetElements()[0], m2k);
        monoM1Local  = MakeEvalMonomial(u[0]->GetElements()[0], -1);
    }
    const DCRTPoly& monoM2k = ctx.haveMonos ? ctx.monoM2k : monoM2kLocal;
    const DCRTPoly& monoM1  = ctx.haveMonos ? ctx.monoM1 : monoM1Local;

    std::vector<Ciphertext<DCRTPoly>> out(d);
    out[0] = proto;
    for (uint32_t j = 1; j < d; ++j)
        out[j] = proto->Clone();

    if (tm) {
        for (uint32_t t = 0; t < d; ++t)
            tm->traffic.gather_bytes += EvalFormCiphertextBytes(u[t]);
        tm->traffic.bytes_read += tm->traffic.gather_bytes;
    }

    auto tGather = Clock::now();
    std::vector<std::vector<std::vector<std::vector<NativeInteger>>>> Upack;
    const bool usePack = (mode == AFCMTMode::AFCMT_RESIDUE);
    if (usePack) {
        Upack.assign(d, std::vector<std::vector<std::vector<NativeInteger>>>(
                            nComp, std::vector<std::vector<NativeInteger>>(nTow)));
        uint64_t layoutBytes = 0;
        for (uint32_t ts = 0; ts < d; ++ts) {
            for (uint32_t comp = 0; comp < nComp; ++comp) {
                for (uint32_t tow = 0; tow < nTow; ++tow) {
                    const NativePoly& src = u[ts]->GetElements()[comp].GetElementAtIndex(tow);
                    auto& dst             = Upack[ts][comp][tow];
                    dst.resize(N);
                    for (uint32_t a = 0; a < k; ++a)
                        for (uint32_t b = 0; b < d; ++b) {
                            const uint32_t nuNat = a + b * k;
                            dst[a * d + b]       = src[idx.bitrev[nuNat]];
                        }
                    layoutBytes += static_cast<uint64_t>(N) * sizeof(uint64_t);
                }
            }
        }
        if (tm)
            tm->traffic.residue_layout_bytes += layoutBytes;
    }
    if (tm) {
        const double lay = MsSince(tGather);
        tm->t_layout_ms += lay;
        tm->t_gather_ms += lay;
    }

    auto tDft = Clock::now();
    for (uint32_t comp = 0; comp < nComp; ++comp) {
        for (uint32_t tow = 0; tow < nTow; ++tow) {
            const NativeInteger q  = out[0]->GetElements()[comp].GetElementAtIndex(tow).GetModulus();
            const NativePoly& wM2k = monoM2k.GetElementAtIndex(tow);
            const NativePoly& wM1  = monoM1.GetElementAtIndex(tow);

            std::vector<const NativePoly*> srcSto(d);
            for (uint32_t t = 0; t < d; ++t)
                srcSto[t] = &u[idx.tstar[t]]->GetElements()[comp].GetElementAtIndex(tow);

            std::vector<NativePoly*> dst(d);
            for (uint32_t j = 0; j < d; ++j)
                dst[j] = &out[j]->GetElements()[comp].GetAllElements()[tow];

            const uint32_t* const* pi = nullptr;
            std::vector<const uint32_t*> piRow;
            if (!usePack) {
                piRow.resize(d);
                for (uint32_t t = 0; t < d; ++t)
                    piRow[t] = idx.pi_sto[t].data();
                pi = piRow.data();
            }
            const bool par = true;

#pragma omp parallel if (par)
            {
                std::vector<NativeInteger> vals(d);
#pragma omp for schedule(static)
                for (int64_t nu64 = 0; nu64 < static_cast<int64_t>(N); ++nu64) {
                    const uint32_t nuNat = static_cast<uint32_t>(nu64);
                    const uint32_t sto   = idx.bitrev[nuNat];
                    if (usePack) {
                        const uint32_t a    = nuNat % k;
                        const uint32_t b    = nuNat / k;
                        uint32_t lane       = b;
                        const uint32_t step = (2u * a + 1u) & (d - 1);
                        for (uint32_t t = 0; t < d; ++t) {
                            vals[t] = Upack[idx.tstar[t]][comp][tow][a * d + lane];
                            lane    = (lane + step) & (d - 1);
                        }
                    }
                    else {
                        for (uint32_t t = 0; t < d; ++t)
                            vals[t] = (*srcSto[t])[pi[t][nuNat]];
                    }
                    const NativeInteger* wcached = nullptr;
                    if (ctx.haveWlen && tow < ctx.wlenInv.size() && ctx.logd > 0)
                        wcached = ctx.wlenInv[tow].data() + static_cast<size_t>(sto) * ctx.logd;
                    TweakValsInplace(vals.data(), d, wM2k[sto], q,
                                     idx.dft_rev.empty() ? nullptr : idx.dft_rev.data(), wcached);
                    NativeInteger phase(1);
                    const NativeInteger w1 = wM1[sto];
                    for (uint32_t j = 0; j < d; ++j) {
                        NativeInteger z = vals[j];
                        z.ModMulEq(phase, q);
                        (*dst[j])[sto] = z;
                        phase.ModMulEq(w1, q);
                    }
                }
            }
        }
    }
    if (tm) {
        tm->t_dft_ms += MsSince(tDft);
        tm->t_spectral_tweak_ms += tm->t_dft_ms;
    }
    return out;
}

void EncodeColumnPoly(DCRTPoly& p, const std::vector<std::vector<int64_t>>& M, uint32_t j, int64_t scale) {
    p.SetFormat(Format::COEFFICIENT);
    const uint32_t d = static_cast<uint32_t>(M.size());
    for (auto& tw : p.GetAllElements()) {
        const NativeInteger q = tw.GetModulus();
        const uint32_t n      = tw.GetLength();
        for (uint32_t i = 0; i < n; ++i)
            tw[i] = 0;
        for (uint32_t i = 0; i < d && i < n; ++i) {
            int64_t v = M[i][j] * scale;
            if (v >= 0)
                tw[i] = NativeInteger(static_cast<uint64_t>(v)).Mod(q);
            else {
                NativeInteger pos(static_cast<uint64_t>(-v));
                tw[i] = q.ModSub(pos.Mod(q), q);
            }
        }
    }
}

void EncodeRowPoly(DCRTPoly& p, const std::vector<std::vector<int64_t>>& M, uint32_t i, int64_t scale) {
    p.SetFormat(Format::COEFFICIENT);
    const uint32_t d = static_cast<uint32_t>(M[0].size());
    for (auto& tw : p.GetAllElements()) {
        const NativeInteger q = tw.GetModulus();
        const uint32_t n      = tw.GetLength();
        for (uint32_t t = 0; t < n; ++t)
            tw[t] = 0;
        for (uint32_t j = 0; j < d && j < n; ++j) {
            int64_t v = M[i][j] * scale;
            if (v >= 0)
                tw[j] = NativeInteger(static_cast<uint64_t>(v)).Mod(q);
            else {
                NativeInteger pos(static_cast<uint64_t>(-v));
                tw[j] = q.ModSub(pos.Mod(q), q);
            }
        }
    }
}

std::vector<Ciphertext<DCRTPoly>> EncryptPolys(const CMTContext& ctx, const std::vector<DCRTPoly>& msgs) {
    const uint32_t batch = ctx.cc->GetEncodingParams()->GetBatchSize();
    auto ptxt            = ctx.cc->MakeCKKSPackedPlaintext(std::vector<double>(batch, 0.0));
    std::vector<Ciphertext<DCRTPoly>> out;
    out.reserve(msgs.size());
    for (const auto& msg : msgs) {
        auto ct      = ctx.cc->Encrypt(ctx.pk, ptxt);
        DCRTPoly add = msg;
        add.SetFormat(Format::EVALUATION);
        const uint32_t want = ct->GetElements()[0].GetNumOfElements();
        const uint32_t have = add.GetNumOfElements();
        if (have > want)
            add.DropLastElements(have - want);
        ct->GetElements()[0] += add;
        out.push_back(ct);
    }
    return out;
}

void BitReverseInPlace(std::vector<NativeInteger>& a) {
    const uint32_t n = static_cast<uint32_t>(a.size());
    for (uint32_t i = 1, j = 0; i < n; ++i) {
        uint32_t bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            std::swap(a[i], a[j]);
    }
}

void CyclicNTT(std::vector<NativeInteger>& a, NativeInteger omega, const NativeInteger& q, bool inverse) {
    const uint32_t n = static_cast<uint32_t>(a.size());
    BitReverseInPlace(a);
    for (uint32_t len = 2; len <= n; len <<= 1) {
        NativeInteger wlen = ModExp(omega, n / len, q);
        if (inverse)
            wlen = wlen.ModInverse(q);
        for (uint32_t i = 0; i < n; i += len) {
            NativeInteger w(1);
            for (uint32_t j = 0; j < len / 2; ++j) {
                const NativeInteger u = a[i + j];
                const NativeInteger v = a[i + j + len / 2].ModMul(w, q);
                a[i + j]              = u.ModAdd(v, q);
                a[i + j + len / 2]    = u.ModSub(v, q);
                w.ModMulEq(wlen, q);
            }
        }
    }
    if (inverse) {
        const NativeInteger ninv = NativeInteger(n).ModInverse(q);
        for (auto& x : a)
            x.ModMulEq(ninv, q);
    }
}

void NegacyclicNTT(std::vector<NativeInteger>& a, NativeInteger psi, const NativeInteger& q, bool inverse) {
    const uint32_t kk         = static_cast<uint32_t>(a.size());
    const NativeInteger omega = psi.ModMul(psi, q);
    if (!inverse) {
        NativeInteger w(1);
        for (uint32_t i = 0; i < kk; ++i) {
            a[i].ModMulEq(w, q);
            w.ModMulEq(psi, q);
        }
        CyclicNTT(a, omega, q, false);
    }
    else {
        CyclicNTT(a, omega, q, true);
        const NativeInteger psiInv = psi.ModInverse(q);
        NativeInteger w(1);
        for (uint32_t i = 0; i < kk; ++i) {
            a[i].ModMulEq(w, q);
            w.ModMulEq(psiInv, q);
        }
    }
}

void ExtractModuleCoords(const DCRTPoly& poly, uint32_t d, uint32_t k, uint32_t tow,
                         std::vector<std::vector<NativeInteger>>& coords) {
    const NativePoly& tw = poly.GetElementAtIndex(tow);
    coords.assign(d, std::vector<NativeInteger>(k));
    for (uint32_t i = 0; i < d; ++i)
        for (uint32_t s = 0; s < k; ++s)
            coords[i][s] = tw[s * d + i];
}

void WriteModuleCoords(DCRTPoly& poly, uint32_t d, uint32_t k, uint32_t tow,
                       const std::vector<std::vector<NativeInteger>>& coords) {
    NativePoly& tw = poly.GetAllElements()[tow];
    for (uint32_t i = 0; i < d; ++i)
        for (uint32_t s = 0; s < k; ++s)
            tw[s * d + i] = coords[i][s];
}

[[maybe_unused]] std::vector<DCRTPoly> PPMMColumns(const std::vector<DCRTPoly>& B, const std::vector<DCRTPoly>& U,
                                                   uint32_t d, uint32_t k) {
    std::vector<DCRTPoly> C(d);
    for (uint32_t j = 0; j < d; ++j) {
        C[j] = B[0];
        C[j].SetFormat(Format::COEFFICIENT);
        for (auto& tw : C[j].GetAllElements()) {
            const uint32_t n = tw.GetLength();
            for (uint32_t i = 0; i < n; ++i)
                tw[i] = 0;
        }
    }
    std::vector<DCRTPoly> Bc = B, Uc = U;
    for (auto& p : Bc)
        p.SetFormat(Format::COEFFICIENT);
    for (auto& p : Uc)
        p.SetFormat(Format::COEFFICIENT);

    const uint32_t nTow = Bc[0].GetNumOfElements();
    for (uint32_t tow = 0; tow < nTow; ++tow) {
        const NativeInteger q    = Bc[0].GetElementAtIndex(tow).GetModulus();
        const NativeInteger root = Bc[0].GetElementAtIndex(tow).GetRootOfUnity();
        const uint32_t Nring     = Bc[0].GetRingDimension();
        const NativeInteger psi  = ModExp(root, Nring / k, q);

        std::vector<std::vector<std::vector<NativeInteger>>> Bhat(d, std::vector<std::vector<NativeInteger>>(d));
        std::vector<std::vector<std::vector<NativeInteger>>> Uhat(d, std::vector<std::vector<NativeInteger>>(d));
        for (uint32_t col = 0; col < d; ++col) {
            std::vector<std::vector<NativeInteger>> coords;
            ExtractModuleCoords(Bc[col], d, k, tow, coords);
            for (uint32_t i = 0; i < d; ++i) {
                Bhat[i][col] = coords[i];
                NegacyclicNTT(Bhat[i][col], psi, q, false);
            }
            ExtractModuleCoords(Uc[col], d, k, tow, coords);
            for (uint32_t i = 0; i < d; ++i) {
                Uhat[i][col] = coords[i];
                NegacyclicNTT(Uhat[i][col], psi, q, false);
            }
        }
        std::vector<std::vector<std::vector<NativeInteger>>> Chat(
            d, std::vector<std::vector<NativeInteger>>(d, std::vector<NativeInteger>(k, NativeInteger(0))));
        for (uint32_t freq = 0; freq < k; ++freq)
            for (uint32_t i = 0; i < d; ++i)
                for (uint32_t l = 0; l < d; ++l) {
                    NativeInteger acc(0);
                    for (uint32_t j = 0; j < d; ++j)
                        acc.ModAddEq(Bhat[i][j][freq].ModMul(Uhat[j][l][freq], q), q);
                    Chat[i][l][freq] = acc;
                }
        for (uint32_t l = 0; l < d; ++l) {
            std::vector<std::vector<NativeInteger>> Ccol(d);
            for (uint32_t i = 0; i < d; ++i) {
                Ccol[i] = Chat[i][l];
                NegacyclicNTT(Ccol[i], psi, q, true);
            }
            WriteModuleCoords(C[l], d, k, tow, Ccol);
        }
    }
    for (auto& p : C)
        p.SetFormat(Format::EVALUATION);
    return C;
}

[[maybe_unused]] std::vector<DCRTPoly> ModuleTranspose(const std::vector<DCRTPoly>& in, uint32_t d, uint32_t k) {
    std::vector<DCRTPoly> out(d);
    std::vector<DCRTPoly> tmp = in;
    for (auto& p : tmp)
        p.SetFormat(Format::COEFFICIENT);
    for (uint32_t j = 0; j < d; ++j) {
        out[j] = tmp[0];
        for (auto& tw : out[j].GetAllElements()) {
            const uint32_t n = tw.GetLength();
            for (uint32_t i = 0; i < n; ++i)
                tw[i] = 0;
        }
    }
    const uint32_t nTow = tmp[0].GetNumOfElements();
    for (uint32_t tow = 0; tow < nTow; ++tow) {
        std::vector<std::vector<std::vector<NativeInteger>>> cols(d);
        for (uint32_t j = 0; j < d; ++j)
            ExtractModuleCoords(tmp[j], d, k, tow, cols[j]);
        for (uint32_t j = 0; j < d; ++j) {
            std::vector<std::vector<NativeInteger>> ncol(d, std::vector<NativeInteger>(k));
            for (uint32_t i = 0; i < d; ++i)
                ncol[i] = cols[i][j];
            WriteModuleCoords(out[j], d, k, tow, ncol);
        }
    }
    for (auto& p : out)
        p.SetFormat(Format::EVALUATION);
    return out;
}

Ciphertext<DCRTPoly> PackCipher(const Ciphertext<DCRTPoly>& proto, DCRTPoly c0, DCRTPoly c1) {
    auto ct = proto->Clone();
    c0.SetFormat(Format::EVALUATION);
    c1.SetFormat(Format::EVALUATION);
    ct->SetElements({std::move(c0), std::move(c1)});
    return ct;
}

// ----- Proposed CCMM kernels (same PPMM arithmetic, fused + OpenMP) -----

struct NegacyclicPlan {
    NativeInteger q, ninv;
    std::vector<NativeInteger> psiPow, psiInvPow;
    std::vector<NativeInteger> wlenFwd, wlenInv;
};

NegacyclicPlan MakeNegacyclicPlan(const NativeInteger& psi, const NativeInteger& q, uint32_t k) {
    NegacyclicPlan p;
    p.q                  = q;
    p.ninv               = NativeInteger(k).ModInverse(q);
    const NativeInteger psiInv = psi.ModInverse(q);
    const NativeInteger omega  = psi.ModMul(psi, q);
    const NativeInteger omegaInv = omega.ModInverse(q);
    p.psiPow.resize(k);
    p.psiInvPow.resize(k);
    NativeInteger w(1), wi(1);
    for (uint32_t i = 0; i < k; ++i) {
        p.psiPow[i]    = w;
        p.psiInvPow[i] = wi;
        w.ModMulEq(psi, q);
        wi.ModMulEq(psiInv, q);
    }
    uint32_t stages = 0;
    for (uint32_t len = 2; len <= k; len <<= 1)
        ++stages;
    p.wlenFwd.resize(stages);
    p.wlenInv.resize(stages);
    uint32_t s = 0;
    for (uint32_t len = 2; len <= k; len <<= 1, ++s) {
        p.wlenFwd[s] = ModExp(omega, k / len, q);
        p.wlenInv[s] = ModExp(omegaInv, k / len, q);
    }
    return p;
}

void CyclicNTTFast(NativeInteger* a, uint32_t n, const NativeInteger& q, const std::vector<NativeInteger>& wlen,
                   bool inverse, const NativeInteger* ninv) {
    for (uint32_t i = 1, j = 0; i < n; ++i) {
        uint32_t bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            std::swap(a[i], a[j]);
    }
    uint32_t s = 0;
    for (uint32_t len = 2; len <= n; len <<= 1, ++s) {
        const NativeInteger wstep = wlen[s];
        for (uint32_t i = 0; i < n; i += len) {
            NativeInteger w(1);
            for (uint32_t j = 0; j < len / 2; ++j) {
                const NativeInteger u = a[i + j];
                const NativeInteger v = a[i + j + len / 2].ModMul(w, q);
                a[i + j]              = u.ModAdd(v, q);
                a[i + j + len / 2]    = u.ModSub(v, q);
                w.ModMulEq(wstep, q);
            }
        }
    }
    if (inverse && ninv) {
        for (uint32_t i = 0; i < n; ++i)
            a[i].ModMulEq(*ninv, q);
    }
}

void NegacyclicNTTFast(NativeInteger* a, uint32_t k, const NegacyclicPlan& p, bool inverse) {
    const NativeInteger& q = p.q;
    if (!inverse) {
        for (uint32_t i = 0; i < k; ++i)
            a[i].ModMulEq(p.psiPow[i], q);
        CyclicNTTFast(a, k, q, p.wlenFwd, false, nullptr);
    }
    else {
        CyclicNTTFast(a, k, q, p.wlenInv, true, &p.ninv);
        for (uint32_t i = 0; i < k; ++i)
            a[i].ModMulEq(p.psiInvPow[i], q);
    }
}

// Frequency-major: hat[f * d * d + i * d + j] so each freq is a dense d×d block.
inline size_t HatFij(uint32_t f, uint32_t i, uint32_t j, uint32_t d) {
    return (static_cast<size_t>(f) * d + i) * d + j;
}

void ExtractAndNTT(const std::vector<DCRTPoly>& cols, uint32_t d, uint32_t k, uint32_t tow,
                   const NegacyclicPlan& plan, std::vector<NativeInteger>& hat) {
    hat.resize(static_cast<size_t>(d) * d * k);
#pragma omp parallel
    {
        std::vector<NativeInteger> tmp(k);
#pragma omp for collapse(2) schedule(static)
        for (int64_t col = 0; col < static_cast<int64_t>(d); ++col) {
            for (int64_t i = 0; i < static_cast<int64_t>(d); ++i) {
                const NativePoly& tw = cols[static_cast<uint32_t>(col)].GetElementAtIndex(tow);
                const uint32_t ii    = static_cast<uint32_t>(i);
                const uint32_t cj    = static_cast<uint32_t>(col);
                for (uint32_t s = 0; s < k; ++s)
                    tmp[s] = tw[s * d + ii];
                NegacyclicNTTFast(tmp.data(), k, plan, false);
                for (uint32_t s = 0; s < k; ++s)
                    hat[HatFij(s, ii, cj, d)] = tmp[s];
            }
        }
    }
}

void INTTAndWrite(std::vector<DCRTPoly>& out, uint32_t d, uint32_t k, uint32_t tow, const NegacyclicPlan& plan,
                  const std::vector<NativeInteger>& hat) {
#pragma omp parallel
    {
        std::vector<NativeInteger> tmp(k);
#pragma omp for collapse(2) schedule(static)
        for (int64_t l = 0; l < static_cast<int64_t>(d); ++l) {
            for (int64_t i = 0; i < static_cast<int64_t>(d); ++i) {
                const uint32_t ii = static_cast<uint32_t>(i);
                const uint32_t ll = static_cast<uint32_t>(l);
                for (uint32_t s = 0; s < k; ++s)
                    tmp[s] = hat[HatFij(s, ii, ll, d)];
                NegacyclicNTTFast(tmp.data(), k, plan, true);
                NativePoly& tw = out[ll].GetAllElements()[tow];
                for (uint32_t s = 0; s < k; ++s)
                    tw[s * d + ii] = tmp[s];
            }
        }
    }
}

// Four d×d products at every frequency: {LB,LA} × {RB,RA}, one memory pass.
void FreqMatMul4(const std::vector<NativeInteger>& LBh, const std::vector<NativeInteger>& LAh,
                 const std::vector<NativeInteger>& RBh, const std::vector<NativeInteger>& RAh,
                 std::vector<NativeInteger>& H00, std::vector<NativeInteger>& H01, std::vector<NativeInteger>& H10,
                 std::vector<NativeInteger>& H11, uint32_t d, uint32_t k, const NativeInteger& q) {
    const size_t n = static_cast<size_t>(d) * d * k;
    H00.resize(n);
    H01.resize(n);
    H10.resize(n);
    H11.resize(n);
    const size_t blk = static_cast<size_t>(d) * d;
#pragma omp parallel for schedule(static)
    for (int64_t freq = 0; freq < static_cast<int64_t>(k); ++freq) {
        const size_t base = static_cast<size_t>(freq) * blk;
        const NativeInteger* LB = LBh.data() + base;
        const NativeInteger* LA = LAh.data() + base;
        const NativeInteger* RB = RBh.data() + base;
        const NativeInteger* RA = RAh.data() + base;
        NativeInteger* C00      = H00.data() + base;
        NativeInteger* C01      = H01.data() + base;
        NativeInteger* C10      = H10.data() + base;
        NativeInteger* C11      = H11.data() + base;
        std::vector<uint64_t> RB_T(d * d);
        std::vector<uint64_t> RA_T(d * d);
        std::vector<uint64_t> LB_I(d * d);
        std::vector<uint64_t> LA_I(d * d);
        for (uint32_t j = 0; j < d; ++j) {
            for (uint32_t l = 0; l < d; ++l) {
                RB_T[l * d + j] = RB[j * d + l].ConvertToInt();
                RA_T[l * d + j] = RA[j * d + l].ConvertToInt();
                LB_I[j * d + l] = LB[j * d + l].ConvertToInt();
                LA_I[j * d + l] = LA[j * d + l].ConvertToInt();
            }
        }
        uint64_t q_val = q.ConvertToInt();
        for (uint32_t i = 0; i < d; ++i) {
            for (uint32_t l = 0; l < d; ++l) {
                unsigned __int128 a00 = 0, a01 = 0, a10 = 0, a11 = 0;
                const uint64_t* lb_ptr = &LB_I[i * d];
                const uint64_t* la_ptr = &LA_I[i * d];
                const uint64_t* rb_ptr = &RB_T[l * d];
                const uint64_t* ra_ptr = &RA_T[l * d];
                
                for (uint32_t j = 0; j < d; ++j) {
                    uint64_t lb = lb_ptr[j];
                    uint64_t la = la_ptr[j];
                    uint64_t rb = rb_ptr[j];
                    uint64_t ra = ra_ptr[j];
                    a00 += (unsigned __int128)lb * rb;
                    a01 += (unsigned __int128)lb * ra;
                    a10 += (unsigned __int128)la * rb;
                    a11 += (unsigned __int128)la * ra;
                }
                const size_t o = static_cast<size_t>(i) * d + l;
                C00[o] = NativeInteger(a00 % q_val);
                C01[o] = NativeInteger(a01 % q_val);
                C10[o] = NativeInteger(a10 % q_val);
                C11[o] = NativeInteger(a11 % q_val);
            }
        }
    }
}

void ZeroCoeffPoly(DCRTPoly& p) {
    p.SetFormat(Format::COEFFICIENT);
    for (auto& tw : p.GetAllElements()) {
        const uint32_t n = tw.GetLength();
        for (uint32_t i = 0; i < n; ++i)
            tw[i] = 0;
    }
}

// Four R_k products from {LB,LA} x {RB,RA} with each operand NTTed once.
void FusedPPMM4(std::vector<DCRTPoly> LBc, std::vector<DCRTPoly> LAc, std::vector<DCRTPoly> RBc,
                std::vector<DCRTPoly> RAc, std::vector<DCRTPoly>& C00, std::vector<DCRTPoly>& C01,
                std::vector<DCRTPoly>& C10, std::vector<DCRTPoly>& C11, uint32_t d, uint32_t k) {
    auto prep = [&](std::vector<DCRTPoly>& C, const DCRTPoly& proto) {
        C.resize(d);
        for (uint32_t j = 0; j < d; ++j) {
            C[j] = proto;
            ZeroCoeffPoly(C[j]);
        }
    };
#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < static_cast<int64_t>(4 * d); ++i) {
        const uint32_t grp = static_cast<uint32_t>(i) / d;
        const uint32_t j   = static_cast<uint32_t>(i) % d;
        if (grp == 0)
            LBc[j].SetFormat(Format::COEFFICIENT);
        else if (grp == 1)
            LAc[j].SetFormat(Format::COEFFICIENT);
        else if (grp == 2)
            RBc[j].SetFormat(Format::COEFFICIENT);
        else
            RAc[j].SetFormat(Format::COEFFICIENT);
    }
    prep(C00, LBc[0]);
    prep(C01, LBc[0]);
    prep(C10, LBc[0]);
    prep(C11, LBc[0]);

    const uint32_t nTow = LBc[0].GetNumOfElements();
    std::vector<NativeInteger> LBh, LAh, RBh, RAh, H00, H01, H10, H11;
    for (uint32_t tow = 0; tow < nTow; ++tow) {
        const NativeInteger q    = LBc[0].GetElementAtIndex(tow).GetModulus();
        const NativeInteger root = LBc[0].GetElementAtIndex(tow).GetRootOfUnity();
        const uint32_t Nring     = LBc[0].GetRingDimension();
        const NativeInteger psi  = ModExp(root, Nring / k, q);
        const auto plan          = MakeNegacyclicPlan(psi, q, k);

        ExtractAndNTT(LBc, d, k, tow, plan, LBh);
        ExtractAndNTT(LAc, d, k, tow, plan, LAh);
        ExtractAndNTT(RBc, d, k, tow, plan, RBh);
        ExtractAndNTT(RAc, d, k, tow, plan, RAh);
        FreqMatMul4(LBh, LAh, RBh, RAh, H00, H01, H10, H11, d, k, q);
        INTTAndWrite(C00, d, k, tow, plan, H00);
        INTTAndWrite(C01, d, k, tow, plan, H01);
        INTTAndWrite(C10, d, k, tow, plan, H10);
        INTTAndWrite(C11, d, k, tow, plan, H11);
    }
#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < static_cast<int64_t>(4 * d); ++i) {
        const uint32_t grp = static_cast<uint32_t>(i) / d;
        const uint32_t j   = static_cast<uint32_t>(i) % d;
        if (grp == 0)
            C00[j].SetFormat(Format::EVALUATION);
        else if (grp == 1)
            C01[j].SetFormat(Format::EVALUATION);
        else if (grp == 2)
            C10[j].SetFormat(Format::EVALUATION);
        else
            C11[j].SetFormat(Format::EVALUATION);
    }
}

std::vector<DCRTPoly> FastModuleTranspose(const std::vector<DCRTPoly>& in, uint32_t d, uint32_t k) {
    std::vector<DCRTPoly> tmp = in;
#pragma omp parallel for schedule(static)
    for (int64_t j = 0; j < static_cast<int64_t>(d); ++j)
        tmp[static_cast<uint32_t>(j)].SetFormat(Format::COEFFICIENT);

    std::vector<DCRTPoly> out(d);
    for (uint32_t j = 0; j < d; ++j) {
        out[j] = tmp[0];
        ZeroCoeffPoly(out[j]);
    }
    const uint32_t nTow = tmp[0].GetNumOfElements();
#pragma omp parallel for collapse(2) schedule(static)
    for (int64_t tow = 0; tow < static_cast<int64_t>(nTow); ++tow) {
        for (int64_t j = 0; j < static_cast<int64_t>(d); ++j) {
            NativePoly& dst = out[static_cast<uint32_t>(j)].GetAllElements()[static_cast<uint32_t>(tow)];
            for (uint32_t s = 0; s < k; ++s)
                for (uint32_t i = 0; i < d; ++i)
                    dst[s * d + i] = tmp[i].GetElementAtIndex(static_cast<uint32_t>(tow))[s * d + static_cast<uint32_t>(j)];
        }
    }
#pragma omp parallel for schedule(static)
    for (int64_t j = 0; j < static_cast<int64_t>(d); ++j)
        out[static_cast<uint32_t>(j)].SetFormat(Format::EVALUATION);
    return out;
}

}  // namespace

uint32_t PiNatural(uint32_t nu, uint32_t t, uint32_t N, uint32_t k) {
    const uint64_t add = static_cast<uint64_t>(k) * static_cast<uint64_t>(t) * (2ull * nu + 1ull);
    return static_cast<uint32_t>((static_cast<uint64_t>(nu) + add) % static_cast<uint64_t>(N));
}

uint32_t LaneAffineGeneral(uint32_t a, uint32_t b, uint32_t t, uint32_t k, uint32_t d) {
    const uint32_t mask  = d - 1;
    const uint32_t alpha = (1u + 2u * k * t) & mask;
    const uint32_t beta  = (t * (2u * a + 1u)) & mask;
    return (alpha * b + beta) & mask;
}

uint32_t LaneAffineSimplified(uint32_t a, uint32_t b, uint32_t t, uint32_t d) {
    const uint32_t mask = d - 1;
    return (b + t * ((2u * a + 1u) & mask)) & mask;
}

uint32_t LaneAffine(uint32_t a, uint32_t b, uint32_t t, uint32_t k, uint32_t d) {
    const uint32_t N = k * d;
    if (!CKLConditionSatisfied(N, d))
        OPENFHE_THROW("LaneAffine simplified form requires d^2 | 2N (equivalently d | 2k)");
    return LaneAffineSimplified(a, b, t, d);
}

uint32_t PiAffine(uint32_t nu, uint32_t t, uint32_t N, uint32_t k, uint32_t d) {
    (void)N;
    const uint32_t a    = nu % k;
    const uint32_t b    = nu / k;
    const uint32_t lane = LaneAffine(a, b, t, k, d);
    return a + k * lane;
}

CMTIndexTable BuildCMTIndexTable(uint32_t N, uint32_t d) {
    if (!CKLConditionSatisfied(N, d))
        OPENFHE_THROW("CKL condition d^2 | 2N is not satisfied");
    CMTIndexTable idx;
    idx.N    = N;
    idx.d    = d;
    idx.k    = N / d;
    idx.logN = 0;
    while ((1u << idx.logN) < N)
        ++idx.logN;
    idx.h.resize(d);
    idx.tstar.resize(d);
    idx.hstar.resize(d);
    idx.bitrev.resize(N);
    idx.alpha.resize(d);
    idx.beta.assign(idx.k, std::vector<uint32_t>(d));
    idx.lane_map.assign(idx.k, std::vector<std::vector<uint32_t>>(d, std::vector<uint32_t>(d)));
    idx.pi_nat.assign(d, std::vector<uint32_t>(N));
    idx.pi_sto.assign(d, std::vector<uint32_t>(N));

    const uint32_t twoN = 2 * N;
    const uint32_t mask = d - 1;
    for (uint32_t nu = 0; nu < N; ++nu)
        idx.bitrev[nu] = ReverseBits(nu, idx.logN);

    uint32_t logd = 0;
    while ((1u << logd) < d)
        ++logd;
    idx.dft_rev.resize(d);
    for (uint32_t i = 0; i < d; ++i)
        idx.dft_rev[i] = ReverseBits(i, logd);

    for (uint32_t t = 0; t < d; ++t) {
        idx.h[t]     = 2 * idx.k * t + 1;
        idx.hstar[t] = ModInverseU32(idx.h[t], twoN);
        idx.tstar[t] = (idx.hstar[t] - 1) / (2 * idx.k);
        idx.alpha[t] = 1u;  // under d|2k, 1+2kt ≡ 1 (mod d)
        for (uint32_t a = 0; a < idx.k; ++a) {
            idx.beta[a][t] = (t * ((2u * a + 1u) & mask)) & mask;
            for (uint32_t b = 0; b < d; ++b)
                idx.lane_map[a][t][b] = (b + idx.beta[a][t]) & mask;
        }
        for (uint32_t nu = 0; nu < N; ++nu) {
            idx.pi_nat[t][nu] = PiNatural(nu, t, N, idx.k);
            idx.pi_sto[t][nu] = idx.bitrev[idx.pi_nat[t][nu]];
        }
    }
    idx.h_for.assign(d, 1u);
    for (uint32_t t = 0; t < d; ++t)
        idx.h_for[idx.tstar[t]] = idx.h[t];
    return idx;
}

std::vector<uint32_t> CMTGaloisIndices(const CMTIndexTable& idx) {
    std::vector<uint32_t> out;
    for (uint32_t t = 0; t < idx.d; ++t)
        if (idx.h[t] != 1)
            out.push_back(idx.h[t]);
    return out;
}

void GenerateCMTEvalKeys(const CMTContext& ctx) {
    auto indices = CMTGaloisIndices(ctx.idx);
    if (!indices.empty())
        ctx.cc->EvalAutomorphismKeyGen(ctx.sk, indices);
    ctx.cc->EvalMultKeyGen(ctx.sk);
}

ClassCountReport CountTransformClasses(uint32_t N, uint32_t d) {
    ClassCountReport r;
    if (!CKLConditionSatisfied(N, d) || d < 2) {
        r.gamma_ok = false;
        r.step_ok  = false;
        return r;
    }
    const uint32_t k = N / d;
    r.expected_gamma       = d / 2;
    r.expected_gamma_count = 2 * k;       // N / (d/2)
    r.expected_step        = d / 2;
    r.expected_step_count  = (2 * k) / d;  // k / (d/2)

    // Gamma class of natural ν is ν mod (d/2) under γ_{ν+d/2}=γ_ν.
    std::vector<uint32_t> gcnt(d / 2, 0);
    for (uint32_t nu = 0; nu < N; ++nu)
        gcnt[nu % (d / 2)]++;
    r.unique_gamma     = 0;
    r.min_gamma_count  = N;
    r.max_gamma_count  = 0;
    for (uint32_t c : gcnt) {
        if (c)
            ++r.unique_gamma;
        r.min_gamma_count = std::min(r.min_gamma_count, c);
        r.max_gamma_count = std::max(r.max_gamma_count, c);
    }
    r.gamma_ok = (r.unique_gamma == r.expected_gamma && r.min_gamma_count == r.expected_gamma_count &&
                  r.max_gamma_count == r.expected_gamma_count);

    std::vector<uint32_t> scnt(d, 0);
    const uint32_t mask = d - 1;
    for (uint32_t a = 0; a < k; ++a)
        scnt[(2u * a + 1u) & mask]++;
    r.unique_step    = 0;
    r.min_step_count = k;
    r.max_step_count = 0;
    for (uint32_t a = 0; a < d; ++a) {
        if (scnt[a] == 0)
            continue;
        ++r.unique_step;
        r.min_step_count = std::min(r.min_step_count, scnt[a]);
        r.max_step_count = std::max(r.max_step_count, scnt[a]);
    }
    r.step_ok = (r.unique_step == r.expected_step && r.min_step_count == r.expected_step_count &&
                 r.max_step_count == r.expected_step_count);
    return r;
}

namespace {
IndexKernelStats FinishIndex(uint64_t checksum, uint64_t nops, double ms) {
    IndexKernelStats s;
    s.checksum       = checksum;
    const double sec = ms / 1000.0;
    s.indices_per_sec = (sec > 0) ? static_cast<double>(nops) / sec : 0;
    s.ns_per_index    = (nops > 0) ? (ms * 1e6) / static_cast<double>(nops) : 0;
    return s;
}
}  // namespace

IndexKernelStats BenchIndexDirect(uint32_t N, uint32_t d, uint32_t iters) {
    const uint32_t k = N / d;
    uint64_t acc     = 0;
    const auto t0    = Clock::now();
    for (uint32_t it = 0; it < iters; ++it)
        for (uint32_t t = 0; t < d; ++t)
            for (uint32_t nu = 0; nu < N; ++nu)
                acc += PiNatural(nu, t, N, k);
    const uint64_t nops = static_cast<uint64_t>(iters) * d * N;
    return FinishIndex(acc, nops, MsSince(t0));
}

IndexKernelStats BenchIndexLookup(uint32_t N, uint32_t d, uint32_t iters) {
    const uint32_t k = N / d;
    std::vector<std::vector<uint32_t>> tab(d, std::vector<uint32_t>(N));
    for (uint32_t t = 0; t < d; ++t)
        for (uint32_t nu = 0; nu < N; ++nu)
            tab[t][nu] = PiNatural(nu, t, N, k);
    uint64_t acc  = 0;
    const auto t0 = Clock::now();
    for (uint32_t it = 0; it < iters; ++it)
        for (uint32_t t = 0; t < d; ++t)
            for (uint32_t nu = 0; nu < N; ++nu)
                acc += tab[t][nu];
    const uint64_t nops = static_cast<uint64_t>(iters) * d * N;
    return FinishIndex(acc, nops, MsSince(t0));
}

IndexKernelStats BenchIndexRecurrence(uint32_t N, uint32_t d, uint32_t iters) {
    const uint32_t k    = N / d;
    const uint32_t mask = d - 1;
    uint64_t acc        = 0;
    const auto t0       = Clock::now();
    for (uint32_t it = 0; it < iters; ++it) {
        for (uint32_t a = 0; a < k; ++a) {
            const uint32_t step = (2u * a + 1u) & mask;
            for (uint32_t b0 = 0; b0 < d; ++b0) {
                uint32_t lane = b0;
                for (uint32_t t = 0; t < d; ++t) {
                    acc += static_cast<uint64_t>(a + k * lane);
                    lane = (lane + step) & mask;
                }
            }
        }
    }
    const uint64_t nops = static_cast<uint64_t>(iters) * d * N;
    return FinishIndex(acc, nops, MsSince(t0));
}

static CMTContext FinishCMTContext(CMTContext ctx, uint32_t N, uint32_t d, uint32_t scaleModSize) {
    ctx.cc->Enable(PKE);
    ctx.cc->Enable(KEYSWITCH);
    ctx.cc->Enable(LEVELEDSHE);
    ctx.cc->Enable(ADVANCEDSHE);
    auto keys         = ctx.cc->KeyGen();
    ctx.pk            = keys.publicKey;
    ctx.sk            = keys.secretKey;
    ctx.idx           = BuildCMTIndexTable(N, d);
    ctx.encodeScale   = 1LL << std::min<int>(40, static_cast<int>(scaleModSize) - 10);
    if (ctx.encodeScale < 1024)
        ctx.encodeScale = 1024;
    GenerateCMTEvalKeys(ctx);

    auto proto = ctx.cc->Encrypt(
        ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(std::vector<double>(ctx.cc->GetEncodingParams()->GetBatchSize(), 0.0)));
    const DCRTPoly& like = proto->GetElements()[0];
    ctx.evalXi.resize(d);
    ctx.evalXi[0] = MakeEvalMonomial(like, 0);
    for (uint32_t i = 1; i < d; ++i)
        ctx.evalXi[i] = MakeEvalMonomial(like, static_cast<int64_t>(i));
    ctx.evalX2k   = MakeEvalMonomial(like, 2 * static_cast<int64_t>(ctx.idx.k));
    ctx.monoM2k   = MakeEvalMonomial(like, -2 * static_cast<int64_t>(ctx.idx.k));
    ctx.monoM1    = MakeEvalMonomial(like, -1);
    ctx.haveMonos = true;

    ctx.logd = 0;
    while ((1u << ctx.logd) < d)
        ++ctx.logd;
    const uint32_t nTow = ctx.evalX2k.GetNumOfElements();
    const uint32_t NN   = ctx.evalX2k.GetRingDimension();
    ctx.wlenFwd.assign(nTow, std::vector<NativeInteger>(static_cast<size_t>(NN) * ctx.logd));
    ctx.wlenInv.assign(nTow, std::vector<NativeInteger>(static_cast<size_t>(NN) * ctx.logd));
    for (uint32_t tow = 0; tow < nTow; ++tow) {
        const NativePoly& wf = ctx.evalX2k.GetElementAtIndex(tow);
        const NativePoly& wi = ctx.monoM2k.GetElementAtIndex(tow);
        const NativeInteger q = wf.GetModulus();
        for (uint32_t sto = 0; sto < NN; ++sto) {
            uint32_t stage = 0;
            for (uint32_t len = 2; len <= d; len <<= 1, ++stage) {
                ctx.wlenFwd[tow][static_cast<size_t>(sto) * ctx.logd + stage] = ModExp(wf[sto], d / len, q);
                ctx.wlenInv[tow][static_cast<size_t>(sto) * ctx.logd + stage] = ModExp(wi[sto], d / len, q);
            }
        }
    }
    ctx.haveWlen = true;
    return ctx;
}

CMTContext MakeCMTContext(uint32_t N, uint32_t d, uint32_t multDepth, uint32_t scaleModSize) {
    if (!CKLConditionSatisfied(N, d))
        OPENFHE_THROW("CKL condition d^2 | 2N is not satisfied");
    CCParams<CryptoContextCKKSRNS> parms;
    parms.SetSecurityLevel(HEStd_NotSet);
    parms.SetRingDim(N);
    parms.SetMultiplicativeDepth(multDepth);
    parms.SetScalingModSize(scaleModSize);
    parms.SetBatchSize(N / 2);
    parms.SetScalingTechnique(FIXEDMANUAL);
    parms.SetKeySwitchTechnique(HYBRID);

    CMTContext ctx;
    ctx.cc = GenCryptoContext(parms);
    return FinishCMTContext(std::move(ctx), N, d, scaleModSize);
}

CMTContext MakeCMTContextForSet(uint32_t d, const CKKSParamSet& ps) {
    CCParams<CryptoContextCKKSRNS> parms;
    parms.SetScalingTechnique(FIXEDMANUAL);
    parms.SetKeySwitchTechnique(HYBRID);
    parms.SetScalingModSize(ps.scaleBits);
    parms.SetFirstModSize(ps.firstModBits);
    parms.SetMultiplicativeDepth(ps.depth);
    if (ps.dnum)
        parms.SetNumLargeDigits(ps.dnum);
    if (ps.sparseHwt)
        parms.SetSecretKeyDist(SPARSE_TERNARY);
    else
        parms.SetSecretKeyDist(UNIFORM_TERNARY);

    uint32_t N = 0;
    if (ps.std128) {
        parms.SetSecurityLevel(HEStd_128_classic);
        if (ps.logN)
            parms.SetRingDim(1u << ps.logN);
    }
    else {
        parms.SetSecurityLevel(HEStd_NotSet);
        if (!ps.logN)
            OPENFHE_THROW("CKKSParamSet: logN required when not STD128");
        N = 1u << ps.logN;
        parms.SetRingDim(N);
        parms.SetBatchSize(N / 2);
    }

    CMTContext ctx;
    ctx.cc = GenCryptoContext(parms);
    if (!N)
        N = ctx.cc->GetRingDimension();
    if (!CKLConditionSatisfied(N, d))
        OPENFHE_THROW("CKL condition d^2 | 2N fails for N=" + std::to_string(N) + " d=" + std::to_string(d));
    return FinishCMTContext(std::move(ctx), N, d, ps.scaleBits);
}

uint64_t EvalFormCiphertextBytes(const Ciphertext<DCRTPoly>& ct) {
    uint64_t n = 0;
    for (const auto& el : ct->GetElements())
        n += static_cast<uint64_t>(el.GetNumOfElements()) * el.GetRingDimension() * sizeof(uint64_t);
    return n;
}

std::vector<Ciphertext<DCRTPoly>> EncryptMatrixColumns(const CMTContext& ctx,
                                                       const std::vector<std::vector<int64_t>>& M) {
    const uint32_t d = ctx.idx.d;
    auto proto =
        ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(
                                    std::vector<double>(ctx.cc->GetEncodingParams()->GetBatchSize(), 0.0)));
    std::vector<DCRTPoly> msgs(d);
    for (uint32_t j = 0; j < d; ++j) {
        msgs[j] = proto->GetElements()[0];
        EncodeColumnPoly(msgs[j], M, j, ctx.encodeScale);
    }
    return EncryptPolys(ctx, msgs);
}

std::vector<Ciphertext<DCRTPoly>> EncryptMatrixRows(const CMTContext& ctx,
                                                    const std::vector<std::vector<int64_t>>& M) {
    const uint32_t d = static_cast<uint32_t>(M.size());
    auto proto =
        ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(
                                    std::vector<double>(ctx.cc->GetEncodingParams()->GetBatchSize(), 0.0)));
    std::vector<DCRTPoly> msgs(d);
    for (uint32_t i = 0; i < d; ++i) {
        msgs[i] = proto->GetElements()[0];
        EncodeRowPoly(msgs[i], M, i, ctx.encodeScale);
    }
    return EncryptPolys(ctx, msgs);
}

std::vector<std::vector<double>> DecryptMatrixColumns(const CMTContext& ctx,
                                                      const std::vector<Ciphertext<DCRTPoly>>& cts) {
    const uint32_t d = ctx.idx.d;
    std::vector<std::vector<double>> M(d, std::vector<double>(cts.size(), 0.0));
    for (uint32_t j = 0; j < cts.size(); ++j) {
        DCRTPoly dec = ctx.cc->GetScheme()->DecryptCore(cts[j], ctx.sk);
        dec.SetFormat(Format::COEFFICIENT);
        Poly large     = dec.CRTInterpolate();
        const auto q   = large.GetModulus();
        const auto half = q >> 1;
        for (uint32_t i = 0; i < d; ++i) {
            auto c = large[i];
            double val;
            if (c > half)
                val = -(q - c).ConvertToDouble();
            else
                val = c.ConvertToDouble();
            M[i][j] = val / static_cast<double>(ctx.encodeScale);
        }
    }
    return M;
}

std::vector<std::vector<int64_t>> CleartextCMT(const std::vector<std::vector<int64_t>>& M) {
    const uint32_t r = static_cast<uint32_t>(M.size());
    const uint32_t c = static_cast<uint32_t>(M[0].size());
    std::vector<std::vector<int64_t>> T(c, std::vector<int64_t>(r, 0));
    for (uint32_t i = 0; i < r; ++i)
        for (uint32_t j = 0; j < c; ++j)
            T[j][i] = M[i][j];
    return T;
}

std::vector<Ciphertext<DCRTPoly>> RunCMT(const CMTContext& ctx, const std::vector<Ciphertext<DCRTPoly>>& cols,
                                         CMTVariant variant, CMTTimings* timings) {
    const auto tAll = Clock::now();
    CMTTimings local;
    CMTTimings* tm = timings ? timings : &local;
    *tm            = CMTTimings{};

    const auto& idx  = ctx.idx;
    const uint32_t d = idx.d;
    const uint32_t N = idx.N;
    const uint32_t k = idx.k;
    if (cols.size() != d)
        OPENFHE_THROW("CMT expects d column ciphertexts");

    std::vector<Ciphertext<DCRTPoly>> aux(d);
#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < static_cast<int64_t>(d); ++i)
        aux[static_cast<uint32_t>(i)] = cols[static_cast<uint32_t>(i)]->Clone();

    const bool spectralPre =
        (variant == CMTVariant::SPECTRAL_TWEAK || variant == CMTVariant::FCMT || variant == CMTVariant::AFCMT);

    // Stay in EVALUATION for the fused path: X^i is a pointwise multiply by a
    // one-hot monomial (one NTT of a sparse poly), not an INTT/rotate/NTT round trip.
    auto t0 = Clock::now();
    if (spectralPre) {
        SetAllFormat(aux, Format::EVALUATION);
#pragma omp parallel for schedule(static)
        for (int64_t i = 1; i < static_cast<int64_t>(d); ++i) {
            const DCRTPoly& mono =
                ctx.haveMonos ? ctx.evalXi[static_cast<uint32_t>(i)] :
                                MakeEvalMonomial(aux[0]->GetElements()[0], i);
            for (auto& el : aux[static_cast<uint32_t>(i)]->GetElements())
                el *= mono;
        }
    }
    else {
#pragma omp parallel for schedule(static)
        for (int64_t i = 1; i < static_cast<int64_t>(d); ++i)
            MultByMonomialInPlace(aux[static_cast<uint32_t>(i)], i, false);
    }
    tm->t_monomial_ms += MsSince(t0);

    t0 = Clock::now();
    ApplyTWEAK(aux, N, k, +1, spectralPre, spectralPre ? d : 0,
               (spectralPre && ctx.haveMonos) ? &ctx.evalX2k : nullptr,
               idx.dft_rev.empty() ? nullptr : idx.dft_rev.data(),
               (spectralPre && ctx.haveWlen) ? &ctx.wlenFwd : nullptr, ctx.logd);
    tm->t_forward_tweak_ms += MsSince(t0);
    tm->t_tweak_fwd_ms += tm->t_forward_tweak_ms;
    if (spectralPre)
        tm->t_spectral_tweak_ms += tm->t_forward_tweak_ms;

    if (!spectralPre) {
#pragma omp parallel for schedule(static)
        for (int64_t t = 0; t < static_cast<int64_t>(d); ++t)
            ScaleByDInvInPlace(aux[static_cast<uint32_t>(t)], d);
    }

    std::vector<Ciphertext<DCRTPoly>> out;
    if (variant == CMTVariant::FCMT || variant == CMTVariant::AFCMT) {
        auto u = BatchKeySwitchNoAuto(ctx.cc, aux, idx.h_for, tm);
        const AFCMTMode mode =
            (variant == CMTVariant::FCMT) ? AFCMTMode::AFCMT_GENERIC : ctx.afcmtMode;
        out = SpectralConsume(ctx, u, tm, mode);
    }
    else {
        std::vector<Ciphertext<DCRTPoly>> scrambled =
            (variant == CMTVariant::CKL_REFERENCE) ? BatchCKLAutomorphism(ctx.cc, aux, idx, tm) :
                                                     BatchHoistedAutomorphism(ctx.cc, aux, idx, tm);
        t0                 = Clock::now();
        const bool specInv = (variant == CMTVariant::SPECTRAL_TWEAK);
        ApplyTWEAK(scrambled, N, k, -1, specInv, 0,
                   (specInv && ctx.haveMonos) ? &ctx.monoM2k : nullptr,
                   idx.dft_rev.empty() ? nullptr : idx.dft_rev.data(),
                   (specInv && ctx.haveWlen) ? &ctx.wlenInv : nullptr, ctx.logd);
        tm->t_tweak_inv_ms += MsSince(t0);
        if (specInv)
            tm->t_spectral_tweak_ms += tm->t_tweak_inv_ms;
        t0 = Clock::now();
#pragma omp parallel for schedule(static)
        for (int64_t i = 0; i < static_cast<int64_t>(d); ++i)
            MultByMonomialInPlace(scrambled[static_cast<uint32_t>(i)], -i, specInv);
        tm->t_monomial_ms += MsSince(t0);
        SetAllFormat(scrambled, Format::EVALUATION);
        out = std::move(scrambled);
    }
    tm->t_cmt_ms = MsSince(tAll);
    return out;
}

std::vector<Ciphertext<DCRTPoly>> RunCCMM(const CMTContext& ctx, const std::vector<Ciphertext<DCRTPoly>>& leftCols,
                                          const std::vector<Ciphertext<DCRTPoly>>& rightCols, CMTVariant variant,
                                          CCMMTimings* timings, bool skipFirstCMT, bool fastKernels) {
    (void)fastKernels;  // always the same PPMM/layout/relin for every variant
    const auto tAll = Clock::now();
    CCMMTimings local;
    CCMMTimings* tm = timings ? timings : &local;
    *tm             = CCMMTimings{};
    const uint32_t d = ctx.idx.d;
    const uint32_t k = ctx.idx.k;

    std::vector<Ciphertext<DCRTPoly>> rightR;
    if (skipFirstCMT) {
        // Row-oriented right operand: CMT_1 is absent by construction.
        tm->t_cmt_ms[0] = 0;
        rightR          = rightCols;
    }
    else {
        auto t0         = Clock::now();
        rightR          = RunCMT(ctx, rightCols, variant, &tm->cmt[0]);
        tm->t_cmt_ms[0] = MsSince(t0);
        tm->n_cmt += 1;
    }

    std::vector<DCRTPoly> LB(d), LA(d), RB(d), RA(d);
    auto t0 = Clock::now();
    for (uint32_t j = 0; j < d; ++j) {
        LB[j] = leftCols[j]->GetElements()[0];
        LA[j] = leftCols[j]->GetElements()[1];
        RB[j] = rightR[j]->GetElements()[0];
        RA[j] = rightR[j]->GetElements()[1];
    }
    RB = FastModuleTranspose(RB, d, k);
    RA = FastModuleTranspose(RA, d, k);
    tm->t_layout_ms += MsSince(t0);

    // Same fused 4-product PPMM for every CMT backend (fair).
    std::vector<DCRTPoly> C00, C01, C10, C11;
    t0 = Clock::now();
    FusedPPMM4(std::move(LB), std::move(LA), std::move(RB), std::move(RA), C00, C01, C10, C11, d, k);
    tm->t_ppmm_ms[0] = MsSince(t0);
    tm->t_ppmm_ms[1] = 0;
    tm->t_ppmm_ms[2] = 0;
    tm->t_ppmm_ms[3] = 0;
    tm->n_ppmm += 4;

    t0  = Clock::now();
    C00 = FastModuleTranspose(C00, d, k);
    C01 = FastModuleTranspose(C01, d, k);
    C10 = FastModuleTranspose(C10, d, k);
    C11 = FastModuleTranspose(C11, d, k);
    std::vector<Ciphertext<DCRTPoly>> pair01(d), pair23(d);
    for (uint32_t j = 0; j < d; ++j) {
        pair01[j] = PackCipher(leftCols[0], C00[j], C01[j]);
        pair23[j] = PackCipher(leftCols[0], C10[j], C11[j]);
    }
    tm->t_layout_ms += MsSince(t0);

    t0       = Clock::now();
    auto D01 = RunCMT(ctx, pair01, variant, &tm->cmt[1]);
    tm->t_cmt_ms[1] = MsSince(t0);
    tm->n_cmt += 1;
    t0       = Clock::now();
    auto D23 = RunCMT(ctx, pair23, variant, &tm->cmt[2]);
    tm->t_cmt_ms[2] = MsSince(t0);
    tm->n_cmt += 1;

    t0 = Clock::now();
    std::vector<Ciphertext<DCRTPoly>> E(d);
#pragma omp parallel for schedule(static)
    for (int64_t j64 = 0; j64 < static_cast<int64_t>(d); ++j64) {
        const uint32_t j = static_cast<uint32_t>(j64);
        auto ct          = D23[j]->Clone();
        DCRTPoly z(ct->GetElements()[0].GetParams(), Format::EVALUATION, true);
        ct->SetElements({z, z, D23[j]->GetElements()[1]});
        ctx.cc->RelinearizeInPlace(ct);
        E[j] = ct;
    }
    tm->t_relin_ms = MsSince(t0);
    tm->n_relin    = d;

    t0 = Clock::now();
    std::vector<Ciphertext<DCRTPoly>> res(d);
#pragma omp parallel for schedule(static)
    for (int64_t j64 = 0; j64 < static_cast<int64_t>(d); ++j64) {
        const uint32_t j = static_cast<uint32_t>(j64);
        auto r           = D01[j]->Clone();
        r->GetElements()[0] += E[j]->GetElements()[0];
        r->GetElements()[1] += E[j]->GetElements()[1];
        r->GetElements()[1] += D23[j]->GetElements()[0];
        res[j] = r;
    }
    tm->t_accumulate_ms = MsSince(t0);

    tm->t_cmt_sum_ms    = tm->t_cmt_ms[0] + tm->t_cmt_ms[1] + tm->t_cmt_ms[2];
    tm->t_ppmm_sum_ms   = tm->t_ppmm_ms[0] + tm->t_ppmm_ms[1] + tm->t_ppmm_ms[2] + tm->t_ppmm_ms[3];
    tm->t_cmt_total_ms  = tm->t_cmt_sum_ms;
    tm->t_ppmm_ms_legacy = tm->t_ppmm_sum_ms;
    tm->t_ccmm_ms       = MsSince(tAll);
    return res;
}

std::vector<Ciphertext<DCRTPoly>> RunCCMMRectangular(const CMTContext& ctx,
                                                     const std::vector<std::vector<int64_t>>& A,
                                                     const std::vector<std::vector<int64_t>>& B, CMTVariant variant,
                                                     CCMMTimings* timings) {
    const uint32_t d = ctx.idx.d;
    if (A.size() != d || B.size() != 2 * d || A[0].size() != 2 * d || B[0].size() != d)
        OPENFHE_THROW("rectangular CCMM expects d x 2d times 2d x d");

    std::vector<std::vector<int64_t>> A0(d, std::vector<int64_t>(d)), A1(d, std::vector<int64_t>(d));
    std::vector<std::vector<int64_t>> B0(d, std::vector<int64_t>(d)), B1(d, std::vector<int64_t>(d));
    for (uint32_t i = 0; i < d; ++i)
        for (uint32_t j = 0; j < d; ++j) {
            A0[i][j] = A[i][j];
            A1[i][j] = A[i][j + d];
            B0[i][j] = B[j][j];  // placeholder overwritten below
        }
    for (uint32_t i = 0; i < d; ++i)
        for (uint32_t j = 0; j < d; ++j) {
            B0[i][j] = B[i][j];
            B1[i][j] = B[i + d][j];
        }

    CCMMTimings t0s, t1s;
    auto L0 = EncryptMatrixColumns(ctx, A0);
    auto R0 = EncryptMatrixColumns(ctx, B0);
    auto P0 = RunCCMM(ctx, L0, R0, variant, &t0s);
    auto L1 = EncryptMatrixColumns(ctx, A1);
    auto R1 = EncryptMatrixColumns(ctx, B1);
    auto P1 = RunCCMM(ctx, L1, R1, variant, &t1s);

    std::vector<Ciphertext<DCRTPoly>> out(d);
    for (uint32_t j = 0; j < d; ++j) {
        out[j] = P0[j]->Clone();
        out[j]->GetElements()[0] += P1[j]->GetElements()[0];
        out[j]->GetElements()[1] += P1[j]->GetElements()[1];
    }
    if (timings) {
        *timings = t0s;
        for (int i = 0; i < 4; ++i)
            timings->t_ppmm_ms[i] += t1s.t_ppmm_ms[i];
        for (int i = 0; i < 3; ++i)
            timings->t_cmt_ms[i] += t1s.t_cmt_ms[i];
        timings->t_cmt_sum_ms += t1s.t_cmt_sum_ms;
        timings->t_ppmm_sum_ms += t1s.t_ppmm_sum_ms;
        timings->t_cmt_total_ms += t1s.t_cmt_total_ms;
        timings->t_relin_ms += t1s.t_relin_ms;
        timings->t_layout_ms += t1s.t_layout_ms;
        timings->t_accumulate_ms += t1s.t_accumulate_ms;
        timings->t_ccmm_ms += t1s.t_ccmm_ms;
        timings->n_cmt += t1s.n_cmt;
        timings->n_ppmm += t1s.n_ppmm;
        timings->n_relin += t1s.n_relin;
    }
    return out;
}

// ---------------------------------------------------------------------------
// Mathematical unit tests A–F
// ---------------------------------------------------------------------------

namespace {

bool EqNative(const NativeInteger& a, const NativeInteger& b) {
    return a == b;
}

MathTestReport Pass(const std::string& name, const std::string& d) {
    return {true, name, d};
}
MathTestReport Fail(const std::string& name, const std::string& d) {
    return {false, name, d};
}

// Small single-modulus ring helpers for Test A.
void PolyMulNegacyclic(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b, std::vector<uint64_t>& c,
                       uint64_t q) {
    const uint32_t n = static_cast<uint32_t>(a.size());
    c.assign(n, 0);
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < n; ++j) {
            const uint64_t prod = (a[i] * b[j]) % q;
            const uint32_t s    = i + j;
            if (s < n)
                c[s] = (c[s] + prod) % q;
            else
                c[s - n] = (c[s - n] + q - prod) % q;
        }
}

MathTestReport TestA_ModuleTensor() {
    // N=16, d=4, k=4, q a small prime. Direct ring product vs the Y-twisted tensor.
    const uint32_t N = 16, d = 4, k = 4;
    const uint64_t q = 7681;
    // ζ = Y^{2k/d} = Y^{2}. Y = X^4, ζ = X^8. Order d=4: (X^8)^4 = X^{32} = (X^{16})^2 = 1.
    std::vector<uint64_t> A(N), W(N), Kdir(N), Kten(N, 0);
    for (uint32_t i = 0; i < N; ++i) {
        A[i] = (3ull * i + 5) % q;
        W[i] = (7ull * i + 2) % q;
    }
    const uint32_t t = 1;
    // σ_t(A): coord i of A multiplied by ζ^{it}. ζ = X^8, ζ^{it} = X^{8 i t}.
    std::vector<uint64_t> sigA(N, 0);
    for (uint32_t i = 0; i < d; ++i) {
        for (uint32_t s = 0; s < k; ++s) {
            const uint32_t src = s * d + i;
            uint32_t dst       = src + (8 * i * t);
            int sign           = 1;
            dst %= (2 * N);
            if (dst >= N) {
                dst -= N;
                sign = -1;
            }
            uint64_t v = A[src];
            if (sign < 0)
                v = (q - v) % q;
            sigA[dst] = (sigA[dst] + v) % q;
        }
    }
    PolyMulNegacyclic(sigA, W, Kdir, q);

    // Tensor: for each output coord r, sum_i ζ^{it} A_i W_{(r-i) mod d} Y^{i>r}
    // A_i(Y) is coeffs A[s*d+i]; W_p(Y) is W[s*d+p].
    auto getCoord = [&](const std::vector<uint64_t>& p, uint32_t coord, uint32_t s) { return p[s * d + coord]; };
    for (uint32_t r = 0; r < d; ++r) {
        for (uint32_t i = 0; i < d; ++i) {
            const uint32_t wp    = (r + d - i) % d;
            const uint32_t wrapY = (i > r) ? 1 : 0;
            // A_i * W_wp as R_k mul (negacyclic length k), then * Y^{wrap} = X^{d wrap}
            std::vector<uint64_t> Ai(k), Ww(k), prod(k, 0);
            for (uint32_t s = 0; s < k; ++s) {
                Ai[s] = getCoord(A, i, s);
                Ww[s] = getCoord(W, wp, s);
            }
            for (uint32_t s1 = 0; s1 < k; ++s1)
                for (uint32_t s2 = 0; s2 < k; ++s2) {
                    uint64_t pr = (Ai[s1] * Ww[s2]) % q;
                    uint32_t ss = s1 + s2;
                    if (ss < k)
                        prod[ss] = (prod[ss] + pr) % q;
                    else
                        prod[ss - k] = (prod[ss - k] + q - pr) % q;
                }
            // multiply by ζ^{it} = X^{8 i t} and place on coordinate r, times Y^{wrap}=X^{4 wrap}
            const uint32_t mon = (8 * i * t + 4 * wrapY) % (2 * N);
            for (uint32_t s = 0; s < k; ++s) {
                uint32_t dst = s * d + r + mon;
                int sign     = 1;
                dst %= (2 * N);
                if (dst >= N) {
                    dst -= N;
                    sign = -1;
                }
                uint64_t v = prod[s];
                if (sign < 0)
                    v = (q - v) % q;
                Kten[dst] = (Kten[dst] + v) % q;
            }
        }
    }
    for (uint32_t i = 0; i < N; ++i)
        if (Kdir[i] != Kten[i])
            return Fail("A_module_tensor", "mismatch at coeff " + std::to_string(i));
    return Pass("A_module_tensor", "N=16 d=4 exact match");
}

}  // namespace

std::vector<MathTestReport> RunMathUnitTests(uint32_t N, uint32_t d, CMTContext* ctxOrNull) {
    std::vector<MathTestReport> reports;
    reports.push_back(TestA_ModuleTensor());

    if (!CKLConditionSatisfied(N, d)) {
        reports.push_back(Fail("param", "invalid (N,d)"));
        return reports;
    }
    const uint32_t k = N / d;

    // Test E: π_direct == π_affine, and general lane == simplified under CKL.
    {
        uint32_t misPi = 0, misLane = 0;
        for (uint32_t t = 0; t < d; ++t)
            for (uint32_t nu = 0; nu < N; ++nu) {
                if (PiNatural(nu, t, N, k) != PiAffine(nu, t, N, k, d))
                    ++misPi;
                const uint32_t a = nu % k, b = nu / k;
                if (LaneAffineGeneral(a, b, t, k, d) != LaneAffineSimplified(a, b, t, d))
                    ++misLane;
            }
        const auto cls = CountTransformClasses(N, d);
        reports.push_back(misPi == 0 && misLane == 0 ?
                              Pass("E_affine_formula", "PiNatural=PiAffine and general=simplified") :
                              Fail("E_affine_formula", "pi=" + std::to_string(misPi) + " lane=" + std::to_string(misLane)));
        reports.push_back(cls.gamma_ok ? Pass("H_gamma_classes",
                                              "unique=" + std::to_string(cls.unique_gamma) +
                                                  " expected=" + std::to_string(cls.expected_gamma) +
                                                  " count=" + std::to_string(cls.min_gamma_count) +
                                                  " expected_count=" + std::to_string(cls.expected_gamma_count)) :
                                         Fail("H_gamma_classes", "unique/count mismatch"));
        reports.push_back(cls.step_ok ? Pass("I_step_classes",
                                             "unique=" + std::to_string(cls.unique_step) +
                                                 " expected=" + std::to_string(cls.expected_step) +
                                                 " count=" + std::to_string(cls.min_step_count) +
                                                 " expected_count=" + std::to_string(cls.expected_step_count)) :
                                        Fail("I_step_classes", "unique/count mismatch"));
    }

    if (!ctxOrNull)
        return reports;
    CMTContext& ctx = *ctxOrNull;
    const auto& idx = ctx.idx;

    // Test B: spectral monomial vs OpenFHE eval multiply.
    {
        DCRTPoly like = ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(std::vector<double>(N / 2, 0.0)))
                            ->GetElements()[0];
        DCRTPoly f(like.GetParams(), Format::COEFFICIENT, true);
        for (auto& tw : f.GetAllElements())
            for (uint32_t i = 0; i < std::min<uint32_t>(17, tw.GetLength()); ++i)
                tw[i] = (i * 3 + 1);
        f.SetFormat(Format::EVALUATION);
        const int64_t u   = 5;
        DCRTPoly mono     = MakeEvalMonomial(f, u);
        DCRTPoly prod     = f;
        prod *= mono;
        DCRTPoly ref = f;
        ref.SetFormat(Format::COEFFICIENT);
        NegacyclicMulXInPlace(ref, u);
        ref.SetFormat(Format::EVALUATION);
        bool ok = true;
        for (uint32_t tow = 0; tow < prod.GetNumOfElements() && ok; ++tow) {
            const auto& a = prod.GetElementAtIndex(tow);
            const auto& b = ref.GetElementAtIndex(tow);
            for (uint32_t i = 0; i < a.GetLength(); ++i)
                if (!EqNative(a[i], b[i])) {
                    ok = false;
                    break;
                }
        }
        reports.push_back(ok ? Pass("B_spectral_monomial", "X^5 f") : Fail("B_spectral_monomial", "eval mismatch"));
    }

    // Test C: coeff TWEAK(-1) == spectral DFT on d dummy streams of one poly family.
    {
        std::vector<Ciphertext<DCRTPoly>> cts(d), ctsC, ctsS;
        auto zeros = ctx.cc->MakeCKKSPackedPlaintext(std::vector<double>(N / 2, 0.0));
        for (uint32_t t = 0; t < d; ++t) {
            cts[t] = ctx.cc->Encrypt(ctx.pk, zeros);
            DCRTPoly msg(cts[t]->GetElements()[0].GetParams(), Format::COEFFICIENT, true);
            for (auto& tw : msg.GetAllElements())
                for (uint32_t i = 0; i < d && i < tw.GetLength(); ++i)
                    tw[i] = NativeInteger((t + 1) * (i + 3));
            msg.SetFormat(Format::EVALUATION);
            cts[t]->GetElements()[0] += msg;
        }
        ctsC = cts;
        ctsS = cts;
        ApplyTWEAK(ctsC, N, k, -1, false);
        ApplyTWEAK(ctsS, N, k, -1, true);
        bool ok = true;
        for (uint32_t t = 0; t < d && ok; ++t) {
            auto a = ctsC[t]->GetElements()[0];
            auto b = ctsS[t]->GetElements()[0];
            a.SetFormat(Format::COEFFICIENT);
            b.SetFormat(Format::COEFFICIENT);
            for (uint32_t tow = 0; tow < a.GetNumOfElements() && ok; ++tow)
                for (uint32_t i = 0; i < a.GetRingDimension(); ++i)
                    if (!EqNative(a.GetElementAtIndex(tow)[i], b.GetElementAtIndex(tow)[i]))
                        ok = false;
        }
        reports.push_back(ok ? Pass("C_inv_tweak_dft", "coeff == spectral") : Fail("C_inv_tweak_dft", "mismatch"));
    }

    // Test D: AutomorphismTransform vs π_t on evaluation form.
    {
        auto ct = ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(std::vector<double>(N / 2, 0.0)));
        DCRTPoly f = ct->GetElements()[0];
        f.SetFormat(Format::COEFFICIENT);
        for (auto& tw : f.GetAllElements())
            for (uint32_t i = 0; i < tw.GetLength(); ++i)
                tw[i] = NativeInteger((i * 9 + 4) % 17 + 1);
        f.SetFormat(Format::EVALUATION);
        bool ok = true;
        std::ostringstream det;
        for (uint32_t t = 0; t < d && ok; ++t) {
            const uint32_t h = idx.h[t];
            std::vector<uint32_t> vec(N);
            PrecomputeAutoMap(N, h, &vec);
            DCRTPoly autoed = f.AutomorphismTransform(h, vec);
            const auto& src = f.GetElementAtIndex(0);
            const auto& dst = autoed.GetElementAtIndex(0);
            for (uint32_t nu = 0; nu < N; ++nu) {
                const uint64_t got = dst[idx.bitrev[nu]].ConvertToInt<uint64_t>();
                const uint64_t expv = src[idx.pi_sto[t][nu]].ConvertToInt<uint64_t>();
                if (got != expv) {
                    ok = false;
                    det << "t=" << t << " nu=" << nu;
                    break;
                }
            }
        }
        reports.push_back(ok ? Pass("D_auto_permutation", "OpenFHE Auto == f(π_t(ν))") :
                               Fail("D_auto_permutation", det.str()));
    }

    // Test F: γ^d = 1 and γ^{d/2} = -1; Test G: γ_{ν+d/2}=γ_ν
    {
        DCRTPoly like = ctx.cc->Encrypt(ctx.pk, ctx.cc->MakeCKKSPackedPlaintext(std::vector<double>(N / 2, 0.0)))
                            ->GetElements()[0];
        DCRTPoly g = MakeEvalMonomial(like, -2 * static_cast<int64_t>(k));
        const auto& tw = g.GetElementAtIndex(0);
        const NativeInteger q = tw.GetModulus();
        bool okF = true, okG = true;
        for (uint32_t nu = 0; nu < N; ++nu) {
            NativeInteger gn = tw[idx.bitrev[nu]];
            NativeInteger acc(1);
            for (uint32_t p = 0; p < d; ++p)
                acc.ModMulEq(gn, q);
            if (acc != NativeInteger(1))
                okF = false;
            if (d >= 2) {
                NativeInteger half(1);
                for (uint32_t p = 0; p < d / 2; ++p)
                    half.ModMulEq(gn, q);
                if (half != q - 1)
                    okF = false;
                const uint32_t n2 = (nu + d / 2) % N;
                if (tw[idx.bitrev[nu]] != tw[idx.bitrev[n2]])
                    okG = false;
            }
        }
        reports.push_back(okF ? Pass("F_gamma_order", "γ^d=1 and γ^{d/2}=-1") : Fail("F_gamma_order", "order failed"));
        reports.push_back(okG ? Pass("G_gamma_half_period", "γ_{ν+d/2}=γ_ν") : Fail("G_gamma_half_period", "period failed"));
    }

    return reports;
}

}  // namespace lbcrypto

namespace lbcrypto {
std::vector<Ciphertext<DCRTPoly>> RunRectangularCCMM(
    const CMTContext& ctx,
    const std::vector<Ciphertext<DCRTPoly>>& leftCols,
    const std::vector<Ciphertext<DCRTPoly>>& rightCols,
    uint32_t B,
    CMTVariant variant,
    CCMMTimings* tm) {
    
    if (!tm) OPENFHE_THROW("RunRectangularCCMM requires timings object");
    
    const uint32_t d = ctx.idx.d;
    auto t0 = Clock::now();
    
    // Accumulators for the raw polynomial pre-CMT state
    std::vector<DCRTPoly> sumC00(d), sumC01(d), sumC10(d), sumC11(d);
    
    for (uint32_t b = 0; b < B; ++b) {
        auto t_cmt = Clock::now();
        std::vector<Ciphertext<DCRTPoly>> rightR = RunCMT(ctx, rightCols, variant, &tm->cmt[0]);
        tm->t_cmt_ms[0] += MsSince(t_cmt);
        tm->n_cmt += 1;
        std::vector<DCRTPoly> LB(d), LA(d), RB(d), RA(d);
        for (uint32_t j = 0; j < d; ++j) {
            LB[j] = leftCols[j]->GetElements()[0];
            LA[j] = leftCols[j]->GetElements()[1];
            RB[j] = rightR[j]->GetElements()[0];
            RA[j] = rightR[j]->GetElements()[1];
        }
        
        t0 = Clock::now();
        std::vector<DCRTPoly> C00(d), C01(d), C10(d), C11(d);
        FusedPPMM4(LB, LA, RB, RA, C00, C01, C10, C11, d, ctx.idx.k);
        tm->t_ppmm_sum_ms += MsSince(t0);
        tm->n_ppmm += 1;
        
        t0 = Clock::now();
        #pragma omp parallel for schedule(static)
        for (int64_t j = 0; j < static_cast<int64_t>(d); ++j) {
            if (b == 0) {
                sumC00[j] = std::move(C00[j]);
                sumC01[j] = std::move(C01[j]);
                sumC10[j] = std::move(C10[j]);
                sumC11[j] = std::move(C11[j]);
            } else {
                sumC00[j] += C00[j];
                sumC01[j] += C01[j];
                sumC10[j] += C10[j];
                sumC11[j] += C11[j];
            }
        }
        tm->t_accumulate_ms += MsSince(t0);
    }
    
    t0 = Clock::now();
    std::vector<Ciphertext<DCRTPoly>> sum01(d), sum23(d);
    #pragma omp parallel for schedule(static)
    for (int64_t j = 0; j < static_cast<int64_t>(d); ++j) {
        sum01[j] = PackCipher(leftCols[0], sumC00[j], sumC01[j]);
        sum23[j] = PackCipher(leftCols[0], sumC10[j], sumC11[j]);
    }
    auto D01 = RunCMT(ctx, sum01, variant, &tm->cmt[1]);
    tm->t_cmt_ms[1] = MsSince(t0);
    tm->n_cmt += 1;
    
    t0 = Clock::now();
    auto D23 = RunCMT(ctx, sum23, variant, &tm->cmt[2]);
    tm->t_cmt_ms[2] = MsSince(t0);
    tm->n_cmt += 1;
    
    t0 = Clock::now();
    std::vector<Ciphertext<DCRTPoly>> res(d);
    std::vector<Ciphertext<DCRTPoly>> E(d);
#pragma omp parallel for schedule(static)
    for (int64_t j64 = 0; j64 < static_cast<int64_t>(d); ++j64) {
        const uint32_t j = static_cast<uint32_t>(j64);
        auto ct = D23[j]->Clone();
        DCRTPoly z(ct->GetElements()[0].GetParams(), Format::EVALUATION, true);
        ct->SetElements({z, z, D23[j]->GetElements()[1]});
        ctx.cc->RelinearizeInPlace(ct);
        E[j] = ct;
    }
    tm->t_relin_ms += MsSince(t0);
    tm->n_relin += 1;
    
    t0 = Clock::now();
    for (int64_t i = 0; i < static_cast<int64_t>(d); ++i) {
        const uint32_t j = static_cast<uint32_t>(i);
        auto ct = leftCols[j]->Clone();
        ct->SetElements({D01[j]->GetElements()[0] + E[j]->GetElements()[0],
                         D01[j]->GetElements()[1] + E[j]->GetElements()[1]});
        res[j] = ct;
    }
    tm->t_accumulate_ms += MsSince(t0);
    
    tm->t_cmt_sum_ms = tm->t_cmt_ms[0] + tm->t_cmt_ms[1] + tm->t_cmt_ms[2];
    tm->t_ccmm_ms = tm->t_cmt_sum_ms + tm->t_ppmm_sum_ms + tm->t_relin_ms + tm->t_layout_ms + tm->t_accumulate_ms;
    
    return res;
}

}
