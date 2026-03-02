#!/usr/bin/env python3
"""A12: Spectral proof attempt for Proposition 3 — Corr(c_i, c_{i+1}) = -(b-1)/b.

Strategy:
  Part A — Exact enumeration of carry correlations for K-bit multiplications.
            Compute Corr(c_j, c_{j+1}) at every interior position j.
  Part B — Build the empirical 1-step transition matrix P(c' | c) by
            marginalizing over digits. Compute its eigenvalues.
  Part C — Spectral decomposition: express Cov(c_j, c_{j+1}) as a sum over
            eigenvalues λ_k and test whether only λ_1=1/b contributes.
  Part D — Multi-base test (b=2,3,5): does -(b-1)/b hold universally?
  Part E — Analytical bounds on higher-eigenvalue contributions.

The goal is to identify why higher-eigenvalue contributions vanish, providing
a pathway to a complete proof.
"""

import sys
import time
import math
from collections import defaultdict
from fractions import Fraction

import numpy as np

np.random.seed(42)


def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def multiply_carries(p, q, base=2):
    """Compute all carry values for p × q in given base."""
    digits_p = []
    x = p
    while x > 0:
        digits_p.append(x % base)
        x //= base

    digits_q = []
    x = q
    while x > 0:
        digits_q.append(x % base)
        x //= base

    Kp, Kq = len(digits_p), len(digits_q)
    D = Kp + Kq - 1

    carries = [0] * (D + 1)
    for j in range(D):
        conv_j = 0
        for i in range(max(0, j - Kq + 1), min(j, Kp - 1) + 1):
            conv_j += digits_p[i] * digits_q[j - i]
        total = conv_j + carries[j]
        carries[j + 1] = total // base

    return carries[:D + 1], D


def exact_carry_correlations(K, base=2):
    """Enumerate all K-digit factor pairs and compute carry correlations.

    Returns (positions, means, variances, covariances, correlations).
    """
    lo = base ** (K - 1)
    hi = base ** K

    D = 2 * K - 1
    carry_sum = np.zeros(D + 1, dtype=np.float64)
    carry_sq = np.zeros(D + 1, dtype=np.float64)
    carry_prod = np.zeros(D, dtype=np.float64)
    count = 0

    for p in range(lo, hi):
        for q in range(lo, hi):
            carries, d = multiply_carries(p, q, base)
            assert d == D
            for j in range(D + 1):
                c = carries[j]
                carry_sum[j] += c
                carry_sq[j] += c * c
            for j in range(D):
                carry_prod[j] += carries[j] * carries[j + 1]
            count += 1

    means = carry_sum / count
    variances = carry_sq / count - means ** 2
    covariances = carry_prod / count - means[:D] * means[1:D + 1]

    correlations = np.zeros(D)
    for j in range(D):
        v1 = variances[j]
        v2 = variances[j + 1]
        if v1 > 1e-30 and v2 > 1e-30:
            correlations[j] = covariances[j] / math.sqrt(v1 * v2)
        else:
            correlations[j] = float('nan')

    return means, variances, covariances, correlations, count


def build_transition_matrix(K, j_target, base=2):
    """Build empirical transition matrix P(c_{j+1} | c_j) at position j_target."""
    lo = base ** (K - 1)
    hi = base ** K

    joint = defaultdict(lambda: defaultdict(int))
    margin = defaultdict(int)

    for p in range(lo, hi):
        for q in range(lo, hi):
            carries, _ = multiply_carries(p, q, base)
            if j_target < len(carries) - 1:
                c_j = carries[j_target]
                c_j1 = carries[j_target + 1]
                joint[c_j][c_j1] += 1
                margin[c_j] += 1

    states = sorted(set(list(margin.keys()) +
                        [c1 for cj in joint.values() for c1 in cj.keys()]))
    n = len(states)
    s2i = {s: i for i, s in enumerate(states)}

    T = np.zeros((n, n))
    pi_dist = np.zeros(n)

    total = sum(margin.values())
    for c_j in states:
        i = s2i[c_j]
        pi_dist[i] = margin.get(c_j, 0) / total
        row_total = margin.get(c_j, 0)
        if row_total > 0:
            for c_j1, cnt in joint.get(c_j, {}).items():
                jj = s2i[c_j1]
                T[i, jj] = cnt / row_total

    return T, pi_dist, states


def spectral_correlation(T, pi_dist, states):
    """Compute the correlation from the spectral decomposition of T.

    Corr(f, Tf) = [sum_x pi(x)f(x)(Tf)(x) - E[f]^2] / Var(f)
    where f(x) = x (the identity function on states).
    """
    n = len(states)
    f = np.array(states, dtype=float)

    Ef = np.sum(pi_dist * f)
    Ef2 = np.sum(pi_dist * f ** 2)
    Var_f = Ef2 - Ef ** 2

    Tf = T @ f
    EfTf = np.sum(pi_dist * f * Tf)
    Cov_fTf = EfTf - Ef ** 2

    corr = Cov_fTf / Var_f if Var_f > 1e-30 else float('nan')

    evals = np.linalg.eigvals(T)
    evals_sorted = sorted(evals.real, reverse=True)

    return corr, Cov_fTf, Var_f, evals_sorted


def spectral_decomposition_correlation(T, pi_dist, states):
    """Decompose Cov(f, Tf) into eigenvalue contributions.

    For transition matrix T with right eigenvectors v_k and eigenvalues λ_k:
    Cov(f, Tf) = Σ_k λ_k · w_k
    where w_k = <(f - E[f]), P_k (f - E[f])>_π
    """
    n = len(states)
    f = np.array(states, dtype=float)
    Ef = np.sum(pi_dist * f)
    f_centered = f - Ef
    Var_f = np.sum(pi_dist * f_centered ** 2)

    D_pi = np.diag(np.sqrt(pi_dist + 1e-30))
    D_pi_inv = np.diag(1.0 / np.sqrt(pi_dist + 1e-30))
    T_sym = D_pi @ T @ D_pi_inv

    evals, evecs_sym = np.linalg.eigh(0.5 * (T_sym + T_sym.T))
    idx = np.argsort(-evals)
    evals = evals[idx]
    evecs_sym = evecs_sym[:, idx]

    evecs_right = D_pi_inv @ evecs_sym

    contributions = []
    for k in range(n):
        v_k = evecs_right[:, k]
        overlap = np.sum(pi_dist * f_centered * v_k)
        norm = np.sum(pi_dist * v_k ** 2)
        if abs(norm) > 1e-20:
            w_k = overlap ** 2 / norm
        else:
            w_k = 0.0
        contributions.append((evals[k], w_k, evals[k] * w_k))

    return contributions, Var_f


def exact_rational_correlation(K, j_target, base=2):
    """Compute exact rational Corr(c_j, c_{j+1}) using Fraction arithmetic."""
    lo = base ** (K - 1)
    hi = base ** K
    count = (hi - lo) ** 2

    sum_cj = Fraction(0)
    sum_cj1 = Fraction(0)
    sum_cj_sq = Fraction(0)
    sum_cj1_sq = Fraction(0)
    sum_cj_cj1 = Fraction(0)

    for p in range(lo, hi):
        for q in range(lo, hi):
            carries, _ = multiply_carries(p, q, base)
            if j_target < len(carries) - 1:
                cj = Fraction(carries[j_target])
                cj1 = Fraction(carries[j_target + 1])
                sum_cj += cj
                sum_cj1 += cj1
                sum_cj_sq += cj * cj
                sum_cj1_sq += cj1 * cj1
                sum_cj_cj1 += cj * cj1

    N = Fraction(count)
    E_cj = sum_cj / N
    E_cj1 = sum_cj1 / N
    V_cj = sum_cj_sq / N - E_cj ** 2
    V_cj1 = sum_cj1_sq / N - E_cj1 ** 2
    Cov = sum_cj_cj1 / N - E_cj * E_cj1

    return Cov, V_cj, V_cj1, E_cj, E_cj1


def main():
    t0 = time.time()
    pr("=" * 72)
    pr("A12: SPECTRAL PROOF ATTEMPT FOR ANTI-CORRELATION LAW")
    pr("     Target: Corr(c_j, c_{j+1}) = -(b-1)/b")
    pr("=" * 72)

    # ═══════════════════════════════════════════════════════════════
    # PART A: EXACT CARRY CORRELATIONS
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART A: EXACT CARRY CORRELATIONS (base 2)")
    pr(f"{'═' * 72}")

    target_b2 = -1 / 2

    for K in range(3, 12):
        if K > 9:
            continue
        means, variances, covs, corrs, count = exact_carry_correlations(K, base=2)
        D = 2 * K - 1

        interior_start = max(2, K // 2)
        interior_end = min(D - 1, K + K // 2)
        interior_corrs = [corrs[j] for j in range(interior_start, interior_end)
                          if not math.isnan(corrs[j])]

        if interior_corrs:
            mean_corr = np.mean(interior_corrs)
            std_corr = np.std(interior_corrs)
            delta = mean_corr - target_b2
            pr(f"\n  K={K:2d} ({count:8d} pairs, D={D:2d}):")
            pr(f"    Interior positions {interior_start}..{interior_end-1}:")
            pr(f"    Mean Corr = {mean_corr:+.8f}  (target = {target_b2:+.4f})")
            pr(f"    Δ from target = {delta:+.6e}")
            pr(f"    Std across positions = {std_corr:.6e}")

            pr(f"    Position-by-position:")
            for j in range(max(1, K - 3), min(D - 1, K + 4)):
                if not math.isnan(corrs[j]):
                    pr(f"      j={j:2d}: Corr = {corrs[j]:+.8f}  "
                       f"E[c_j] = {means[j]:.4f}  Var(c_j) = {variances[j]:.4f}")

    # ═══════════════════════════════════════════════════════════════
    # PART B: TRANSITION MATRIX AND EIGENVALUES
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART B: EMPIRICAL TRANSITION MATRIX (base 2)")
    pr(f"{'═' * 72}")

    for K in [6, 7, 8]:
        D = 2 * K - 1
        j_mid = K - 1

        T, pi_dist, states = build_transition_matrix(K, j_mid, base=2)
        corr, cov, var_f, evals = spectral_correlation(T, pi_dist, states)

        pr(f"\n  K={K}, j={j_mid} (mid-position):")
        pr(f"    States: {states[:10]}{'...' if len(states) > 10 else ''}")
        pr(f"    Corr(c_j, c_{{j+1}}) from T: {corr:+.8f}  (target: -0.5)")
        pr(f"    Eigenvalues (top 6): {['%.6f' % e for e in evals[:6]]}")
        pr(f"    Expected: 1.0, 0.5, 0.25, 0.125, ...")

    # ═══════════════════════════════════════════════════════════════
    # PART C: SPECTRAL DECOMPOSITION OF CORRELATION
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART C: SPECTRAL DECOMPOSITION — eigenvalue contributions")
    pr(f"{'═' * 72}")
    pr("  Cov(f, Tf) = Σ_k λ_k · w_k")
    pr("  If only λ₁=1/b contributes → Corr = (1/b)·w₁/Var = -(b-1)/b")

    for K in [7, 8, 9]:
        D = 2 * K - 1
        for j_target in [K - 1, K]:
            if K == 9 and j_target > K - 1:
                continue
            T, pi_dist, states = build_transition_matrix(K, j_target, base=2)
            contribs, var_f = spectral_decomposition_correlation(
                T, pi_dist, states)

            total_cov = sum(c[2] for c in contribs)
            corr_total = total_cov / var_f if var_f > 1e-30 else float('nan')

            pr(f"\n  K={K}, j={j_target}:")
            pr(f"    Var(c) = {var_f:.6f}")
            pr(f"    Total Corr = {corr_total:+.8f}")
            pr(f"    {'λ_k':>10s} {'w_k':>12s} {'λ_k·w_k':>12s} "
               f"{'contrib/Var':>12s} {'cumul':>12s}")
            cumul = 0.0
            for lam, w, lw in contribs:
                if abs(lw) > 1e-12:
                    cumul += lw / var_f
                    pr(f"    {lam:10.6f} {w:12.6f} {lw:12.8f} "
                       f"{lw/var_f:+12.8f} {cumul:+12.8f}")

    # ═══════════════════════════════════════════════════════════════
    # PART D: MULTI-BASE TEST
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART D: MULTI-BASE TEST — Corr = -(b-1)/b?")
    pr(f"{'═' * 72}")

    for base in [2, 3, 5]:
        target = -(base - 1) / base
        max_K = {2: 9, 3: 6, 5: 4}[base]

        pr(f"\n  Base {base}: target Corr = {target:+.6f}")
        for K in range(3, max_K + 1):
            means, variances, covs, corrs, count = exact_carry_correlations(K, base)
            D = 2 * K - 1
            j_mid = K - 1
            c_mid = corrs[j_mid] if j_mid < len(corrs) and not math.isnan(corrs[j_mid]) else float('nan')
            delta = c_mid - target if not math.isnan(c_mid) else float('nan')
            pr(f"    K={K:2d} ({count:6d} pairs): Corr(j={j_mid}) = {c_mid:+.8f}  "
               f"Δ = {delta:+.4e}")

    # ═══════════════════════════════════════════════════════════════
    # PART E: EXACT RATIONAL ANALYSIS
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART E: EXACT RATIONAL ARITHMETIC — small K")
    pr(f"{'═' * 72}")
    pr("  Testing if Cov/Var has exact form -(b-1)/b = -1/2")

    for K in range(3, 8):
        D = 2 * K - 1
        j_mid = K - 1
        Cov, V_cj, V_cj1, E_cj, E_cj1 = exact_rational_correlation(K, j_mid, base=2)

        if V_cj > 0 and V_cj1 > 0:
            ratio = Cov * Cov / (V_cj * V_cj1)
            corr_approx = float(Cov) / math.sqrt(float(V_cj) * float(V_cj1))
            target_cov = -V_cj / 2

            pr(f"\n  K={K}, j={j_mid}:")
            pr(f"    E[c_j] = {E_cj} = {float(E_cj):.6f}")
            pr(f"    Var(c_j) = {V_cj} = {float(V_cj):.6f}")
            pr(f"    Cov(c_j, c_{{j+1}}) = {Cov} = {float(Cov):.8f}")
            pr(f"    Corr = {corr_approx:+.8f}")
            pr(f"    Cov / Var(c_j) = {float(Cov / V_cj):+.8f}  (target = -0.5)")
            pr(f"    Cov + Var/2 = {Cov + V_cj/2} = {float(Cov + V_cj/2):.8e}")
            if Cov + V_cj / 2 != 0:
                pr(f"    → Higher-eigenvalue contribution is NON-ZERO: "
                   f"{float(Cov + V_cj/2):.6e}")
            else:
                pr(f"    → Cov = -Var/2 EXACTLY! Only λ₁=1/2 contributes.")

    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("SUMMARY")
    pr(f"{'═' * 72}")
    pr(f"  Runtime: {time.time() - t0:.1f}s")
    pr(f"\n  Key question: does Cov(c_j, c_{{j+1}}) = -(1/2)·Var(c_j) EXACTLY")
    pr(f"  at finite K, or only asymptotically as K → ∞?")
    pr(f"  If exact at finite K: look for algebraic identity in the spectral")
    pr(f"  decomposition. If asymptotic: bound the rate of convergence.")
    pr("=" * 72)


if __name__ == '__main__':
    main()
