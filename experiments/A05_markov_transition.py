#!/usr/bin/env python3
"""
A05: Carry Markov Chain — Transition Matrix

The multiplication p × q in base l produces carries that propagate
left-to-right. At position k with n contributing digit pairs:
  conv_k = Σ_{i+j=k} g_i · h_j     (convolution)
  carry_{k+1} = floor((conv_k + carry_k) / l)
  digit_k = (conv_k + carry_k) mod l

This is a Markov chain on carry states. We build the exact transition
matrix T and compute its spectrum.

Key questions:
  A) What are the eigenvalues of T? (Diaconis-Fulman for addition: 1, 1/b, 1/b², ...)
  B) How does the spectrum depend on n (number of digit pairs) and base l?
  C) Does the stationary distribution match the empirical carry distribution?
  D) What is the mixing time?
"""

import sys, os, time, math
import numpy as np
from scipy.special import comb as binom_coeff

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from carry_utils import primes_up_to

np.random.seed(42)


def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def product_distribution(b):
    """Distribution of g*h for g, h uniform on {0,...,b-1}.
    Returns array where P[k] = Prob(g*h = k)."""
    max_val = (b - 1) ** 2
    dist = np.zeros(max_val + 1)
    for g in range(b):
        for h in range(b):
            dist[g * h] += 1
    dist /= b * b
    return dist


def nfold_convolution(single_dist, n):
    """N-fold convolution of a discrete distribution using FFT."""
    if n == 0:
        return np.array([1.0])
    if n == 1:
        return single_dist.copy()
    result = single_dist.copy()
    for _ in range(n - 1):
        result = np.convolve(result, single_dist)
    return result


def build_transition_matrix(b, conv_dist):
    """Build carry transition matrix T(c'|c) for base b.
    T[c', c] = Prob(carry goes from c to c').
    conv_dist[v] = Prob(convolution = v)."""
    v_max = len(conv_dist) - 1
    c_max = (v_max + v_max) // b + 1
    c_max = min(c_max, v_max)

    T = np.zeros((c_max + 1, c_max + 1))
    for c in range(c_max + 1):
        for v in range(len(conv_dist)):
            if conv_dist[v] < 1e-30:
                continue
            c_next = (v + c) // b
            if c_next <= c_max:
                T[c_next, c] += conv_dist[v]

    col_sums = T.sum(axis=0)
    active = col_sums > 0.99
    n_active = np.sum(active)
    T_trim = T[np.ix_(active, active)]
    return T_trim, n_active


def main():
    t0 = time.time()
    pr("=" * 72)
    pr("A05: CARRY MARKOV CHAIN — TRANSITION MATRIX")
    pr("=" * 72)

    # ═══════════════════════════════════════════════════════════════
    # PART A: EIGENVALUE SPECTRUM OF T
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART A: EIGENVALUE SPECTRUM")
    pr(f"{'═' * 72}\n")
    pr("  Diaconis-Fulman (ADDITION in base b): eigenvalues are 1, 1/b, 1/b², ...")
    pr("  Our chain (MULTIPLICATION in base b): eigenvalues = ?\n")

    for b in [2, 3, 5, 7]:
        pd = product_distribution(b)
        pr(f"  Base b={b}:")
        pr(f"    Product dist g·h: support {np.where(pd > 0)[0].tolist()}")
        pr(f"    E[g·h] = {np.sum(np.arange(len(pd)) * pd):.4f} "
           f"(= (b-1)²/4 = {(b-1)**2/4:.4f})")

        for n in [4, 8, 16, 32, 64]:
            conv = nfold_convolution(pd, n)
            T, n_states = build_transition_matrix(b, conv)

            try:
                eigs = np.linalg.eigvals(T)
            except Exception:
                continue
            eigs_sorted = sorted(eigs.real, reverse=True)

            df_eigs = [1 / b**k for k in range(min(6, n_states))]

            pr(f"    n={n:3d} ({n_states:3d} states): "
               f"λ = [{', '.join(f'{e:.4f}' for e in eigs_sorted[:6])}]")
            if n == 4:
                pr(f"           D-F reference: "
                   f"[{', '.join(f'{e:.4f}' for e in df_eigs)}]")
        pr()

    # ═══════════════════════════════════════════════════════════════
    # PART B: STATIONARY DISTRIBUTION
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART B: STATIONARY DISTRIBUTION")
    pr(f"{'═' * 72}\n")

    for b in [2, 3, 5]:
        pd = product_distribution(b)
        pr(f"  Base b={b}:")

        for n in [16, 32, 64]:
            conv = nfold_convolution(pd, n)
            T, ns = build_transition_matrix(b, conv)

            eigs, vecs = np.linalg.eig(T)
            idx_1 = np.argmin(np.abs(eigs - 1.0))
            stat = np.abs(vecs[:, idx_1])
            stat /= stat.sum()

            mu = n * (b - 1)**2 / (4 * b)
            mean_carry = np.sum(np.arange(ns) * stat)
            var_carry = np.sum(np.arange(ns)**2 * stat) - mean_carry**2

            pr(f"    n={n:3d}: <c> = {mean_carry:.3f} "
               f"(theory bulk ≈ n(b-1)²/(4b) = {mu:.3f}), "
               f"Var = {var_carry:.3f}, "
               f"states={ns}")
        pr()

    # ═══════════════════════════════════════════════════════════════
    # PART C: SPECTRAL GAP AND MIXING TIME
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART C: SPECTRAL GAP AND MIXING TIME")
    pr(f"{'═' * 72}\n")
    pr("  Spectral gap δ = 1 - |λ₂| determines mixing time ~ 1/δ.\n")

    for b in [2, 3, 5, 7, 11]:
        pd = product_distribution(b)
        pr(f"  Base b={b}:")

        for n in [8, 16, 32, 64]:
            conv = nfold_convolution(pd, n)
            T, ns = build_transition_matrix(b, conv)

            eigs = np.linalg.eigvals(T)
            eigs_abs = sorted(np.abs(eigs), reverse=True)
            lambda2 = eigs_abs[1] if len(eigs_abs) > 1 else 0
            gap = 1 - lambda2
            mix_time = 1 / max(gap, 1e-10)

            pr(f"    n={n:3d}: |λ₂|={lambda2:.6f}, gap={gap:.6f}, "
               f"t_mix≈{mix_time:.1f}")
        pr()

    # ═══════════════════════════════════════════════════════════════
    # PART D: COMPARISON WITH EMPIRICAL CARRY DISTRIBUTION
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART D: MARKOV STATIONARY vs EMPIRICAL CARRIES")
    pr(f"{'═' * 72}\n")

    import random
    from collections import Counter
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
    from carry_utils import random_prime, to_digits

    random.seed(42)

    def empirical_carry_dist(bits, base, n_samples):
        counts = Counter()
        total = 0
        for _ in range(n_samples):
            p = random_prime(bits)
            q = random_prime(bits)
            gd = to_digits(p, base)
            hd = to_digits(q, base)
            D = len(to_digits(p * q, base))
            conv = [0] * (len(gd) + len(hd) - 1)
            for i, a in enumerate(gd):
                for j, b_val in enumerate(hd):
                    conv[i + j] += a * b_val
            carries = [0] * (D + 1)
            for i in range(D):
                v = conv[i] if i < len(conv) else 0
                carries[i + 1] = (v + carries[i]) // base
            for c in carries[2:-2]:
                counts[c] += 1
                total += 1
        return {c: n / total for c, n in counts.items()}, total

    for base in [2, 3, 5]:
        bits = 32
        pr(f"  Base {base}, {bits}-bit semiprimes (1000 samples):")

        emp_dist, n_total = empirical_carry_dist(bits, base, 1000)

        n_pairs = 2 * bits // max(1, int(math.log2(base)))
        pd = product_distribution(base)
        conv_dist = nfold_convolution(pd, n_pairs)
        T, ns = build_transition_matrix(base, conv_dist)
        eigs, vecs = np.linalg.eig(T)
        idx_1 = np.argmin(np.abs(eigs - 1.0))
        stat = np.abs(vecs[:, idx_1])
        stat /= stat.sum()

        max_c = max(max(emp_dist.keys()), ns - 1)
        pr(f"    {'carry':>6s}  {'empirical':>10s}  {'Markov':>10s}  {'ratio':>8s}")
        for c in range(min(max_c + 1, 12)):
            e = emp_dist.get(c, 0)
            m = stat[c] if c < ns else 0
            r = e / m if m > 1e-10 else float('inf')
            pr(f"    {c:6d}  {e:10.5f}  {m:10.5f}  {r:8.3f}")
        pr()

    # ═══════════════════════════════════════════════════════════════
    # PART E: FULL EIGENVALUE LIST FOR BASE 2
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART E: FULL EIGENVALUE SPECTRUM — BASE 2")
    pr(f"{'═' * 72}\n")

    pd2 = product_distribution(2)
    for n in [8, 16, 32, 64, 128]:
        conv = nfold_convolution(pd2, n)
        T, ns = build_transition_matrix(2, conv)
        eigs = np.sort(np.linalg.eigvals(T).real)[::-1]
        pr(f"  n={n:4d} ({ns:4d} states):")
        pr(f"    All eigenvalues: {np.array2string(eigs[:20], precision=6, separator=', ')}")
        if ns > 20:
            pr(f"    ... ({ns-20} more)")

        ratios = [eigs[i] / eigs[i + 1] if abs(eigs[i + 1]) > 1e-10 else float('inf')
                  for i in range(min(5, ns - 1))]
        pr(f"    Ratios λ_k/λ_{'{k+1}'}: {[f'{r:.3f}' for r in ratios]}")
        pr()

    # ═══════════════════════════════════════════════════════════════
    # SYNTHESIS
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("SYNTHESIS")
    pr(f"{'═' * 72}")
    pr("""
  Key findings:
    1. Eigenvalues of the MULTIPLICATION carry chain:
       Are they 1, 1/b, 1/b², ... (like Diaconis-Fulman addition)?
       Or something different?

    2. Spectral gap: how fast does the chain mix?
       Fast mixing → bulk carries reach stationarity quickly
       → boundary layer (top few positions) dominates α_k corrections

    3. Stationary distribution: does it match empirical carries?
       If yes: Markov model is correct for the bulk.
       prior experiments already showed this; we confirm with exact eigenanalysis.

    4. Eigenvalue structure as n → ∞:
       If eigenvalues cluster or develop structure, this may connect
       to the transfer operator formalism .
""")
    pr(f"\n  Total runtime: {time.time() - t0:.1f}s")
    pr("=" * 72)


if __name__ == '__main__':
    main()
