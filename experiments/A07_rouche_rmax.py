#!/usr/bin/env python3
"""
A07: Rouché Analysis at |z| = 2

Goal: Determine whether Rouché's theorem can prove r_max ≤ 2 for all
carry quotient polynomials.

The monic quotient polynomial (negated Q):
  p(z) = z^{D-1} + c_{D-1} z^{D-2} + ... + c_2 z + c_1

where c_k = carry_k ≥ 0. For r_max ≤ 2 via Rouché, we need:

  |z^{D-1}| > |c_{D-1} z^{D-2} + ... + c_1|   on |z| = 2

i.e.  2^{D-1} > Σ_{k=1}^{D-1} c_k · 2^{k-1}

Equivalently: Σ c_k · 2^{k-1} < 2^{D-1}, i.e. Σ c_k / 2^{D-k} < 1.

This is equivalent to: the carry polynomial evaluated at z=2, minus the
leading term, divided by 2^{D-1}, is less than 1.

Parts:
  A) Compute the Rouché ratio R = (Σ c_k · 2^{k-1}) / 2^{D-1} for many semiprimes
  B) Find worst cases: when R is closest to 1
  C) Analytic bound on R using carry statistics
  D) Refined Rouché on |z| = 2+ε for ε → 0
  E) Pellet's theorem: separate zero-counting on annuli
"""

import sys, os, time, random, math
import numpy as np
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from carry_utils import random_prime, to_digits, carry_poly_int, quotient_poly_int


random.seed(42)
np.random.seed(42)


def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def compute_carries(p, q, base=2):
    gd = to_digits(p, base)
    hd = to_digits(q, base)
    conv_len = len(gd) + len(hd) - 1
    conv = [0] * conv_len
    for i, a in enumerate(gd):
        for j, b_val in enumerate(hd):
            conv[i + j] += a * b_val
    max_len = max(conv_len, len(to_digits(p * q, base))) + 2
    carries = [0] * (max_len + 1)
    for k in range(max_len):
        conv_k = conv[k] if k < conv_len else 0
        carries[k + 1] = (conv_k + carries[k]) // base
    last_nz = 0
    for k in range(max_len, 0, -1):
        if carries[k] != 0:
            last_nz = k
            break
    return carries[:last_nz + 1]


def rouche_ratio(carries):
    """Compute R = (Σ_{k=1}^{D-1} c_k · 2^{k-1}) / 2^{D-1}."""
    D = len(carries) - 1
    if D < 3:
        return None
    total = sum(carries[k] * (2 ** (k - 1)) for k in range(1, D))
    return total / (2.0 ** (D - 1))


def rouche_ratio_exact(carries):
    """Exact integer ratio: Σ c_k · 2^{k-1} vs 2^{D-1}."""
    D = len(carries) - 1
    if D < 3:
        return None, None
    numer = sum(carries[k] * (2 ** (k - 1)) for k in range(1, D))
    denom = 2 ** (D - 1)
    return numer, denom


def main():
    t0 = time.time()
    pr("=" * 72)
    pr("P1-01: ROUCHÉ ANALYSIS AT |z| = 2")
    pr("=" * 72)

    # ═══════════════════════════════════════════════════════════════
    # PART A: ROUCHÉ RATIO DISTRIBUTION
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART A: ROUCHÉ RATIO R = Σ c_k·2^{k-1} / 2^{D-1}")
    pr(f"{'═' * 72}")
    pr("If R < 1 for ALL semiprimes, Rouché proves r_max ≤ 2.\n")

    results_by_bits = {}

    for bits in [8, 12, 16, 20, 24, 32]:
        n_test = 30000
        ratios = []
        rouche_fails = 0

        for _ in range(n_test):
            p = random_prime(bits)
            q = random_prime(bits)
            if p == q:
                continue
            carries = compute_carries(p, q, 2)
            R = rouche_ratio(carries)
            if R is None:
                continue
            ratios.append(R)
            if R >= 1.0:
                rouche_fails += 1

        arr = np.array(ratios)
        results_by_bits[bits] = arr
        pr(f"  {bits:3d}-bit ({len(ratios)} samples):")
        pr(f"    R < 1 (Rouché OK): {len(ratios) - rouche_fails} "
           f"({100*(1 - rouche_fails/len(ratios)):.2f}%)")
        pr(f"    R ≥ 1 (Rouché FAILS): {rouche_fails}")
        pr(f"    max(R)  = {arr.max():.6f}")
        pr(f"    mean(R) = {arr.mean():.6f}")
        pr(f"    med(R)  = {np.median(arr):.6f}")

    # ═══════════════════════════════════════════════════════════════
    # PART B: WORST CASES — WHEN R IS LARGEST
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART B: WORST CASES — LARGEST ROUCHÉ RATIOS")
    pr(f"{'═' * 72}\n")

    worst_cases = []
    for _ in range(100000):
        bits = random.randint(8, 32)
        p = random_prime(bits)
        q = random_prime(bits)
        if p == q:
            continue
        carries = compute_carries(p, q, 2)
        R = rouche_ratio(carries)
        if R is None:
            continue
        D = len(carries) - 1
        worst_cases.append((R, p, q, D, carries[:]))

    worst_cases.sort(key=lambda x: -x[0])
    pr(f"  {'R':>10s} | {'p':>10s} | {'q':>10s} | {'D':>4s} | carry profile")
    pr(f"  {'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*4}-+-{'-'*30}")
    for R, p, q, D, carries in worst_cases[:15]:
        nonzero = [(k, carries[k]) for k in range(1, len(carries)) if carries[k] > 0]
        tag = " ** FAILS" if R >= 1.0 else ""
        pr(f"  {R:10.6f} | {p:10d} | {q:10d} | {D:4d} | {nonzero[-4:]}{tag}")

    n_fail = sum(1 for R, _, _, _, _ in worst_cases if R >= 1.0)
    pr(f"\n  Total Rouché failures: {n_fail}/{len(worst_cases)} "
       f"({100*n_fail/len(worst_cases):.3f}%)")

    # ═══════════════════════════════════════════════════════════════
    # PART C: ANALYTIC BOUND ON R USING CARRY PROFILE
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART C: ANALYTIC BOUND ON R")
    pr(f"{'═' * 72}")
    pr("""
  The carry c_k satisfies: c_k ≤ (b-1)·min(k, D-k) for base b=2.
  In practice, c_k → (b-1)/(2b) = 1/4 in the interior (Markov stationarity).

  Upper bound on R:
    R = Σ_{k=1}^{D-1} c_k / 2^{D-k}
      ≤ c_{D-1}/2 + c_{D-2}/4 + ... + c_1/2^{D-1}

  If all c_k ≤ C_max, then R ≤ C_max · Σ_{j=1}^{D-1} 2^{-j} < C_max.
  For base 2: c_k ∈ {0, 1, ...}, but boundary carries can be 2 or more.
""")

    max_carry_by_position = {}
    for _ in range(100000):
        bits = random.randint(8, 32)
        p = random_prime(bits)
        q = random_prime(bits)
        carries = compute_carries(p, q, 2)
        D = len(carries) - 1
        for k in range(1, D + 1):
            pos_from_top = D - k
            if pos_from_top not in max_carry_by_position:
                max_carry_by_position[pos_from_top] = 0
            max_carry_by_position[pos_from_top] = max(
                max_carry_by_position[pos_from_top], carries[k])

    pr("  Max carry by position from top (D-k):")
    for pos in sorted(max_carry_by_position.keys())[:15]:
        pr(f"    position D-{pos}: max carry = {max_carry_by_position[pos]}")

    # ═══════════════════════════════════════════════════════════════
    # PART D: REFINED ROUCHÉ ON |z| = 2+ε
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART D: ROUCHÉ ON |z| = 2+ε")
    pr(f"{'═' * 72}")
    pr("If R(2+ε) < 1 for some ε > 0, then r_max < 2+ε.\n")

    for epsilon in [0.01, 0.05, 0.1, 0.5, 1.0]:
        r = 2.0 + epsilon
        n_test = 50000
        fails = 0
        max_ratio = 0.0

        for _ in range(n_test):
            bits = random.randint(8, 32)
            p = random_prime(bits)
            q = random_prime(bits)
            carries = compute_carries(p, q, 2)
            D = len(carries) - 1
            if D < 3:
                continue
            lower = sum(carries[k] * r ** (k - 1) for k in range(1, D))
            leading = r ** (D - 1)
            ratio = lower / leading
            max_ratio = max(max_ratio, ratio)
            if ratio >= 1.0:
                fails += 1

        pr(f"  ε = {epsilon:.2f} (|z| = {r:.2f}): "
           f"fails = {fails}/{n_test}, max ratio = {max_ratio:.6f}")

    # ═══════════════════════════════════════════════════════════════
    # PART E: PELLET'S THEOREM — ZERO COUNTING IN ANNULI
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART E: EIGENVALUE DISTRIBUTION IN ANNULI")
    pr(f"{'═' * 72}\n")

    annuli_counts = {r: [] for r in [0.5, 1.0, 1.5, 2.0, 2.5]}

    for _ in range(20000):
        bits = random.randint(12, 32)
        p = random_prime(bits)
        q = random_prime(bits)
        C = carry_poly_int(p, q, 2)
        Q = quotient_poly_int(C, 2)
        if len(Q) < 3:
            continue
        lead = float(Q[-1])
        if abs(lead) < 1e-30:
            continue
        n = len(Q) - 1
        M = np.zeros((n, n))
        for i in range(n - 1):
            M[i + 1, i] = 1.0
        for i in range(n):
            M[i, n - 1] = -float(Q[i]) / lead
        try:
            ev = np.linalg.eigvals(M)
        except Exception:
            continue

        mods = np.abs(ev)
        for r in annuli_counts:
            annuli_counts[r].append(int(np.sum(mods > r)))

    pr(f"  Fraction of eigenvalues with |λ| > r:")
    for r in sorted(annuli_counts.keys()):
        arr = np.array(annuli_counts[r])
        pr(f"    |λ| > {r:.1f}: mean = {arr.mean():.3f}, "
           f"max = {arr.max()}, P(any) = {np.mean(arr > 0):.4f}")

    # ═══════════════════════════════════════════════════════════════
    # SYNTHESIS
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("SYNTHESIS — ROUCHÉ ANALYSIS")
    pr(f"{'═' * 72}")

    any_fail = any(
        1 for R, _, _, _, _ in worst_cases if R >= 1.0
    )

    if any_fail:
        pr("""
  RESULT: Rouché at |z| = 2 FAILS for some semiprimes.
  This means the leading term z^{D-1} does NOT always dominate on |z| = 2.
  However, the actual eigenvalues may still be ≤ 2 — Rouché is sufficient
  but not necessary. A different proof technique is needed.
""")
    else:
        pr("""
  RESULT: Rouché at |z| = 2 SUCCEEDS for all tested semiprimes!
  This means: 2^{D-1} > Σ c_k · 2^{k-1} for all carry profiles tested.

  To make this a theorem, we need to prove the bound analytically using
  the carry recursion constraints.
""")

    pr(f"  Runtime: {time.time() - t0:.1f}s")
    pr("=" * 72)


if __name__ == '__main__':
    main()
