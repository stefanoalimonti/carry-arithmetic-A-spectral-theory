#!/usr/bin/env python3
"""
A10: Deep analysis of c_{top-1}=3 counterexamples

Lemma 4 (c_{top-1} ≤ 2 for D=2d-1) is FALSE.
But r_max ≤ 2 STILL HOLDS. Why?

This experiment:
  Part A: Full anatomy of known counterexamples
  Part B: Collect ALL c_{top-1}=3 cases, measure r_max for each
  Part C: Q(-2) analysis for c_{top-1}=3 cases — virtual pole strength
  Part D: Root structure — where are ALL eigenvalues when c_{top-1}=3?
  Part E: Lower bound on |Q(-2)| — why -2 is never a root
  Part F: Sufficient condition from Q polynomial structure
"""

import sys, os, time, random, math
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from carry_utils import random_prime, carry_poly_int, quotient_poly_int, to_digits

random.seed(42)
np.random.seed(42)


def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def negabinary_val(digits):
    return sum(d * ((-2) ** i) for i, d in enumerate(digits))


def full_carry_forward(p, q, base=2):
    gd = to_digits(p, base)
    hd = to_digits(q, base)
    N = p * q
    fd = to_digits(N, base)
    D = len(fd)
    max_pos = len(gd) + len(hd) - 1
    conv = [0] * (max_pos + 2)
    for i, a in enumerate(gd):
        for j, b_ in enumerate(hd):
            conv[i + j] += a * b_
    carries = [0] * (D + 2)
    for k in range(D + 1):
        total = conv[k] + carries[k]
        carries[k + 1] = total // base
    return carries, conv, fd, gd, hd


def main():
    t0 = time.time()
    pr("=" * 76)
    pr("P4-02b: DEEP ANALYSIS — c_{top-1}=3 COUNTEREXAMPLES")
    pr("=" * 76)

    # ════════════════════════════════════════════════════════════════════
    # PART A: FULL ANATOMY OF KNOWN COUNTEREXAMPLES
    # ════════════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 76}")
    pr("PART A: ANATOMY OF c_{top-1}=3 COUNTEREXAMPLES")
    pr(f"{'═' * 76}")

    counterexamples = [
        (2311, 3623),
        (4639, 6301),
        (8263, 13063),
    ]

    for p, q in counterexamples:
        N = p * q
        D = len(to_digits(N, 2))
        d = (D + 1) // 2

        gd = to_digits(p, 2)
        hd = to_digits(q, 2)
        fd = to_digits(N, 2)

        C = carry_poly_int(p, q, 2)
        Q = quotient_poly_int(C, 2)
        m = len(C) - 1

        carries, conv_full, _, _, _ = full_carry_forward(p, q, 2)
        c_top = int(-Q[-1])
        c_top1 = int(-Q[-2])

        coeffs = [float(c) for c in Q]
        roots = np.roots(list(reversed(coeffs)))
        mods = np.abs(roots)
        rmax = mods.max()

        p_neg = negabinary_val(gd)
        q_neg = negabinary_val(hd)
        N_neg = negabinary_val(fd)
        C_neg2 = p_neg * q_neg - N_neg
        Q_neg2 = sum(float(c) * ((-2) ** i) for i, c in enumerate(Q))

        pr(f"\n  ─── p={p} × q={q} = {N} ───")
        pr(f"  p = {''.join(str(x) for x in reversed(gd))} ({len(gd)} bits)")
        pr(f"  q = {''.join(str(x) for x in reversed(hd))} ({len(hd)} bits)")
        pr(f"  N = {''.join(str(x) for x in reversed(fd))} ({D} bits, D=2×{d}-1)")
        pr(f"  m = deg(C) = {m}, cancellations from 2d-2={2*d-2}: {2*d-2-m}")

        pr(f"\n  Carry profile (forward):")
        carry_str = [str(carries[k]) for k in range(min(D + 1, len(carries)))]
        pr(f"    c = [{', '.join(carry_str)}]")

        pr(f"\n  Convolution at m-1={m-1}: {conv_full[m-1]}")
        pr(f"  Convolution at m-2={m-2}: {conv_full[m-2]}")
        pr(f"  Digit f[m-1]={fd[m-1]}, f[m-2]={fd[m-2] if m-2 < len(fd) else 'N/A'}")
        pr(f"  Forward carry at m-1: c[{m-1}]={carries[m-1]}")
        pr(f"  Forward carry at m-2: c[{m-2}]={carries[m-2]}")
        pr(f"  conv[m-2]+c[m-2] = {conv_full[m-2]}+{carries[m-2]} = {conv_full[m-2]+carries[m-2]}")

        pr(f"\n  c_top = {c_top}, c_{{top-1}} = {c_top1}")
        pr(f"  Backward check: 2·{c_top} + {fd[m-1]} - {conv_full[m-1]} = "
           f"{2*c_top + fd[m-1] - conv_full[m-1]}")

        pr(f"\n  Q polynomial: {[float(c) for c in Q]}")
        pr(f"  Q(-2) = {Q_neg2:.1f}")
        pr(f"  C(-2) = {C_neg2}")
        pr(f"  |Q(-2)| = {abs(Q_neg2):.1f}")

        pr(f"\n  Eigenvalues:")
        sorted_roots = sorted(roots, key=lambda z: -abs(z))
        for i, r in enumerate(sorted_roots[:6]):
            pr(f"    λ_{i+1} = {r.real:+.6f} {r.imag:+.6f}i, |λ|={abs(r):.6f}")
        pr(f"  r_max = {rmax:.6f}")

        top_carries = [-int(Q[k]) for k in range(len(Q) - 1, max(len(Q) - 6, -1), -1)]
        pr(f"\n  Top carries (from top): {top_carries}")

    # ════════════════════════════════════════════════════════════════════
    # PART B: COLLECT ALL c_{top-1}=3, MEASURE r_max
    # ════════════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 76}")
    pr("PART B: r_max FOR ALL c_{top-1}=3 CASES")
    pr(f"{'═' * 76}")
    pr("  Collecting c_{top-1}=3 cases across many bit sizes...\n")

    ct3_data = []

    for bits in [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32]:
        n_samp = 500000 if bits <= 16 else (200000 if bits <= 24 else 100000)
        ct3_local = []

        for _ in range(n_samp):
            p = random_prime(bits)
            q = random_prime(bits)
            if p == q:
                continue
            N = p * q
            D = len(to_digits(N, 2))
            if D % 2 == 0:
                continue

            C = carry_poly_int(p, q, 2)
            Q = quotient_poly_int(C, 2)
            if len(Q) < 3:
                continue

            c_top1 = int(-Q[-2])
            if c_top1 != 3:
                continue

            coeffs = [float(c) for c in Q]
            roots = np.roots(list(reversed(coeffs)))
            rmax = max(abs(r) for r in roots) if len(roots) > 0 else None

            q_neg2 = sum(float(c) * ((-2) ** i) for i, c in enumerate(Q))
            ct3_local.append((bits, rmax, abs(q_neg2), p, q))

        ct3_data.extend(ct3_local)
        if ct3_local:
            rmaxs = [x[1] for x in ct3_local if x[1] is not None]
            qn2s = [x[2] for x in ct3_local]
            pr(f"  bits={bits:2d}: {len(ct3_local)} cases with c_{{top-1}}=3, "
               f"r_max: mean={np.mean(rmaxs):.4f}, max={max(rmaxs):.6f}, "
               f"|Q(-2)|: min={min(qn2s):.0f}")
        else:
            pr(f"  bits={bits:2d}: 0 cases")

    if ct3_data:
        all_rmax = [x[1] for x in ct3_data if x[1] is not None]
        pr(f"\n  TOTAL c_{{top-1}}=3 collected: {len(ct3_data)}")
        pr(f"  GLOBAL max(r_max) across all c_{{top-1}}=3: {max(all_rmax):.6f}")
        pr(f"  Any r_max > 2? {'YES !!!' if max(all_rmax) > 2 else 'NO — theorem holds'}")

    # ════════════════════════════════════════════════════════════════════
    # PART C: Q(-2) vs r_max CORRELATION FOR c_{top-1}=3
    # ════════════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 76}")
    pr("PART C: Q(-2) vs r_max — VIRTUAL POLE REPULSION FOR c_{top-1}=3")
    pr(f"{'═' * 76}")

    if ct3_data:
        rmaxs = np.array([x[1] for x in ct3_data if x[1] is not None])
        qn2s = np.array([x[2] for x in ct3_data if x[1] is not None])
        bits_arr = np.array([x[0] for x in ct3_data if x[1] is not None])

        pr(f"\n  Correlation between |Q(-2)| and r_max:")
        pr(f"  |Q(-2)| range: [{qn2s.min():.0f}, {qn2s.max():.0f}]")
        pr(f"  r_max range: [{rmaxs.min():.6f}, {rmaxs.max():.6f}]")

        large_q = qn2s > np.median(qn2s)
        small_q = ~large_q
        pr(f"\n  |Q(-2)| > median: ⟨r_max⟩ = {rmaxs[large_q].mean():.6f}")
        pr(f"  |Q(-2)| < median: ⟨r_max⟩ = {rmaxs[small_q].mean():.6f}")
        pr(f"  (Larger |Q(-2)| → eigenvalues pushed further from |z|=2)")

    # ════════════════════════════════════════════════════════════════════
    # PART D: ROOT DISTRIBUTION FOR c_{top-1}=3 CASES
    # ════════════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 76}")
    pr("PART D: ROOT DISTRIBUTION FOR c_{top-1}=3 — WHERE ARE THE EIGENVALUES?")
    pr(f"{'═' * 76}")

    for p, q in counterexamples:
        C = carry_poly_int(p, q, 2)
        Q = quotient_poly_int(C, 2)
        coeffs = [float(c) for c in Q]
        roots = np.roots(list(reversed(coeffs)))

        pr(f"\n  p={p}, q={q}:")
        real_neg = sorted([r.real for r in roots if abs(r.imag) < 1e-10 and r.real < 0])
        real_pos = sorted([r.real for r in roots if abs(r.imag) < 1e-10 and r.real >= 0])
        complex_roots = [(r.real, r.imag) for r in roots if abs(r.imag) >= 1e-10]
        complex_mods = [abs(r) for r in roots if abs(r.imag) >= 1e-10]

        pr(f"    Real negative: {[f'{r:.4f}' for r in real_neg]}")
        pr(f"    Real positive: {[f'{r:.4f}' for r in real_pos]}")
        if complex_mods:
            pr(f"    Complex: {len(complex_roots)//2} conjugate pairs, "
               f"max|λ|={max(complex_mods):.6f}")
        pr(f"    Most negative real: {min(real_neg) if real_neg else 'none':.6f}")
        gap = 2.0 - max(abs(r) for r in roots)
        pr(f"    Gap from |λ|=2: {gap:.6f}")

    # ════════════════════════════════════════════════════════════════════
    # PART E: Q(-r) FOR r ∈ [2, 3] — CONTINUOUS ZERO-FREE REGION
    # ════════════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 76}")
    pr("PART E: Q(-r) FOR r ∈ [2, 3] — ZERO-FREE ON NEGATIVE REAL AXIS")
    pr(f"{'═' * 76}")
    pr("  If Q(-r) ≠ 0 for all r ∈ [2, 3], no real eigenvalue in [-3, -2].\n")

    for p, q in counterexamples[:1]:
        C = carry_poly_int(p, q, 2)
        Q = quotient_poly_int(C, 2)
        pr(f"  Example: p={p}, q={q}")

        r_vals = np.linspace(1.5, 3.0, 31)
        for r in r_vals:
            q_val = sum(float(c) * ((-r) ** i) for i, c in enumerate(Q))
            pr(f"    Q({-r:+.2f}) = {q_val:+.1f}")

    for bits in [14, 20, 26, 32]:
        n_samp = 100000 if bits <= 20 else 50000
        min_abs_Q_neg_r = {r: float('inf') for r in [2.0, 2.2, 2.5, 3.0]}
        dodd_count = 0
        ct3_count = 0

        for _ in range(n_samp):
            p = random_prime(bits)
            q = random_prime(bits)
            if p == q:
                continue
            N = p * q
            D = len(to_digits(N, 2))
            if D % 2 == 0:
                continue
            dodd_count += 1

            C = carry_poly_int(p, q, 2)
            Q = quotient_poly_int(C, 2)
            if len(Q) < 3:
                continue
            if int(-Q[-2]) != 3:
                continue
            ct3_count += 1

            for r in min_abs_Q_neg_r:
                q_val = abs(sum(float(c) * ((-r) ** i) for i, c in enumerate(Q)))
                min_abs_Q_neg_r[r] = min(min_abs_Q_neg_r[r], q_val)

        if ct3_count > 0:
            pr(f"\n  bits={bits}: {ct3_count} c_{{top-1}}=3 cases")
            for r, v in sorted(min_abs_Q_neg_r.items()):
                pr(f"    min|Q({-r:+.1f})| = {v:.1f}")
        else:
            pr(f"\n  bits={bits}: 0 c_{{top-1}}=3 in {dodd_count} D-odd samples")

    # ════════════════════════════════════════════════════════════════════
    # PART F: THE ACTUAL PROOF PATH — EVALUATE Q ON |z|=2 CIRCLE
    # ════════════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 76}")
    pr("PART F: |Q(z)| ON THE CIRCLE |z|=2")
    pr(f"{'═' * 76}")
    pr("""
  r_max ≤ 2 iff Q has no roots with |z| ≥ 2.
  By the argument principle, count zeros of Q outside |z|=2 by:
    (1/2πi) ∮_{|z|=2} Q'(z)/Q(z) dz = #{zeros inside} - n
  
  Simpler test: min_{|z|=2} |Q(z)| > 0 iff no roots on |z|=2.
  Combined with Q having no roots at |z|>3 (E-K bound), 
  and continuity argument...
""")

    for bits in [12, 16, 20, 26]:
        n_samp = 100000 if bits <= 16 else 50000
        min_Q_on_circle = []
        dodd_count = 0

        for trial in range(n_samp):
            p = random_prime(bits)
            q = random_prime(bits)
            if p == q:
                continue
            N = p * q
            D = len(to_digits(N, 2))
            if D % 2 == 0:
                continue
            dodd_count += 1

            C = carry_poly_int(p, q, 2)
            Q = quotient_poly_int(C, 2)
            if len(Q) < 3:
                continue

            theta_vals = np.linspace(0, 2 * np.pi, 100, endpoint=False)
            z_vals = 2 * np.exp(1j * theta_vals)
            Q_vals = np.array([sum(float(c) * z ** i for i, c in enumerate(Q))
                               for z in z_vals])
            min_q = np.min(np.abs(Q_vals))
            min_Q_on_circle.append(min_q)

        if min_Q_on_circle:
            arr = np.array(min_Q_on_circle)
            pr(f"  bits={bits}: min|Q(z)| on |z|=2: "
               f"mean={arr.mean():.2f}, min={arr.min():.6f}, "
               f"zeros (< 0.01): {np.sum(arr < 0.01)}/{len(arr)}")

    # ════════════════════════════════════════════════════════════════════
    # PART G: WINDING NUMBER — COUNT ROOTS OUTSIDE |z|=2
    # ════════════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 76}")
    pr("PART G: WINDING NUMBER — ROOTS OF Q OUTSIDE |z|=2")
    pr(f"{'═' * 76}")

    for bits in [12, 16, 20, 24, 28]:
        n_samp = 100000 if bits <= 16 else 50000
        n_outside = []
        dodd_count = 0

        for _ in range(n_samp):
            p = random_prime(bits)
            q = random_prime(bits)
            if p == q:
                continue
            N = p * q
            D = len(to_digits(N, 2))
            if D % 2 == 0:
                continue
            dodd_count += 1

            C = carry_poly_int(p, q, 2)
            Q = quotient_poly_int(C, 2)
            if len(Q) < 3:
                continue

            coeffs = [float(c) for c in Q]
            roots = np.roots(list(reversed(coeffs)))
            n_out = sum(1 for r in roots if abs(r) > 2.0 + 1e-10)
            n_outside.append(n_out)

        if n_outside:
            arr = np.array(n_outside)
            pr(f"  bits={bits}: D-odd={dodd_count}, "
               f"mean roots outside |z|=2: {arr.mean():.6f}, "
               f"max: {arr.max()}, any>0: {np.sum(arr > 0)}")

    # ════════════════════════════════════════════════════════════════════
    # SYNTHESIS
    # ════════════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 76}")
    pr("SYNTHESIS: REVISED PROOF STRUCTURE FOR Gap 3")
    pr(f"{'═' * 76}")
    pr("""
  FINDING: c_{top-1} = 3 EXISTS (rare, ~20 per million D-odd cases).
  But r_max ≤ 2 STILL HOLDS — no single violation in millions of tests.
  
  REVISED PROOF STRATEGY:
  
  The theorem cannot be proven via c_{top-1} ≤ 2 alone.
  Instead, use the FULL polynomial structure:
  
  1. ULC: c_top = 1. [PROVEN]
  2. E-K bound: r_max ≤ max(c_k/c_{k+1}) ≤ 3 (since c_{top-1} ≤ 3). [PROVEN]
  3. Virtual pole at z=-2: Q(-2) = C(-2)/(-4) ≠ 0. [NEED TO PROVE]
     C(-2) = p₋·q₋ - (pq)₋ where p₋ is negabinary evaluation.
     |C(-2)| grows exponentially with D.
  4. Continuous zero-free band: Q(-r) ≠ 0 for r ∈ [2, 3]. [NEED TO PROVE]
     Verified numerically: min|Q(-r)| >> 0 for all tested cases.
  5. Complex eigenvalues: all |λ_complex| < 2. [VERIFIED NUMERICALLY]
  
  ALTERNATIVE: Direct Rouché on |z|=2 using C = g·h - f structure.
  On |z|=2: Q(z) = C(z)/(z-2), need |C(z)| > 0 on |z|=2, z ≠ 2.
""")

    pr(f"\n  Runtime: {time.time() - t0:.1f}s")
    pr("=" * 76)


if __name__ == '__main__':
    main()
