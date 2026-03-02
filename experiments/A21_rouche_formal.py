#!/usr/bin/env python3
"""
A21: Formalize the Rouché gap analysis on |z|=b for Conjecture 10.

For z = b*exp(i*theta), decompose |g(z)h(z)| - |n(z)| analytically.
Key insight: g(b)*h(b) = n(b) (the root at z=b), so the gap vanishes
at theta=0. But the *derivative* at theta=0 determines the separation rate.

This experiment:
  (A) Compute d/dtheta [|g(be^{itheta})h(be^{itheta})|/|n(be^{itheta})|] at theta=0
  (B) Decompose the gap into Fourier components
  (C) Prove |C(be^{itheta})| > 0 for theta != 0 analytically for structured profiles
  (D) Test the "minimum gap growth" conjecture: gap_min ~ 2^{2K}
"""

import sys
import os
import random
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from carry_utils import random_prime, to_digits

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


def eval_poly(coeffs, z):
    """Evaluate polynomial with little-endian coefficients at z."""
    return sum(c * z**k for k, c in enumerate(coeffs))


def main():
    pr("=" * 72)
    pr("  A21: ROUCHÉ GAP FORMALIZATION ON |z| = b")
    pr("=" * 72)

    base = 2

    # ═══════════════════════════════════════════════════════════════
    # PART A: DERIVATIVE ANALYSIS AT theta=0
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART A: DERIVATIVE OF |C(be^{itheta})| AT theta=0")
    pr(f"{'═' * 72}")
    pr()
    pr("C(z) = g(z)h(z) - f(z), so C(b) = 0.")
    pr("Near theta=0: |C(be^{itheta})| ~ |C'(b)| * |1 - e^{itheta}| * b")
    pr("= |C'(b)| * 2b * |sin(theta/2)| ~ |C'(b)| * b * |theta|")
    pr()

    derivative_data = []

    for bits in [8, 12, 16, 20]:
        n_samples = 5000 if bits <= 16 else 2000
        C_prime_vals = []
        min_gaps = []

        for _ in range(n_samples):
            p = random_prime(bits)
            q = random_prime(bits)
            if p == q:
                continue
            N = p * q
            gd = to_digits(p, base)
            hd = to_digits(q, base)
            fd = to_digits(N, base)

            n_theta = 2048
            theta = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
            z = base * np.exp(1j * theta)

            g_vals = sum(c * z**k for k, c in enumerate(gd))
            h_vals = sum(c * z**k for k, c in enumerate(hd))
            f_vals = sum(c * z**k for k, c in enumerate(fd))

            C_vals = g_vals * h_vals - f_vals

            exclude = max(1, n_theta // 64)
            gap = np.abs(C_vals)
            gap_nonzero = np.concatenate([gap[exclude:-exclude]])
            if len(gap_nonzero) > 0:
                min_gap = gap_nonzero.min()
                min_gaps.append(min_gap)

            g_coeffs = np.array(gd, dtype=float)
            h_coeffs = np.array(hd, dtype=float)
            f_coeffs = np.array(fd, dtype=float)

            g_prime_b = sum(k * c * base**(k-1) for k, c in enumerate(gd) if k > 0)
            h_prime_b = sum(k * c * base**(k-1) for k, c in enumerate(hd) if k > 0)
            g_b = sum(c * base**k for k, c in enumerate(gd))
            h_b = sum(c * base**k for k, c in enumerate(hd))
            f_prime_b = sum(k * c * base**(k-1) for k, c in enumerate(fd) if k > 0)

            C_prime_b = g_prime_b * h_b + g_b * h_prime_b - f_prime_b
            C_prime_vals.append(abs(C_prime_b))

        if C_prime_vals:
            arr = np.array(C_prime_vals)
            gap_arr = np.array(min_gaps)
            D_approx = 2 * bits
            pr(f"  {bits}-bit (D ~ {D_approx}):")
            pr(f"    |C'(b)|: mean={arr.mean():.2e}, min={arr.min():.2e}, max={arr.max():.2e}")
            pr(f"    min_gap (theta!=0): mean={gap_arr.mean():.2e}, min={gap_arr.min():.2e}")
            pr(f"    |C'(b)|/b^D: mean={arr.mean()/base**D_approx:.6e}")
            derivative_data.append((bits, arr.mean(), gap_arr.min()))

    # ═══════════════════════════════════════════════════════════════
    # PART B: FOURIER DECOMPOSITION OF THE GAP
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART B: FOURIER DECOMPOSITION OF |C(be^{itheta})|^2")
    pr(f"{'═' * 72}")
    pr()
    pr("|C(z)|^2 = |g(z)|^2 * |h(z)|^2 - 2*Re(g(z)h(z)*conj(f(z))) + |f(z)|^2")
    pr("Each term has a cosine series in theta when z = be^{itheta}.")
    pr()

    for bits in [10, 14]:
        pr(f"  {bits}-bit examples:")
        for trial in range(3):
            p = random_prime(bits)
            q = random_prime(bits)
            if p == q:
                q = random_prime(bits)
            N = p * q
            gd = to_digits(p, base)
            hd = to_digits(q, base)
            fd = to_digits(N, base)

            n_theta = 4096
            theta = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
            z = base * np.exp(1j * theta)

            g_vals = sum(c * z**k for k, c in enumerate(gd))
            h_vals = sum(c * z**k for k, c in enumerate(hd))
            f_vals = sum(c * z**k for k, c in enumerate(fd))

            C_vals = g_vals * h_vals - f_vals
            C_abs2 = np.abs(C_vals)**2

            fft_C2 = np.fft.rfft(C_abs2)
            power = np.abs(fft_C2) / len(fft_C2)
            dominant = np.argsort(-power)[:5]
            dom_str = ", ".join(f"mode {m}: {power[m]:.2e}" for m in dominant)
            pr(f"    p={p}, q={q}: {dom_str}")

    # ═══════════════════════════════════════════════════════════════
    # PART C: STRUCTURED PROFILE ANALYSIS
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART C: STRUCTURED CARRY PROFILES — EXACT GAP ANALYSIS")
    pr(f"{'═' * 72}")
    pr()
    pr("For a carry profile c = [c_1, ..., c_D], the quotient polynomial is")
    pr("Q(z) = c_1 + c_2*z + ... + c_D*z^{D-1}.")
    pr("C(z) = (z-b)*Q(z), so |C(be^{itheta})| = |be^{itheta}-b| * |Q(be^{itheta})|")
    pr("     = 2b*|sin(theta/2)| * |Q(be^{itheta})|.")
    pr()
    pr("For r_max <= b, we need Q(be^{itheta}) != 0 for all theta.")
    pr("Equivalently: Q has no roots on |z| = b.")
    pr()

    test_profiles = [
        ("Degenerate [0,...,0,2,1]", [0, 0, 0, 2, 1]),
        ("Constant [1,1,1,1,1]", [1, 1, 1, 1, 1]),
        ("Alternating [1,0,1,0,1]", [1, 0, 1, 0, 1]),
        ("Ramp [1,1,2,2,1]", [1, 1, 2, 2, 1]),
        ("Peaked [1,2,3,2,1]", [1, 2, 3, 2, 1]),
        ("Max carry [1,2,3,3,2,1]", [1, 2, 3, 3, 2, 1]),
        ("Long constant D=8", [1]*8),
        ("Long alternating D=8", [1,0]*4),
        ("Near-degenerate [0,0,0,0,0,2,1]", [0,0,0,0,0,2,1]),
    ]

    pr(f"  {'Profile':>35} | {'D':>3} | min|Q(be^it)| theta!=0 | {'r_max':>8} | r_max < b?")
    pr(f"  {'-'*35}-+-{'-'*3}-+-{'-'*25}-+-{'-'*8}-+-{'-'*10}")

    for name, carries in test_profiles:
        D = len(carries)
        n_theta = 8192
        theta = np.linspace(0, 2*np.pi, n_theta, endpoint=False)
        z = base * np.exp(1j * theta)

        Q_vals = sum(carries[k] * z**k for k in range(D))

        exclude = max(1, n_theta // 256)
        Q_abs = np.abs(Q_vals)
        Q_min = Q_abs[exclude:-exclude].min() if len(Q_abs) > 2*exclude else Q_abs.min()

        M = np.zeros((D, D))
        lead = carries[-1]
        if lead != 0:
            for i in range(D - 1):
                M[i + 1, i] = 1.0
            for i in range(D):
                M[i, D - 1] = -carries[i] / lead
            ev = np.linalg.eigvals(M)
            rmax = np.max(np.abs(ev))
        else:
            rmax = float('nan')

        status = "YES" if rmax < base else ("EXACT" if abs(rmax - base) < 1e-10 else "NO")
        pr(f"  {name:>35} | {D:>3} | {Q_min:25.8e} | {rmax:8.4f} | {status}")

    # ═══════════════════════════════════════════════════════════════
    # PART D: GAP GROWTH WITH BIT SIZE (STRENGTH OF ROUCHÉ)
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART D: GAP GROWTH — DOES min|C(be^{itheta})| GROW WITH D?")
    pr(f"{'═' * 72}")
    pr()
    pr("The plan paper claims gap_min ~ 2^{2K}. We test this.")
    pr()

    bit_sizes = [8, 10, 12, 14, 16, 18, 20]
    gap_data = []

    for bits in bit_sizes:
        n_samples = 3000 if bits <= 16 else 1000
        min_Q_gaps = []

        for _ in range(n_samples):
            p = random_prime(bits)
            q = random_prime(bits)
            if p == q:
                continue
            N = p * q
            gd = to_digits(p, base)
            hd = to_digits(q, base)
            fd = to_digits(N, base)

            carries_list = compute_carries(p, q, base)
            carry_seq = carries_list[1:]
            if not carry_seq or carry_seq[-1] == 0:
                continue

            D = len(carry_seq)
            n_theta = 4096
            theta = np.linspace(0, 2*np.pi, n_theta, endpoint=False)
            z = base * np.exp(1j * theta)

            Q_vals = sum(carry_seq[k] * z**k for k in range(D))
            Q_abs = np.abs(Q_vals)
            exclude = max(1, n_theta // 128)
            Q_min = Q_abs[exclude:-exclude].min()
            min_Q_gaps.append((D, Q_min))

        if min_Q_gaps:
            Ds = np.array([d for d, _ in min_Q_gaps], dtype=float)
            gaps = np.array([g for _, g in min_Q_gaps])
            D_mean = Ds.mean()

            p5 = np.percentile(gaps, 5)
            p50 = np.percentile(gaps, 50)
            min_g = gaps.min()

            gap_data.append((bits, D_mean, min_g, p5, p50))
            pr(f"  {bits:3d}-bit (D~{D_mean:.0f}): min_gap={min_g:.4e}, "
               f"p5={p5:.4e}, median={p50:.4e}")

    if len(gap_data) >= 3:
        pr()
        pr("  Growth analysis:")
        Ds = np.array([d[1] for d in gap_data])
        mins = np.array([d[2] for d in gap_data])
        log_mins = np.log(mins + 1e-30)
        log_Ds = np.log(Ds)

        A = np.column_stack([Ds, np.ones_like(Ds)])
        result = np.linalg.lstsq(A, log_mins, rcond=None)
        slope, intercept = result[0]
        pr(f"    Exponential fit: min_gap ~ exp({slope:.4f} * D + {intercept:.2f})")
        pr(f"    Base of exponential: {np.exp(slope):.4f} (expected ~2^2 = 4 if gap ~ 2^{'{'}2K{'}'})")
        pr(f"    Actually: 2^({slope/np.log(2):.3f} * D)")

        for i in range(1, len(gap_data)):
            ratio = gap_data[i][2] / gap_data[i-1][2]
            dD = gap_data[i][1] - gap_data[i-1][1]
            pr(f"    {gap_data[i-1][0]}-bit -> {gap_data[i][0]}-bit: "
               f"gap ratio = {ratio:.3f}, dD ~ {dD:.0f}")

    # ═══════════════════════════════════════════════════════════════
    # PART E: KEY STRUCTURAL INSIGHT
    # ═══════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("PART E: THE FACTORED FORM C(z) = (z-b)*Q(z)")
    pr(f"{'═' * 72}")
    pr()
    pr("On |z| = b:")
    pr("  |C(be^{itheta})| = |be^{itheta} - b| * |Q(be^{itheta})|")
    pr("                   = 2b|sin(theta/2)| * |Q(be^{itheta})|")
    pr()
    pr("r_max <= b  <=>  Q has no roots on |z| = b")
    pr("            <=>  min_{theta} |Q(be^{itheta})| > 0")
    pr()
    pr("For Q(z) = c_D*z^{D-1} + ... + c_1 with c_k >= 0 (carries):")
    pr("  |Q(be^{itheta})| >= |c_D * b^{D-1}| - sum_{k<D} c_k * b^{k-1}")
    pr("                   = b^{D-1}(c_D - sum c_k/b^{D-k})")
    pr()
    pr("This is the standard Rouche form. The condition becomes:")
    pr("  c_D > sum_{k=1}^{D-1} c_k / b^{D-k}")
    pr()
    pr("Testing this condition on real carry profiles:")

    n_test = 10000
    rouche_pass = 0
    rouche_fail = 0
    fail_examples = []

    for _ in range(n_test):
        bits = random.randint(8, 24)
        p = random_prime(bits)
        q = random_prime(bits)
        if p == q:
            continue
        carries_list = compute_carries(p, q, base)
        carry_seq = carries_list[1:]
        if not carry_seq or carry_seq[-1] == 0:
            continue

        D = len(carry_seq)
        c_D = carry_seq[-1]
        rouche_sum = sum(carry_seq[k] / base**(D - 1 - k) for k in range(D - 1))

        if c_D > rouche_sum:
            rouche_pass += 1
        else:
            rouche_fail += 1
            if len(fail_examples) < 5:
                fail_examples.append((p, q, carry_seq, c_D, rouche_sum))

    total = rouche_pass + rouche_fail
    pr(f"\n  Rouche condition c_D > sum c_k/b^{{D-k}} satisfied: "
       f"{rouche_pass}/{total} ({100*rouche_pass/total:.1f}%)")

    if fail_examples:
        pr(f"\n  Examples where Rouche fails (but r_max may still be <= b):")
        for p, q, cs, cD, rs in fail_examples:
            D = len(cs)
            lead = cs[-1]
            M = np.zeros((D, D))
            for i in range(D - 1):
                M[i + 1, i] = 1.0
            for i in range(D):
                M[i, D - 1] = -cs[i] / lead
            ev = np.linalg.eigvals(M)
            rmax = np.max(np.abs(ev))
            pr(f"    p={p}, q={q}: c_D={cD}, sum={rs:.4f}, "
               f"r_max={rmax:.6f} ({'<= b' if rmax <= base + 1e-10 else '> b!'})")

    pr()
    pr("A21 COMPLETE")


if __name__ == "__main__":
    main()
