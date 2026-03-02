#!/usr/bin/env python3
"""
A11: Analytical approach to Conjecture 8 (r_max <= b) via Rouché.

STRATEGY: On |z| = b, for z ≠ b, we have g(z)h(z) ≠ n(z) because:
  |g(z)| < g(b) and |h(z)| < h(b) and |n(z)| < n(b)   [strict for z ≠ b]
But this alone doesn't prove g(z)h(z) ≠ n(z).

The Rouché approach: if on |z| = R we can show |g(z)h(z)| > |n(z)|
(or vice versa), then g(z)h(z) - n(z) has the same number of roots
inside |z| < R as g(z)h(z).

This experiment:
  1. Compute min_{θ} |g(2e^{iθ})h(2e^{iθ}) - n(2e^{iθ})| for many semiprimes
  2. If always > 0 for θ ≠ 0 → Rouché approach on |z| = 2 is viable
  3. Test Fujiwara bound: 2 max(|q_{n-2}/q_{n-1}|, |q_{n-3}/q_{n-1}|^{1/2}, ...)
  4. Near-extremal analysis: identify carry profiles giving r_max close to b
  5. Test at multiple bit sizes and multiple bases
"""

import sys, os, random
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from carry_utils import random_prime, to_digits

random.seed(42)
np.random.seed(42)

BASE = 2


def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def carry_analysis(p, q, base=2):
    """Return carry polynomial Q, eigenvalues, and r_max for p*q."""
    N = p * q
    gd = to_digits(p, base)
    hd = to_digits(q, base)
    fd = to_digits(N, base)

    conv = [0] * (len(gd) + len(hd) - 1)
    for i, a in enumerate(gd):
        for j, b_val in enumerate(hd):
            conv[i + j] += a * b_val

    D_max = max(len(conv), len(fd))
    carries = [0] * (D_max + 2)
    for k in range(D_max):
        conv_k = conv[k] if k < len(conv) else 0
        carries[k + 1] = (conv_k + carries[k]) // base

    D_carry = 0
    for j in range(len(carries) - 1, 0, -1):
        if carries[j] != 0:
            D_carry = j
            break

    carry_seq = carries[1:D_carry + 1]
    if not carry_seq or carry_seq[-1] == 0:
        return None

    gd_arr = np.array(gd, dtype=float)
    hd_arr = np.array(hd, dtype=float)
    fd_arr = np.array(fd, dtype=float)

    D = len(carry_seq)
    lead = carry_seq[-1]
    M = np.zeros((D, D), dtype=float)
    for i in range(D - 1):
        M[i + 1, i] = 1.0
    for i in range(D):
        M[i, D - 1] = -carry_seq[i] / lead

    try:
        ev = np.linalg.eigvals(M)
        rmax = np.max(np.abs(ev))
    except Exception:
        return None

    return {
        'p': p, 'q': q, 'N': N,
        'gd': gd, 'hd': hd, 'fd': fd,
        'carries': carry_seq, 'eigenvalues': ev,
        'rmax': rmax, 'D': D, 'c_top': carry_seq[-1],
        'c_topm1': carry_seq[-2] if len(carry_seq) >= 2 else 0,
    }


def min_separation_on_circle(gd, hd, fd, R, n_theta=4096):
    """Compute min_{θ≠0} |g(Re^{iθ})h(Re^{iθ}) - n(Re^{iθ})| on |z|=R."""
    theta = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    z = R * np.exp(1j * theta)

    g_vals = np.zeros(n_theta, dtype=complex)
    for k, gk in enumerate(gd):
        g_vals += gk * z**k

    h_vals = np.zeros(n_theta, dtype=complex)
    for k, hk in enumerate(hd):
        h_vals += hk * z**k

    n_vals = np.zeros(n_theta, dtype=complex)
    for k, nk in enumerate(fd):
        n_vals += nk * z**k

    diff = g_vals * h_vals - n_vals
    diff_mod = np.abs(diff)
    diff_mod[0] = np.inf  # exclude θ = 0 (z = b)
    return np.min(diff_mod), theta[np.argmin(diff_mod)]


def fujiwara_bound(coeffs):
    """Fujiwara bound for roots of monic polynomial.
    coeffs = [a_0, a_1, ..., a_{n-1}] for x^n + a_{n-1}x^{n-1} + ... + a_0."""
    n = len(coeffs)
    if n == 0:
        return 0
    terms = []
    for k in range(n):
        ak = abs(coeffs[k])
        if ak > 0:
            terms.append(2 * ak ** (1.0 / (n - k)))
    terms.append(abs(coeffs[0]) ** (1.0 / n))
    return max(terms) if terms else 0


# =============================================================
# PART A: r_max distribution and near-extremal analysis
# =============================================================
pr("=" * 70)
pr("PART A: r_max distribution across bit sizes")
pr("=" * 70)

for bits in [8, 12, 16, 20, 24]:
    results = []
    n_samples = 2000

    for _ in range(n_samples):
        p = random_prime(bits)
        q = random_prime(bits)
        if p == q:
            continue
        res = carry_analysis(p, q, BASE)
        if res is not None:
            results.append(res)

    if not results:
        continue

    rmax_vals = np.array([r['rmax'] for r in results])
    rmax_ratio = rmax_vals / BASE

    pr(f"\n--- {bits}-bit primes (N = {len(results)} semiprimes) ---")
    pr(f"  r_max / b: mean = {np.mean(rmax_ratio):.4f}, "
       f"max = {np.max(rmax_ratio):.6f}, "
       f"median = {np.median(rmax_ratio):.4f}")
    pr(f"  r_max range: [{np.min(rmax_vals):.4f}, {np.max(rmax_vals):.6f}]")
    pr(f"  Fraction r_max > b-0.01: {np.mean(rmax_vals > BASE - 0.01):.6f}")
    pr(f"  Fraction r_max > b-0.001: {np.mean(rmax_vals > BASE - 0.001):.6f}")

    top5 = sorted(results, key=lambda r: r['rmax'], reverse=True)[:5]
    pr(f"\n  Top 5 near-extremal cases:")
    for t in top5:
        c = t['carries']
        pr(f"    r_max = {t['rmax']:.6f} (r/b = {t['rmax']/BASE:.6f}), "
           f"D = {t['D']}, c_top = {t['c_top']}, c_top-1 = {t['c_topm1']}, "
           f"max_carry = {max(c)}")

# =============================================================
# PART B: Rouché test on |z| = b
# =============================================================
pr("\n" + "=" * 70)
pr("PART B: Rouché test — min |g(z)h(z) - n(z)| on |z| = b")
pr("=" * 70)

for bits in [8, 12, 16, 20]:
    min_seps = []
    n_samples = 500

    for _ in range(n_samples):
        p = random_prime(bits)
        q = random_prime(bits)
        if p == q:
            continue
        res = carry_analysis(p, q, BASE)
        if res is None:
            continue
        min_sep, min_theta = min_separation_on_circle(
            res['gd'], res['hd'], res['fd'], BASE
        )
        min_seps.append((min_sep, min_theta, res))

    if not min_seps:
        continue

    seps = np.array([s[0] for s in min_seps])
    pr(f"\n--- {bits}-bit primes (N = {len(min_seps)} tested) ---")
    pr(f"  min |g·h - n| on |z|=b: min = {np.min(seps):.6f}, "
       f"median = {np.median(seps):.2f}, mean = {np.mean(seps):.2f}")
    pr(f"  ALL > 0: {np.all(seps > 0)}")

    worst = sorted(min_seps, key=lambda x: x[0])[:3]
    for sep, th, res in worst:
        pr(f"  Closest approach: |diff| = {sep:.6f} at θ = {th:.4f} "
           f"(r_max = {res['rmax']:.4f}, D = {res['D']})")

# =============================================================
# PART C: Fujiwara bound test
# =============================================================
pr("\n" + "=" * 70)
pr("PART C: Fujiwara bound comparison")
pr("=" * 70)

for bits in [12, 16, 20, 24]:
    fuji_bounds = []
    actual_rmax = []
    n_samples = 1000

    for _ in range(n_samples):
        p = random_prime(bits)
        q = random_prime(bits)
        if p == q:
            continue
        res = carry_analysis(p, q, BASE)
        if res is None:
            continue

        c = res['carries']
        lead = c[-1]
        monic_coeffs = [-c[k] / lead for k in range(len(c) - 1)]
        fb = fujiwara_bound(monic_coeffs)

        fuji_bounds.append(fb)
        actual_rmax.append(res['rmax'])

    fuji_bounds = np.array(fuji_bounds)
    actual_rmax = np.array(actual_rmax)

    pr(f"\n--- {bits}-bit primes (N = {len(fuji_bounds)} tested) ---")
    pr(f"  Fujiwara bound: mean = {np.mean(fuji_bounds):.2f}, "
       f"max = {np.max(fuji_bounds):.2f}, "
       f"median = {np.median(fuji_bounds):.2f}")
    pr(f"  Actual r_max:   mean = {np.mean(actual_rmax):.4f}, "
       f"max = {np.max(actual_rmax):.6f}")
    pr(f"  Ratio Fujiwara/actual: mean = {np.mean(fuji_bounds/actual_rmax):.2f}")
    pr(f"  Fujiwara <= b: {np.mean(fuji_bounds <= BASE):.4f} "
       f"(fraction that Fujiwara alone proves)")

# =============================================================
# PART D: Multi-base test
# =============================================================
pr("\n" + "=" * 70)
pr("PART D: Multi-base r_max / b")
pr("=" * 70)

for base in [2, 3, 5, 7, 10]:
    results = []
    bits = 16
    for _ in range(1000):
        p = random_prime(bits)
        q = random_prime(bits)
        if p == q:
            continue
        res = carry_analysis(p, q, base)
        if res is not None:
            results.append(res)

    if not results:
        continue

    rmax_vals = np.array([r['rmax'] for r in results])
    pr(f"  Base {base:2d}: max(r_max/b) = {np.max(rmax_vals)/base:.6f}, "
       f"mean(r_max/b) = {np.mean(rmax_vals)/base:.4f}, "
       f"N = {len(results)}")

# =============================================================
# PART E: Structure of near-extremal carry profiles
# =============================================================
pr("\n" + "=" * 70)
pr("PART E: Near-extremal carry profile structure (24-bit)")
pr("=" * 70)

results = []
for _ in range(5000):
    p = random_prime(24)
    q = random_prime(24)
    if p == q:
        continue
    res = carry_analysis(p, q, BASE)
    if res is not None:
        results.append(res)

top20 = sorted(results, key=lambda r: r['rmax'], reverse=True)[:20]
pr(f"\nTop 20 r_max values out of {len(results)} semiprimes:")
for i, t in enumerate(top20):
    c = t['carries']
    c_top5 = c[-5:] if len(c) >= 5 else c
    pr(f"  {i+1:2d}. r_max = {t['rmax']:.6f} (ratio {t['rmax']/BASE:.6f}), "
       f"D = {t['D']}, top carries = {c_top5}")

    if i < 3:
        ev = t['eigenvalues']
        extremal = ev[np.argmax(np.abs(ev))]
        pr(f"      Extremal eigenvalue: {extremal:.6f} "
           f"(angle = {np.angle(extremal)*180/np.pi:.1f}°)")

pr("\n" + "=" * 70)
pr("SUMMARY")
pr("=" * 70)
pr("Key questions for a proof of Conjecture 8:")
pr("1. Is min_{θ≠0} |g(be^{iθ})h(be^{iθ}) - n(be^{iθ})| > 0 ALWAYS?")
pr("   If yes → Rouché on |z|=b proves all roots inside |z|<b.")
pr("   The remaining root at z=b is divided out in Q(x)=C(x)/(x-b).")
pr("2. What structural property of g,h,n (binary digit polynomials)")
pr("   prevents g(z)h(z) = n(z) on |z|=b, z≠b?")
pr("3. Does the Fujiwara bound give r_max ≤ b for a positive fraction?")
