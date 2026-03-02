#!/usr/bin/env python3
"""
A20: Compute explicit constant C in Theorem 1 perturbation bound.

Theorem 1 states ||T - T_bar||_{ell^1} <= C * rho, where the paper
gives C <= (N+1) * v_max^N. This bound is extremely loose.

For multiplication, V has m=2 independent uniform-mod-b components,
so rho=0 and the bound is vacuously satisfied. The bound is meant for
general input distributions. We test with:
  (a) Biased Bernoulli sums (p != 1/2)
  (b) Truncated geometric distributions
  (c) Shifted/skewed distributions
to measure C_actual = ||T - T_bar|| / rho for rho > 0, and search for
a tighter form than (N+1) * v_max^N.
"""

import numpy as np
from fractions import Fraction
import sys

def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def make_biased_binom(n, p_num, p_den, base=2):
    """Binomial(n, p_num/p_den) distribution as exact Fractions.
    Returns dict {v: Fraction(prob)}."""
    from math import comb
    p = Fraction(p_num, p_den)
    q = 1 - p
    dist = {}
    for k in range(n + 1):
        prob = Fraction(comb(n, k)) * p**k * q**(n - k)
        if prob > 0:
            dist[k] = prob
    return dist


def make_geometric_trunc(v_max, r_num, r_den):
    """Truncated geometric on {0, ..., v_max} with ratio r_num/r_den.
    Returns dict {v: Fraction(prob)}."""
    r = Fraction(r_num, r_den)
    raw = {}
    for v in range(v_max + 1):
        raw[v] = r ** v
    total = sum(raw.values())
    return {v: p / total for v, p in raw.items()}


def make_skewed_uniform(v_max, base=2):
    """Uniform on {0, ..., v_max} with doubled weight on multiples of base.
    Returns dict {v: Fraction(prob)}."""
    weights = {}
    for v in range(v_max + 1):
        weights[v] = Fraction(2 if v % base == 0 else 1)
    total = sum(weights.values())
    return {v: w / total for v, w in weights.items()}


def build_T_and_Tbar(dist, base=2):
    """Build T and T_bar (averaged over c mod b classes).
    Returns (dim, T, T_bar, rho, v_max, N)."""
    v_max = max(dist.keys())
    N = v_max // (base - 1)
    dim = N + 1

    T = [[Fraction(0)] * dim for _ in range(dim)]
    for c_in in range(dim):
        for v, p in dist.items():
            c_out = (v + c_in) // base
            if c_out < dim:
                T[c_in][c_out] += p

    T_bar = [[Fraction(0)] * dim for _ in range(dim)]
    for c_in in range(dim):
        c0 = c_in % base
        count_in_class = sum(1 for c_ref in range(dim) if c_ref % base == c0)
        avg_row = [Fraction(0)] * dim
        for c_ref in range(dim):
            if c_ref % base == c0:
                for j_col in range(dim):
                    avg_row[j_col] += T[c_ref][j_col]
        T_bar[c_in] = [x / count_in_class for x in avg_row]

    fourier_coeffs = []
    for m in range(1, base):
        omega_re = np.cos(2 * np.pi * m / base)
        omega_im = np.sin(2 * np.pi * m / base)
        coeff_re = sum(float(p) * np.cos(2 * np.pi * m * v / base)
                       for v, p in dist.items())
        coeff_im = sum(float(p) * np.sin(2 * np.pi * m * v / base)
                       for v, p in dist.items())
        fourier_coeffs.append(np.sqrt(coeff_re**2 + coeff_im**2))
    rho = max(fourier_coeffs) if fourier_coeffs else 0.0

    return dim, T, T_bar, rho, v_max, N


def row_norm(A, B, dim):
    """max_row sum_col |A-B|."""
    max_rs = Fraction(0)
    for i in range(dim):
        rs = Fraction(0)
        for j in range(dim):
            rs += abs(A[i][j] - B[i][j])
        if rs > max_rs:
            max_rs = rs
    return max_rs


def main():
    pr("=" * 72)
    pr("  A20: THEOREM 1 CONSTANT C — EXPLICIT COMPUTATION")
    pr("=" * 72)
    pr()
    pr("Theorem 1: ||T - T_bar||_{row} <= C * rho")
    pr("Paper claims: C <= (N+1) * v_max^N")
    pr()

    test_cases = []

    for n in [3, 4, 5, 6, 8, 10]:
        for p_num, p_den in [(1, 3), (1, 4), (2, 5), (3, 4)]:
            test_cases.append((f"Binom({n},{p_num}/{p_den})", 2,
                               lambda n=n, pn=p_num, pd=p_den: make_biased_binom(n, pn, pd, 2)))

    for v_max in [4, 6, 8, 10]:
        for r_num, r_den in [(1, 3), (2, 3), (3, 4)]:
            test_cases.append((f"Geom({v_max},{r_num}/{r_den})", 2,
                               lambda vm=v_max, rn=r_num, rd=r_den: make_geometric_trunc(vm, rn, rd)))

    for v_max in [4, 6, 8]:
        test_cases.append((f"Skewed({v_max})", 2,
                           lambda vm=v_max: make_skewed_uniform(vm, 2)))

    for n in [4, 6, 8]:
        for p_num, p_den in [(1, 4), (1, 2)]:
            test_cases.append((f"Binom({n},{p_num}/{p_den},b=3)", 3,
                               lambda n=n, pn=p_num, pd=p_den: make_biased_binom(n, pn, pd, 3)))

    pr(f"{'Name':>28} | {'b':>2} | {'dim':>4} | {'v_max':>5} | {'N':>3} | "
       f"{'rho':>10} | {'||T-Tbar||':>12} | {'C_actual':>12} | {'C_paper':>14} | {'slack':>10}")
    pr("-" * 130)

    all_results = []

    for name, base, dist_fn in test_cases:
        try:
            dist = dist_fn()
            dim, T, T_bar, rho, v_max, N = build_T_and_Tbar(dist, base)

            if rho < 1e-14:
                continue

            norm_val = row_norm(T, T_bar, dim)
            C_actual = float(norm_val) / rho
            C_paper = (N + 1) * float(v_max) ** N
            slack = C_paper / C_actual if C_actual > 0 else float('inf')

            pr(f"{name:>28} | {base:>2} | {dim:>4} | {v_max:>5} | {N:>3} | "
               f"{rho:>10.6f} | {float(norm_val):>12.8f} | {C_actual:>12.4f} | {C_paper:>14.1f} | {slack:>10.1f}x")

            all_results.append({
                'name': name, 'base': base, 'dim': dim, 'v_max': v_max,
                'N': N, 'rho': rho, 'norm': float(norm_val),
                'C_actual': C_actual, 'C_paper': C_paper, 'slack': slack
            })
        except Exception as e:
            pr(f"{name:>28} | ERROR: {e}")

    pr()
    pr("=" * 72)
    pr("  ANALYSIS")
    pr("=" * 72)
    pr()

    if all_results:
        Ns = np.array([r['N'] for r in all_results], dtype=float)
        Cs = np.array([r['C_actual'] for r in all_results])
        vmaxs = np.array([r['v_max'] for r in all_results], dtype=float)
        bases = np.array([r['base'] for r in all_results], dtype=float)
        slacks = np.array([r['slack'] for r in all_results])

        pr(f"  Slack factor (C_paper / C_actual):")
        pr(f"    min: {slacks.min():.1f}x")
        pr(f"    max: {slacks.max():.1f}x")
        pr(f"    median: {np.median(slacks):.1f}x")
        pr()

        b2 = [r for r in all_results if r['base'] == 2]
        if len(b2) > 3:
            pr("  Base 2 tighter bound search:")
            Ns_2 = np.array([r['N'] for r in b2], dtype=float)
            Cs_2 = np.array([r['C_actual'] for r in b2])

            candidates = {
                'N+1': Ns_2 + 1,
                '(N+1)^2': (Ns_2 + 1)**2,
                '(N+1)^3': (Ns_2 + 1)**3,
                'N*v_max': Ns_2 * np.array([r['v_max'] for r in b2], dtype=float),
                '(N+1)*v_max': (Ns_2 + 1) * np.array([r['v_max'] for r in b2], dtype=float),
            }
            for cname, cvals in candidates.items():
                valid = all(cvals[i] >= Cs_2[i] for i in range(len(b2)))
                max_ratio = max(Cs_2[i] / cvals[i] for i in range(len(b2)))
                min_ratio = min(Cs_2[i] / cvals[i] for i in range(len(b2)))
                pr(f"    {cname:>20}: valid={valid}, C/bound in [{min_ratio:.4f}, {max_ratio:.4f}]")
            pr()

        pr("  Per-row analysis (largest ||T-Tbar|| contributors):")
        for r in all_results[:5]:
            dist = None
            for name, base, dist_fn in test_cases:
                if name == r['name'] and base == r['base']:
                    dist = dist_fn()
                    break
            if dist is None:
                continue
            dim_c, T_c, Tb_c, _, _, _ = build_T_and_Tbar(dist, r['base'])
            max_row = -1
            max_row_val = Fraction(0)
            for i in range(dim_c):
                rs = sum(abs(T_c[i][j] - Tb_c[i][j]) for j in range(dim_c))
                if rs > max_row_val:
                    max_row_val = rs
                    max_row = i
            pr(f"    {r['name']}: max perturbation at row c={max_row} "
               f"(c mod {r['base']} = {max_row % r['base']}), "
               f"||row|| = {float(max_row_val):.8f}")

    pr()
    pr("A20 COMPLETE")


if __name__ == "__main__":
    main()
