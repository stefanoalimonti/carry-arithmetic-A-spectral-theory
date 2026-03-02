#!/usr/bin/env python3
"""
A19: Determine the correct polynomial exponent in Proposition 3.

Proposition 3 claims |lambda_k(j) - 1/b^k| = O(j^{k-3} * b^{-(j-1)}).
We verify this by computing exact eigenvalues lambda_k(j) for k=3,4,5
and j=2,...,10, then fitting |lambda_k(j) - 1/b^k| ~ C * j^alpha * (1/b)^j.

Result: alpha ≈ k-3 for all k tested, confirming the tight exponent.
"""

import numpy as np
from fractions import Fraction
from itertools import product as iproduct
import sys

def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def compute_mult_conv_dist_exact(j, base=2):
    """Exact distribution of conv_j for base-b multiplication.
    MSB digits g_0, h_0 are uniform on {1, ..., b-1} (nonzero).
    Remaining digits g_i, h_i are uniform on {0, ..., b-1}.
    Returns dict {v: Fraction(prob)}."""
    digits_g = list(range(1, base)) + [list(range(base))] * j
    digits_h = list(range(1, base)) + [list(range(base))] * j

    g_choices = [list(range(1, base))] + [list(range(base))] * j
    h_choices = [list(range(1, base))] + [list(range(base))] * j

    dist = {}
    total = 0
    for g_digits in iproduct(*g_choices):
        for h_digits in iproduct(*h_choices):
            v = sum(g_digits[i] * h_digits[j - i] for i in range(j + 1))
            dist[v] = dist.get(v, 0) + 1
            total += 1

    return {v: Fraction(cnt, total) for v, cnt in dist.items()}


def build_transfer_matrix_exact(dist, base=2):
    """Build carry transition matrix using Fraction arithmetic.
    Returns (dim, T) where T[c_in][c_out] is a Fraction."""
    v_max = max(dist.keys())
    c_max = v_max // (base - 1) + 1
    dim = c_max + 1

    T = [[Fraction(0)] * dim for _ in range(dim)]
    for c_in in range(dim):
        for v, p in dist.items():
            c_out = (v + c_in) // base
            if c_out < dim:
                T[c_in][c_out] += p
    return dim, T


def eigenvalues_from_exact_matrix(dim, T_frac):
    """Convert exact Fraction matrix to float and compute eigenvalues."""
    T_float = np.array([[float(T_frac[i][j]) for j in range(dim)]
                         for i in range(dim)])
    eigs = np.linalg.eigvals(T_float)
    return sorted(np.abs(eigs), reverse=True)


def main():
    pr("=" * 72)
    pr("  A19: PROPOSITION 3 EXPONENT DETERMINATION")
    pr("=" * 72)
    pr()
    pr("Goal: verify that |lambda_k(j) - 1/b^k| ~ j^{k-3} * (1/b)^j")
    pr()

    base = 2
    j_values = list(range(2, 13))
    k_targets = [3, 4, 5]

    results = {k: [] for k in k_targets}

    for j in j_values:
        pr(f"  j = {j}: computing exact conv distribution...", end=" ")
        if j > 10:
            pr("(skipped, too large for exact enumeration)")
            continue

        dist = compute_mult_conv_dist_exact(j, base)
        dim, T_frac = build_transfer_matrix_exact(dist, base)
        eigs = eigenvalues_from_exact_matrix(dim, T_frac)

        pr(f"dim={dim}, top eigenvalues: ", end="")
        for idx, k in enumerate(k_targets):
            if k < len(eigs):
                lam_k = eigs[k]
                target = 1.0 / base**k
                err = abs(lam_k - target)
                results[k].append((j, err))
                pr(f"λ_{k}={lam_k:.8f} (err={err:.2e})", end="  ")
        pr()

    pr()
    pr("=" * 72)
    pr("  FITTING: |λ_k(j) - 1/b^k| ~ C * j^α * (1/b)^j")
    pr("=" * 72)
    pr()

    for k in k_targets:
        data = [(j, err) for j, err in results[k] if err > 1e-15]
        if len(data) < 3:
            pr(f"  k={k}: insufficient data points ({len(data)})")
            continue

        js = np.array([d[0] for d in data], dtype=float)
        errs = np.array([d[1] for d in data], dtype=float)

        log_errs = np.log(errs)
        log_base_term = -js * np.log(base)

        residuals = log_errs - log_base_term
        log_js = np.log(js)

        A = np.column_stack([log_js, np.ones_like(log_js)])
        result = np.linalg.lstsq(A, residuals, rcond=None)
        alpha_fit, log_C = result[0]

        pr(f"  k = {k}:")
        pr(f"    Fitted α = {alpha_fit:.4f}")
        pr(f"    Expected k-3 = {k-3}, k-2 = {k-2}")
        pr(f"    Fitted C = {np.exp(log_C):.6f}")

        closer_to = "k-3" if abs(alpha_fit - (k-3)) < abs(alpha_fit - (k-2)) else "k-2"
        pr(f"    --> Closer to: {closer_to} (distance to k-3: {abs(alpha_fit-(k-3)):.3f}, to k-2: {abs(alpha_fit-(k-2)):.3f})")
        pr()

    pr()
    pr("=" * 72)
    pr("  RAW DATA TABLE")
    pr("=" * 72)
    pr()
    pr(f"{'j':>4} | ", end="")
    for k in k_targets:
        pr(f"|λ_{k} - 1/{base}^{k}|        ", end="")
    pr()
    pr("-" * 60)
    for j in j_values:
        if j > 10:
            continue
        pr(f"{j:>4} | ", end="")
        for k in k_targets:
            matches = [err for jj, err in results[k] if jj == j]
            if matches:
                pr(f"{matches[0]:20.12e}  ", end="")
            else:
                pr(f"{'N/A':>20s}  ", end="")
        pr()

    pr()
    pr("=" * 72)
    pr("  RATIO TEST: err(j+1)/err(j) should approach 1/b = 0.5")
    pr("=" * 72)
    pr()
    for k in k_targets:
        data = [(j, err) for j, err in results[k] if err > 1e-15]
        if len(data) < 2:
            continue
        pr(f"  k = {k}:")
        for i in range(1, len(data)):
            j_prev, err_prev = data[i-1]
            j_curr, err_curr = data[i]
            if err_prev > 0:
                ratio = err_curr / err_prev
                pr(f"    j={j_prev}→{j_curr}: ratio = {ratio:.6f}")
        pr()

    pr("A19 COMPLETE")


if __name__ == "__main__":
    main()
