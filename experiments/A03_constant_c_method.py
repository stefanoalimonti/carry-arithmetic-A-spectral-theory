#!/usr/bin/env python3
"""
A03B: CONSTANT c — FIXED METHODOLOGY

The correct measurement: AGGREGATE hits and expected across all semiprimes
for each test prime l, THEN take the ratio. This avoids the noise from
individual 0/1 hits.

  rho(l) = (total_hits over all N) / (total_expected over all N)
"""

import sys, os, time, random, math
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from carry_utils import (random_prime, carry_poly_int, quotient_poly_int,
                         poly_roots_mod, primes_up_to)

random.seed(2024)
np.random.seed(2024)


def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def main():
    t0 = time.time()
    pr("=" * 72)
    pr("A03B: CONSTANT c — AGGREGATED MEASUREMENT")
    pr("=" * 72)

    configs = [
        (16, 3000),
        (24, 2000),
        (32, 1000),
    ]

    for BIT_SIZE, N_SEMI in configs:
        pr(f"\n{'═' * 72}")
        pr(f"  {BIT_SIZE}-BIT SEMIPRIMES ({N_SEMI} samples)")
        pr(f"{'═' * 72}")

        semiprimes = []
        while len(semiprimes) < N_SEMI:
            p = random_prime(BIT_SIZE)
            q = random_prime(BIT_SIZE)
            if p != q:
                semiprimes.append((p, q))

        pr(f"  Precomputing Q polynomials...")
        Q_data = []
        for p, q in semiprimes:
            C = carry_poly_int(p, q, 2)
            Q = quotient_poly_int(C, 2)
            if len(Q) >= 2:
                Q_data.append((p, q, Q))
        pr(f"  {len(Q_data)} valid (deg ≈ {np.mean([len(q)-1 for _,_,q in Q_data]):.0f})")

        test_primes = [l for l in primes_up_to(200) if l > 2]

        rho_by_l = {}
        for l in test_primes:
            total_hits = 0
            total_expected = 0.0
            n_with_roots = 0

            for p, q, Q in Q_data:
                roots = poly_roots_mod(Q, l)
                if not roots:
                    continue
                n_with_roots += 1
                n_roots = len(roots)
                pm = p % l
                qm = q % l

                hits = int(pm in roots)
                if qm != pm:
                    hits += int(qm in roots)
                    n_distinct = 2
                else:
                    n_distinct = 1

                total_hits += hits
                total_expected += n_roots / l * n_distinct

            if total_expected > 10:
                rho = total_hits / total_expected
                se = rho * math.sqrt(1.0 / total_hits + 1.0 / total_expected) if total_hits > 0 else 0.1
                rho_by_l[l] = {
                    'rho': rho,
                    'se': se,
                    'hits': total_hits,
                    'expected': total_expected,
                    'n_with_roots': n_with_roots,
                }

        target = 0.5

        pr(f"\n  {'l':>5s}  {'rho':>8s}  {'hits':>6s}  {'expected':>8s}  "
           f"{'excess':>8s}  {'c_est':>8s}")
        l_arr = []
        rho_arr = []
        se_arr = []
        for l in sorted(rho_by_l.keys()):
            r = rho_by_l[l]
            excess = r['rho'] - target
            c_est = l * excess
            l_arr.append(l)
            rho_arr.append(r['rho'])
            se_arr.append(r['se'])
            if l in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                     47, 53, 59, 67, 71, 79, 89, 97, 101, 127, 151, 191]:
                pr(f"  {l:5d}  {r['rho']:8.4f}  {r['hits']:6d}  "
                   f"{r['expected']:8.1f}  {excess:+8.4f}  {c_est:8.3f}")

        l_arr = np.array(l_arr, dtype=float)
        rho_arr = np.array(rho_arr)
        se_arr = np.array(se_arr)
        excess_arr = rho_arr - target

        # Fit: excess = c/l + d/l^2
        mask = l_arr >= 5
        l_fit = l_arr[mask]
        ex_fit = excess_arr[mask]

        A = np.column_stack([1.0 / l_fit, 1.0 / l_fit ** 2])
        try:
            coeffs = np.linalg.lstsq(A, ex_fit, rcond=None)[0]
            c_fit = coeffs[0]
            d_fit = coeffs[1]
            pred = A @ coeffs
            resid = ex_fit - pred
            rmse = np.sqrt(np.mean(resid ** 2))
            pr(f"\n  Fit rho = 0.5 + c/l + d/l² (l ≥ 5):")
            pr(f"    c = {c_fit:.4f}")
            pr(f"    d = {d_fit:.2f}")
            pr(f"    RMSE = {rmse:.5f}")
        except Exception:
            c_fit = 0.0

        # Simple c estimate from medium primes
        mask_med = (l_arr >= 11) & (l_arr <= 67)
        if np.sum(mask_med) > 3:
            c_estimates = l_arr[mask_med] * excess_arr[mask_med]
            c_simple = np.mean(c_estimates)
            c_se = np.std(c_estimates) / np.sqrt(len(c_estimates))
            pr(f"\n  Simple estimate (l=11..67): c = {c_simple:.4f} ± {c_se:.4f}")

            dist_quarter = abs(c_simple - 0.25) / max(c_se, 0.001)
            dist_pi = abs(c_simple - math.pi ** 2 / 36) / max(c_se, 0.001)
            pr(f"    Distance from 1/4:   {dist_quarter:.1f}σ")
            pr(f"    Distance from π²/36: {dist_pi:.1f}σ")

        # Large l: convergence to 0.5
        mask_large = l_arr >= 50
        if np.sum(mask_large) > 5:
            large_excess = excess_arr[mask_large]
            pr(f"\n  Large l (≥50): mean excess = {np.mean(large_excess):+.5f}")

    # ════════════════════════════════════════════════════════════════
    # MULTI-BASE MEASUREMENT
    # ════════════════════════════════════════════════════════════════
    pr(f"\n{'═' * 72}")
    pr("MULTI-BASE c(b)")
    pr(f"{'═' * 72}")

    semi_20 = []
    while len(semi_20) < 1500:
        p = random_prime(20)
        q = random_prime(20)
        if p != q:
            semi_20.append((p, q))

    for base in [2, 3, 5, 7, 10]:
        target = (base - 1.0) / base
        test_l = [l for l in primes_up_to(67) if l > base]

        c_vals = []
        for l in test_l:
            total_hits = 0
            total_exp = 0.0
            for p, q in semi_20:
                C = carry_poly_int(p, q, base)
                Q = quotient_poly_int(C, base)
                if len(Q) < 2:
                    continue
                roots = poly_roots_mod(Q, l)
                if not roots:
                    continue
                pm = p % l
                qm = q % l
                hits = int(pm in roots)
                if qm != pm:
                    hits += int(qm in roots)
                    nf = 2
                else:
                    nf = 1
                total_hits += hits
                total_exp += len(roots) / l * nf

            if total_exp > 5:
                rho = total_hits / total_exp
                c_vals.append(l * (rho - target))

        if c_vals:
            c_mean = np.mean(c_vals)
            c_se = np.std(c_vals) / np.sqrt(len(c_vals))
            pr(f"  base {base:2d}: (b-1)/b={target:.4f}  "
               f"c={c_mean:.4f}±{c_se:.4f}  "
               f"c/(b-1)={c_mean/(base-1):.4f}")

    pr(f"\n{'═' * 72}")
    pr("SYNTHESIS")
    pr(f"{'═' * 72}")
    pr(f"""
  The constant c in rho(l) = (b-1)/b + c/l determines whether
  the carry barrier has a connection to zeta(2) or is purely
  combinatorial.
  
  pi²/36 = {math.pi**2/36:.6f} (zeta connection)
  1/4    = {0.25:.6f} (combinatorial)
""")

    pr(f"\nTotal runtime: {time.time() - t0:.1f}s")
    pr("=" * 72)


if __name__ == "__main__":
    main()
