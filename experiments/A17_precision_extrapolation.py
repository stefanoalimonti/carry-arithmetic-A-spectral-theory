"""
A17_precision_extrapolation.py — High-precision extrapolation of c₁(K)

Uses EXACT integer values of sum_cm1 and n_ulc for K=4..16 to compute
c₁(K) as exact fractions, then applies:
1. Wynn-epsilon (iterated Shanks) with mpmath arbitrary precision
2. Levin u-transform
3. PSLQ on the extrapolated value
"""

from fractions import Fraction
from mpmath import mp, mpf, pi as mpi, log as mlog, pslq, nstr, euler, catalan
import math

mp.dps = 80

exact_data = {
    4:  (41, 41),
    5:  (205, 195),
    6:  (956, 881),
    7:  (4174, 3759),
    8:  (17717, 15635),
    9:  (73301, 63853),
    10: (299404, 258613),
    11: (1212075, 1041063),
    12: (4882521, 4178907),
    13: (19607966, 16745419),
    14: (78609132, 67045167),
    15: (314829727, 268306599),
    16: (1260189799, 1073488167),
}

print("A17: HIGH-PRECISION EXTRAPOLATION")
print("=" * 70)

c1_frac = {}
c1_mp = {}
for K in sorted(exact_data):
    s, n = exact_data[K]
    c1_frac[K] = Fraction(s - n, n)
    c1_mp[K] = mpf(s) / mpf(n) - 1
    diff = c1_mp[K] - mpi/18
    print(f"K={K:2d}: c₁ = {nstr(c1_mp[K], 18):>22s}  Δ = {nstr(diff, 8)}")

print(f"\nπ/18 = {nstr(mpi/18, 20)}")

target = mpi / 18
Ks = sorted(c1_mp.keys())
vals = [c1_mp[K] for K in Ks]

print(f"\n{'='*70}")
print("PART 1: Wynn-epsilon algorithm")
print(f"{'='*70}")

def wynn_epsilon(seq):
    """Wynn epsilon algorithm for sequence extrapolation."""
    n = len(seq)
    eps = [[mpf(0)] * (n + 1) for _ in range(n + 1)]

    for i in range(n):
        eps[i][1] = seq[i]

    for k in range(2, n + 1):
        for i in range(n - k + 1):
            diff = eps[i+1][k-1] - eps[i][k-1]
            if abs(diff) < mpf('1e-70'):
                eps[i][k] = mpf('1e70')
            else:
                eps[i][k] = eps[i+1][k-2] + 1/diff

    results = []
    for k in range(1, n + 1, 2):
        row = []
        for i in range(n - k + 1):
            row.append(eps[i][k])
        results.append(row)

    return eps, results


eps, results = wynn_epsilon(vals)

print(f"\nWynn-epsilon diagonal (best estimates):")
for level in range(len(results)):
    if len(results[level]) > 0:
        best = results[level][-1]
        diff = best - target
        print(f"  Level {level}: {nstr(best, 20)}  Δ = {nstr(diff, 10)}"
              f"  (from {len(results[level])} values)")

print(f"\n{'='*70}")
print("PART 2: Richardson extrapolation chain")
print(f"{'='*70}")

def richardson_chain(seq, K_list):
    """Iterated Richardson extrapolation assuming c₁ ~ c∞ + A·r^K."""
    current = list(seq)
    current_K = list(K_list)
    level = 0

    print(f"\nLevel 0: {len(current)} values")
    for i in range(max(0, len(current)-4), len(current)):
        diff = current[i] - target
        print(f"  K={current_K[i]}: {nstr(current[i], 18)}  Δ = {nstr(diff, 8)}")

    while len(current) >= 3:
        level += 1
        new_vals = []
        new_K = []
        for i in range(len(current) - 2):
            a, b, c = current[i], current[i+1], current[i+2]
            denom = a - 2*b + c
            if abs(denom) > mpf('1e-60'):
                s = (a*c - b*b) / denom
                new_vals.append(s)
                new_K.append(current_K[i+1])

        current = new_vals
        current_K = new_K
        if len(current) == 0:
            break

        print(f"\nLevel {level}: {len(current)} values")
        for i in range(max(0, len(current)-4), len(current)):
            diff = current[i] - target
            print(f"  K~{current_K[i]}: {nstr(current[i], 18)}  Δ = {nstr(diff, 10)}")

    if len(current) > 0:
        return current[-1]
    return None


best_rich = richardson_chain(vals, Ks)

print(f"\n{'='*70}")
print("PART 3: Best estimate comparison")
print(f"{'='*70}")

wynn_best = None
for level in range(len(results)-1, -1, -1):
    if len(results[level]) > 0:
        candidate = results[level][-1]
        if abs(candidate - target) < mpf(1):
            wynn_best = candidate
            break

if wynn_best:
    print(f"\nWynn-epsilon best: {nstr(wynn_best, 20)}")
    print(f"  Δ from π/18:    {nstr(wynn_best - target, 12)}")

if best_rich:
    print(f"\nRichardson best:   {nstr(best_rich, 20)}")
    print(f"  Δ from π/18:    {nstr(best_rich - target, 12)}")

print(f"\n{'='*70}")
print("PART 4: Geometric extrapolation with variable r(K)")
print(f"{'='*70}")

deltas = {K: c1_mp[K] - target for K in Ks}
print(f"\nΔ(K) and convergence ratios:")
for i in range(1, len(Ks)):
    K = Ks[i]
    Km1 = Ks[i-1]
    r = deltas[K] / deltas[Km1] if abs(deltas[Km1]) > 1e-20 else 0
    print(f"  K={K:2d}: Δ = {nstr(deltas[K], 10):>16s}  r = {nstr(r, 8)}")

print(f"\nFit r(K) = a + b/K:")
K_arr = [mpf(K) for K in Ks[4:]]
r_arr = [deltas[Ks[i]] / deltas[Ks[i-1]] for i in range(5, len(Ks))]

if len(r_arr) >= 3:
    sum_K = sum(1/k for k in K_arr[:len(r_arr)])
    sum_r = sum(r_arr)
    sum_Kr = sum(r/k for r, k in zip(r_arr, K_arr))
    sum_K2 = sum(1/(k*k) for k in K_arr[:len(r_arr)])
    n = len(r_arr)
    det = n * sum_K2 - sum_K * sum_K
    if abs(det) > 1e-20:
        a_fit = (sum_r * sum_K2 - sum_K * sum_Kr) / det
        b_fit = (n * sum_Kr - sum_K * sum_r) / det
        print(f"  r(K) ≈ {nstr(a_fit, 8)} + {nstr(b_fit, 6)}/K")
        print(f"  r(∞) = {nstr(a_fit, 8)}")
        r_inf = a_fit

        print(f"\n  Extrapolation using r(K) model:")
        for K_last in [14, 15, 16]:
            r_eff = a_fit + b_fit / mpf(K_last)
            c1_inf = c1_mp[K_last] + deltas[K_last] * r_eff / (1 - r_eff)
            print(f"    From K={K_last}: c₁(∞) = {nstr(c1_inf, 18)}, "
                  f"Δ = {nstr(c1_inf - target, 10)}")

print(f"\n{'='*70}")
print("PART 5: PSLQ on best estimates")
print(f"{'='*70}")

estimates = []
if wynn_best and abs(wynn_best - target) < mpf('0.01'):
    estimates.append(("Wynn", wynn_best))
if best_rich and abs(best_rich - target) < mpf('0.01'):
    estimates.append(("Richardson", best_rich))
estimates.append(("π/18 (exact)", target))

for name, est in estimates:
    print(f"\nPSLQ on {name} = {nstr(est, 18)}:")
    bases = [
        ("{1, π, ln2}", [est, mpi, mlog(2)]),
        ("{1, π, ln2, π²}", [est, mpi, mlog(2), mpi**2]),
        ("{1, π, ln2, π·ln2, ln²2}", [est, mpi, mlog(2), mpi*mlog(2), mlog(2)**2]),
    ]
    for bname, basis in bases:
        result = pslq(basis)
        if result:
            print(f"  {bname}: {result}")
        else:
            print(f"  {bname}: no relation")

print(f"\n{'='*70}")
print("PART 6: Denominator structure")
print(f"{'='*70}")

print(f"\n{'K':>4s} {'n_ulc':>12s} {'excess':>12s} {'n_ulc/2^(2K-2)':>16s} {'excess/n_ulc':>16s}")
for K in Ks:
    s, n = exact_data[K]
    excess = s - n
    frac_ulc = n / (1 << (2*(K-1)))
    c1_val = excess / n
    print(f"{K:4d} {n:12d} {excess:12d} {frac_ulc:16.10f} {c1_val:16.12f}")

print(f"\nFactorization of n_ulc:")
def factorize(n):
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors

for K in Ks:
    _, n = exact_data[K]
    f = factorize(n)
    print(f"  K={K:2d}: n_ulc={n} = {'·'.join(str(x) for x in f)}")

print(f"\n{'='*70}")
print("PART 7: OEIS search preparation")
print(f"{'='*70}")

print(f"\nn_ulc sequence: {', '.join(str(exact_data[K][1]) for K in Ks)}")
print(f"\nexcess sequence: {', '.join(str(exact_data[K][0]-exact_data[K][1]) for K in Ks)}")
print(f"\nsum_cm1 sequence: {', '.join(str(exact_data[K][0]) for K in Ks)}")
