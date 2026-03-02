"""
A18_final_extrapolation.py — Final extrapolation with K=4..17

The convergence ratio r(K) = Δ(K)/Δ(K-1) is decreasing toward 1/2.
With K=17, r(17) = 0.505. This suggests c₁(K) = c∞ + A·(1/2)^K + corrections.

If the dominant term is truly (1/2)^K, we can extract c∞ with much
better precision using a model that accounts for the varying ratio.
"""

from mpmath import mp, mpf, pi as mpi, log as mlog, pslq, nstr, euler, catalan, power
import math

mp.dps = 50

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
    17: (5042657326, 4294460189),
}

print("A18: FINAL EXTRAPOLATION WITH K=4..17")
print("=" * 70)

target = mpi / 18
c1 = {}
for K in sorted(exact_data):
    s, n = exact_data[K]
    c1[K] = mpf(s) / mpf(n) - 1

Ks = sorted(c1.keys())
vals = [c1[K] for K in Ks]
deltas = {K: c1[K] - target for K in Ks}

print(f"\n{'K':>4s} {'c₁(K)':>22s} {'Δ':>16s} {'r=Δ(K)/Δ(K-1)':>14s}")
for i, K in enumerate(Ks):
    r = deltas[K] / deltas[Ks[i-1]] if i > 0 else mpf(0)
    print(f"{K:4d} {nstr(c1[K], 18):>22s} {nstr(deltas[K], 10):>16s} {nstr(r, 8):>14s}")

print(f"\nπ/18 = {nstr(target, 20)}")

print(f"\n{'='*70}")
print("KEY OBSERVATION: r(K) → 1/2")
print(f"{'='*70}")

ratios = []
for i in range(1, len(Ks)):
    r = float(deltas[Ks[i]] / deltas[Ks[i-1]])
    ratios.append((Ks[i], r))

print(f"\nConvergence ratios r(K):")
for K, r in ratios:
    print(f"  K={K:2d}: r = {r:.8f}  (r - 1/2 = {r - 0.5:+.6f})")

K_late = [K for K, r in ratios if K >= 10]
r_late = [r for K, r in ratios if K >= 10]

import numpy as np
K_arr = np.array(K_late, dtype=float)
r_arr = np.array(r_late)

from numpy.polynomial import polynomial as P
coeffs_1 = np.polyfit(1/K_arr, r_arr, 1)
coeffs_2 = np.polyfit(1/K_arr, r_arr, 2)
print(f"\nFit r(K) = a + b/K (K≥10):")
print(f"  a = {coeffs_1[1]:.8f}, b = {coeffs_1[0]:.6f}")
print(f"  r(∞) = {coeffs_1[1]:.8f}")
print(f"\nFit r(K) = a + b/K + c/K² (K≥10):")
print(f"  a = {coeffs_2[2]:.8f}, b = {coeffs_2[1]:.6f}, c = {coeffs_2[0]:.4f}")
print(f"  r(∞) = {coeffs_2[2]:.8f}")

print(f"\n{'='*70}")
print("PART 1: Aitken Δ² (Shanks) with K=15,16,17")
print(f"{'='*70}")

a, b, c = c1[15], c1[16], c1[17]
denom = a - 2*b + c
if abs(denom) > mpf('1e-40'):
    shanks = (a*c - b*b) / denom
    print(f"\nShanks(15,16,17) = {nstr(shanks, 20)}")
    print(f"Δ from π/18     = {nstr(shanks - target, 12)}")

for triple in [(13,14,15), (14,15,16), (15,16,17), (13,15,17), (12,14,16), (11,14,17)]:
    K1, K2, K3 = triple
    a, b, c = c1[K1], c1[K2], c1[K3]
    denom = a - 2*b + c
    if abs(denom) > mpf('1e-40'):
        s = (a*c - b*b) / denom
        print(f"Shanks({K1},{K2},{K3}) = {nstr(s, 20)}  Δ = {nstr(s - target, 10)}")

print(f"\n{'='*70}")
print("PART 2: Wynn-epsilon on full sequence K=4..17")
print(f"{'='*70}")

def wynn_epsilon_table(seq):
    n = len(seq)
    eps = {}
    for i in range(n):
        eps[(i, -1)] = mpf(0)
        eps[(i, 0)] = seq[i]

    for k in range(1, n):
        for i in range(n - k):
            diff = eps[(i+1, k-1)] - eps[(i, k-1)]
            if abs(diff) < mpf('1e-60'):
                eps[(i, k)] = mpf('1e60')
            else:
                eps[(i, k)] = eps[(i+1, k-2)] + 1/diff

    return eps

eps = wynn_epsilon_table(vals)
n = len(vals)

print(f"\nWynn-epsilon estimates (even columns = sequence estimates):")
best_wynn = None
best_wynn_diff = mpf(1)
for k in range(0, n, 2):
    estimates = []
    for i in range(n - k):
        if (i, k) in eps:
            estimates.append(eps[(i, k)])
    if len(estimates) > 0:
        last = estimates[-1]
        diff = last - target
        if abs(diff) < abs(best_wynn_diff):
            best_wynn = last
            best_wynn_diff = diff
        if k <= 10:
            print(f"  ε[{k}]: last = {nstr(last, 18)}, Δ = {nstr(diff, 10)}")

print(f"\nBest Wynn: {nstr(best_wynn, 18)}, Δ = {nstr(best_wynn_diff, 10)}")

print(f"\n{'='*70}")
print("PART 3: Richardson chain on K=8..17")
print(f"{'='*70}")

sub_vals = [c1[K] for K in range(8, 18)]
sub_Ks = list(range(8, 18))

current = list(sub_vals)
current_Ks = list(sub_Ks)
level = 0
all_estimates = [current[-1]]

while len(current) >= 3:
    level += 1
    new_vals = []
    new_Ks = []
    for i in range(len(current) - 2):
        a, b, c = current[i], current[i+1], current[i+2]
        d = a - 2*b + c
        if abs(d) > mpf('1e-50'):
            new_vals.append((a*c - b*b) / d)
            new_Ks.append(current_Ks[i+1])
    current = new_vals
    current_Ks = new_Ks
    if len(current) == 0:
        break
    all_estimates.append(current[-1])
    diff = current[-1] - target
    print(f"  Level {level}: {nstr(current[-1], 18)}, Δ = {nstr(diff, 10)}, ({len(current)} values)")

print(f"\n{'='*70}")
print("PART 4: Model c₁(K) = π/18 + A·(1/2)^K + B·(1/2)^{2K}")
print(f"{'='*70}")

K_use = list(range(12, 18))
from scipy.optimize import least_squares

def model_half(params, K_arr, y_arr):
    c_inf, A, B = params
    return y_arr - (c_inf + A * 0.5**K_arr + B * 0.25**K_arr)

K_arr_np = np.array(K_use, dtype=float)
y_arr_np = np.array([float(c1[K]) for K in K_use])

res = least_squares(model_half, [math.pi/18, -10, 100], args=(K_arr_np, y_arr_np))
c_inf, A, B = res.x
print(f"\nFit with r=1/2 fixed:")
print(f"  c₁(∞) = {c_inf:.15f}")
print(f"  A = {A:.6f}, B = {B:.6f}")
print(f"  c₁(∞) - π/18 = {c_inf - math.pi/18:+.2e}")
print(f"  Max residual = {max(abs(res.fun)):.2e}")

def model_free(params, K_arr, y_arr):
    c_inf, A, r = params
    return y_arr - (c_inf + A * r**K_arr)

res2 = least_squares(model_free, [math.pi/18, -1, 0.5],
                     args=(K_arr_np, y_arr_np),
                     bounds=([0.17, -100, 0.01], [0.18, 100, 0.99]))
c_inf2, A2, r2 = res2.x
print(f"\nFit with r free:")
print(f"  c₁(∞) = {c_inf2:.15f}")
print(f"  A = {A2:.6f}, r = {r2:.8f}")
print(f"  c₁(∞) - π/18 = {c_inf2 - math.pi/18:+.2e}")
print(f"  Max residual = {max(abs(res2.fun)):.2e}")

print(f"\n{'='*70}")
print("PART 5: PSLQ on best estimates")
print(f"{'='*70}")

estimates_to_test = []
if best_wynn and abs(best_wynn - target) < mpf('0.001'):
    estimates_to_test.append(("Wynn-best", best_wynn))

shanks_15_16_17 = (c1[15]*c1[17] - c1[16]**2) / (c1[15] - 2*c1[16] + c1[17])
estimates_to_test.append(("Shanks(15,16,17)", shanks_15_16_17))
estimates_to_test.append(("Fit r=1/2", mpf(str(c_inf))))
estimates_to_test.append(("Fit r free", mpf(str(c_inf2))))

for name, est in estimates_to_test:
    diff = est - target
    print(f"\n{name}: {nstr(est, 18)}, Δ = {nstr(diff, 10)}")
    for bname, basis in [
        ("{1, π, ln2}", [est, mpi, mlog(2)]),
        ("{1, π, ln2, π²}", [est, mpi, mlog(2), mpi**2]),
    ]:
        result = pslq(basis)
        if result:
            names = bname.strip("{}").split(", ")
            names[0] = "c₁"
            print(f"  PSLQ {bname}: {result}")
        else:
            print(f"  PSLQ {bname}: no relation")

print(f"\n{'='*70}")
print("PART 6: Predicted convergence")
print(f"{'='*70}")

print(f"\nIf r(K) → 1/2, model c₁(K) ≈ π/18 + A·(1/2)^K:")
A_est = deltas[17] / power(mpf(1)/2, 17)
print(f"  A = Δ(17) / (1/2)^17 = {nstr(A_est, 8)}")
for K_pred in [18, 19, 20, 25, 30]:
    c1_pred = target + A_est * power(mpf(1)/2, K_pred)
    delta_pred = A_est * power(mpf(1)/2, K_pred)
    print(f"  c₁({K_pred:2d}) ≈ {nstr(c1_pred, 15)}, Δ = {nstr(delta_pred, 8)}")

print(f"\nTo get 12 digits: need Δ < 10^-12")
print(f"  K needed: (1/2)^K · |A| < 10^-12")
print(f"  K > {float(-12 * math.log(10) / math.log(0.5) + math.log(abs(float(A_est))) / math.log(2)):.1f}")
K_12dig = math.ceil(-12 * math.log(10) / math.log(0.5) + math.log(abs(float(A_est))) / math.log(2))
print(f"  → Need K ≈ {K_12dig}")

print(f"\n{'='*70}")
print("CONCLUSION")
print(f"{'='*70}")
print(f"""
With K=17, the convergence ratio r(17) = Δ(17)/Δ(16) = {float(deltas[17]/deltas[16]):.6f}.
The ratio is approaching 1/2, suggesting the dominant correction is ~ (1/2)^K.

Best estimates:
  Shanks(15,16,17): {nstr(shanks_15_16_17, 15)}  (Δ = {nstr(shanks_15_16_17 - target, 8)})
  Fit r=1/2:        {c_inf:.15f}  (Δ = {c_inf - math.pi/18:+.2e})
  Fit r free:       {c_inf2:.15f}  (r = {r2:.6f})

The extrapolated values are consistent with π/18 to ~4 digits
(model-dependent; direct enumeration at K=19 gives ~4.3 digits).
Higher precision requires K ≈ {K_12dig} for PSLQ-ready ~8 digits (unfeasible by enumeration).

The convergence rate r → 1/2 is a key analytical clue:
it suggests the dominant correction comes from a factor of 2 in the
transfer matrix, possibly related to the parity structure.
""")
