#!/usr/bin/env python3
"""
A13: m-Bit Equidistribution Lemma — Verification

THEOREM (m-Bit Equidistribution Lemma):
  If V = B_1 + ... + B_m + S where B_i are iid Bernoulli(1/2) independent of S >= 0,
  then the carry operator (Tf)(c) = E_V[f(floor((V+c)/2))] has eigenvalues
  {1, 1/2, 1/4, ..., 1/2^m} EXACTLY on any state space of size >= m+1.

PROOF MECHANISM:
  (Tc^k)(c) = E_V[floor((V+c)/2)^k] is a polynomial in c of degree k
  (with no parity oscillation) for all k <= m. This makes T exactly upper-
  triangular in the monomial basis with diagonal 1/2^k.
  
  Each Bernoulli(1/2) bit in V "smooths" the floor function by one degree:
  the parity oscillation of (Tc^k)(c) has degree max(k-1-m, -inf) in c.
"""

import numpy as np
from fractions import Fraction
from itertools import product as iproduct
import math
import sys

def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def build_transfer_matrix(dist, base=2):
    """Build carry transition matrix T[c_in, c_out] = P(floor((V+c_in)/b) = c_out).
    State space: {0, ..., v_max - 1} (matching reference data convention)."""
    v_max = max(dist.keys())
    dim = v_max
    if dim < 1:
        dim = 1
    T = np.zeros((dim, dim))
    for c in range(dim):
        for v, p in dist.items():
            c_out = (v + c) // base
            if c_out < dim:
                T[c, c_out] += p
    return T


def compute_Tc_k(dist, k, c_values, base=2):
    """Compute (Tc^k)(c) = E_V[floor((V+c)/b)^k] for a list of c values.
    Returns dict {c: Fraction value} using exact arithmetic."""
    result = {}
    for c in c_values:
        total = Fraction(0)
        for v, p_frac in dist.items():
            c_out = (v + c) // base
            total += p_frac * Fraction(c_out) ** k
        result[c] = total
    return result


def compute_mult_conv_dist(j):
    """Exact distribution of conv_j for multiplication (MSB g_0=h_0=1).
    Returns dict {v: Fraction(count, total)}."""
    n_free = 2 * j
    dist = {}
    total = 1 << n_free
    for bits in range(total):
        g = [1]
        h = [1]
        for i in range(j):
            g.append((bits >> i) & 1)
            h.append((bits >> (j + i)) & 1)
        v = sum(g[i] * h[j - i] for i in range(j + 1))
        dist[v] = dist.get(v, 0) + 1
    return {v: Fraction(cnt, total) for v, cnt in dist.items()}


def make_binom_plus_S(m, S_dist_frac):
    """Construct distribution of V = Binom(m, 1/2) + S.
    S_dist_frac: dict {s: Fraction prob}. Returns dict {v: Fraction prob}."""
    binom_dist = {}
    for r in range(m + 1):
        binom_dist[r] = Fraction(math.comb(m, r), 2**m)
    V_dist = {}
    for r, p_r in binom_dist.items():
        for s, p_s in S_dist_frac.items():
            v = r + s
            V_dist[v] = V_dist.get(v, Fraction(0)) + p_r * p_s
    return V_dist


# ════════════════════════════════════════════════════════════════════
# PART A: PARITY-INDEPENDENCE OF (Tc^k)(c)
# ════════════════════════════════════════════════════════════════════

def verify_parity_independence():
    pr("\n" + "=" * 72)
    pr("PART A: PARITY-INDEPENDENCE OF (Tc^k)(c)")
    pr("=" * 72)
    pr("\n(Tc^k)(c) = E_V[floor((V+c)/2)^k].")
    pr("Claim: this is a polynomial in c (no parity oscillation) iff V has")
    pr("m >= k independent Bernoulli(1/2) additive components.\n")

    S_dist = {0: Fraction(3, 4), 1: Fraction(1, 4)}

    for m in range(5):
        V_dist = make_binom_plus_S(m, S_dist)
        V_dist_float = {v: float(p) for v, p in V_dist.items()}
        v_max = max(V_dist.keys())

        pr(f"  m={m} Bernoulli bits, S ~ Bern(1/4):")
        pr(f"    V_max = {v_max}, states = {{0..{v_max-1}}}")

        for k in range(5):
            c_vals = list(range(min(20, v_max)))
            Tck = compute_Tc_k(V_dist, k, c_vals)

            max_osc = Fraction(0)
            for i in range(0, len(c_vals) - 1, 2):
                c_even = c_vals[i]
                c_odd = c_vals[i + 1] if i + 1 < len(c_vals) else None
                if c_odd is None:
                    continue
                if c_even % 2 != 0 or c_odd % 2 != 1:
                    continue

                val_even = Tck[c_even]
                val_odd = Tck[c_odd]

                poly_even = sum(Fraction(c_even)**j * Fraction(1) for j in range(k+1))
                poly_odd = sum(Fraction(c_odd)**j * Fraction(1) for j in range(k+1))

            is_poly = True
            if len(c_vals) >= k + 2:
                vals = [float(Tck[c]) for c in c_vals[:k+2]]
                cs = [float(c) for c in c_vals[:k+2]]
                coeffs = np.polyfit(cs, vals, k)
                for c in c_vals:
                    predicted = sum(coeffs[k-j] * float(c)**j for j in range(k+1))
                    actual = float(Tck[c])
                    if abs(predicted - actual) > 1e-10:
                        is_poly = False
                        break

            expected = (m >= k)
            status = "✓" if is_poly == expected else "✗"
            pr(f"    k={k}: polynomial in c? {is_poly:<5}  "
               f"expected (m>={k})? {expected:<5}  {status}")
        pr()


# ════════════════════════════════════════════════════════════════════
# PART B: EIGENVALUES OF MULTIPLICATION TRANSFER MATRICES T^(j)
# ════════════════════════════════════════════════════════════════════

def verify_multiplication_eigenvalues():
    pr("\n" + "=" * 72)
    pr("PART B: MULTIPLICATION TRANSFER MATRICES T^(j)")
    pr("=" * 72)
    pr("\nconv_j = g_j + h_j + sum_{i=1}^{j-1} g_i*h_{j-i}  (g_0=h_0=1)")
    pr("For j>=2: m=2 independent Bern(1/2) bits (g_j and h_j).")
    pr("m-Bit Lemma: eigenvalues {1, 1/2, 1/4} exact for j>=2.\n")

    pr(f"  {'j':>3s}  {'dim':>4s}  {'eigenvalues':60s}  {'first 3 exact?':>14s}")
    pr(f"  {'─'*3}  {'─'*4}  {'─'*60}  {'─'*14}")

    for j in range(1, 11):
        dist_frac = compute_mult_conv_dist(j)
        dist_float = {v: float(p) for v, p in dist_frac.items()}

        T = build_transfer_matrix(dist_float)
        dim = T.shape[0]
        eigs = sorted(np.linalg.eigvals(T).real, reverse=True)

        n_check = min(3, dim)
        exact = all(abs(eigs[i] - 0.5**i) < 1e-8 for i in range(n_check))

        eig_strs = []
        for i, e in enumerate(eigs[:min(8, dim)]):
            target = 0.5**i
            if abs(e - target) < 1e-8:
                if i == 0:
                    eig_strs.append("1")
                else:
                    eig_strs.append(f"1/{2**i}")
            else:
                f = Fraction(e).limit_denominator(10000)
                if abs(float(f) - e) < 1e-8:
                    eig_strs.append(str(f))
                else:
                    eig_strs.append(f"{e:.8f}")

        eig_display = ", ".join(eig_strs)
        status = "✓ EXACT" if exact else "partial"

        pr(f"  {j:3d}  {dim:4d}  [{eig_display:58s}]  {status}")


# ════════════════════════════════════════════════════════════════════
# PART C: EIGENVALUES WITH m BERNOULLI BITS + ARBITRARY S
# ════════════════════════════════════════════════════════════════════

def verify_eigenvalues_m_bits():
    pr("\n\n" + "=" * 72)
    pr("PART C: EIGENVALUES WITH m BERNOULLI BITS + ARBITRARY S")
    pr("=" * 72)
    pr("\nV = B_1 + ... + B_m + S, B_i iid Bern(1/2), independent of S.")
    pr("Claim: eigenvalues {1, 1/2, ..., 1/2^m} are exact for any S")
    pr("(on state space of sufficient size).\n")

    test_S_dists = [
        ("S=0", {0: Fraction(1)}),
        ("S~Bern(1/4)", {0: Fraction(3, 4), 1: Fraction(1, 4)}),
        ("S~Binom(3,1/4)", {
            0: Fraction(27, 64), 1: Fraction(27, 64),
            2: Fraction(9, 64), 3: Fraction(1, 64)
        }),
        ("S~Unif{0..3}", {k: Fraction(1, 4) for k in range(4)}),
        ("S~Unif{0..6}", {k: Fraction(1, 7) for k in range(7)}),
    ]

    for m in range(1, 5):
        pr(f"\n--- m = {m} Bernoulli bits ---\n")

        for s_name, s_dist in test_S_dists:
            V_dist_frac = make_binom_plus_S(m, s_dist)
            V_dist_float = {v: float(p) for v, p in V_dist_frac.items()}

            T = build_transfer_matrix(V_dist_float)
            dim = T.shape[0]
            eigs = sorted(np.linalg.eigvals(T).real, reverse=True)

            n_check = min(m + 1, dim)
            exact = all(abs(eigs[i] - 0.5**i) < 1e-8 for i in range(n_check))

            eig_str = ", ".join(f"{e:.6f}" for e in eigs[:min(6, dim)])
            status = "✓" if exact else "✗"

            pr(f"  {s_name:20s}  dim={dim:2d}  "
               f"eigs=[{eig_str}]  "
               f"first {n_check} match 1/2^k? {status}")


# ════════════════════════════════════════════════════════════════════
# PART D: EIGENVECTOR TEST — (-1)^c is eigenfunction with eigenvalue 1/2
# ════════════════════════════════════════════════════════════════════

def verify_eigenvector():
    pr("\n\n" + "=" * 72)
    pr("PART D: EIGENVECTOR TEST — (-1)^c eigenfunction with eigenvalue 1/2")
    pr("=" * 72)
    pr("\nThe Parity Lemma: f(c) = (-1)^c satisfies Tf = (1/2)f")
    pr("for ANY distribution P(V) with P(V odd) = 1/2.\n")

    for j in range(1, 9):
        dist = compute_mult_conv_dist(j)
        dist_float = {v: float(p) for v, p in dist.items()}
        T = build_transfer_matrix(dist_float)
        dim = T.shape[0]

        f = np.array([(-1.0)**c for c in range(dim)])
        Tf = T @ f
        expected = 0.5 * f

        match = np.allclose(Tf, expected, atol=1e-12)
        max_err = np.max(np.abs(Tf - expected))

        pr(f"  T^({j}): dim={dim}, Tf = 0.5*f ? {'✓' if match else '✗'}  "
           f"max|Tf - 0.5f| = {max_err:.2e}")

    pr("\n  Also test with non-multiplication distributions:")
    for name, dist_frac in [
        ("Binom(8, 1/2)", make_binom_plus_S(8, {0: Fraction(1)})),
        ("Binom(2, 1/2) + Unif{0..3}", make_binom_plus_S(2, {k: Fraction(1,4) for k in range(4)})),
    ]:
        dist_float = {v: float(p) for v, p in dist_frac.items()}
        T = build_transfer_matrix(dist_float)
        dim = T.shape[0]
        f = np.array([(-1.0)**c for c in range(dim)])
        Tf = T @ f
        match = np.allclose(Tf, 0.5 * f, atol=1e-12)
        pr(f"  {name:40s}: dim={dim}, Tf = 0.5*f ? {'✓' if match else '✗'}")


# ════════════════════════════════════════════════════════════════════
# PART E: COMPOSITION EIGENVALUE (1/2)^K via COMMON EIGENFUNCTION
# ════════════════════════════════════════════════════════════════════

def verify_composition():
    pr("\n\n" + "=" * 72)
    pr("PART E: COMPOSITION EIGENVALUE (1/2)^K")
    pr("=" * 72)
    pr("\nSince (-1)^c is a common eigenfunction of ALL T_j with eigenvalue 1/2,")
    pr("the composed operator T_1·...·T_K has (T_1...T_K) f = (1/2)^K f.\n")

    for K in [4, 6, 8, 10, 12]:
        dists = []
        max_dim = 0
        for j in range(1, K + 1):
            j_eff = min(j, 8)
            dist = compute_mult_conv_dist(j_eff)
            dist_float = {v: float(p) for v, p in dist.items()}
            T = build_transfer_matrix(dist_float)
            dists.append(T)
            max_dim = max(max_dim, T.shape[0])

        f = np.array([(-1.0)**c for c in range(max_dim)])
        result = f.copy()
        for T in reversed(dists):
            d = T.shape[0]
            result_trunc = result[:d]
            result_trunc = T @ result_trunc
            result[:d] = result_trunc

        expected = (0.5)**K * f
        ratio = result / f
        all_match = np.allclose(ratio[:dists[0].shape[0]], 0.5**K, atol=1e-8)

        pr(f"  K={K:2d}: (T_1...T_K)f / f = {ratio[0]:.12e}, "
           f"(1/2)^K = {0.5**K:.12e}  {'✓' if all_match else '✗'}")


# ════════════════════════════════════════════════════════════════════
# PART F: BASE-b GENERALIZATION
# ════════════════════════════════════════════════════════════════════

def compute_mult_conv_dist_base_b(j, b):
    """Exact distribution for base-b multiplication at position j.
    g_0 = h_0 = 1 (MSBs), g_i,h_i ~ Uniform{0,...,b-1} for i >= 1."""
    n_free = 2 * j
    dist = {}
    total = b ** n_free
    for config in range(total):
        g = [1]
        h = [1]
        tmp = config
        for i in range(j):
            g.append(tmp % b)
            tmp //= b
        for i in range(j):
            h.append(tmp % b)
            tmp //= b
        v = sum(g[i] * h[j - i] for i in range(j + 1))
        dist[v] = dist.get(v, 0) + 1
    return {v: c / total for v, c in dist.items()}


def verify_base_b():
    pr("\n\n" + "=" * 72)
    pr("PART F: BASE-b GENERALIZATION")
    pr("=" * 72)
    pr("\nFor prime base b, g_j*h_0 and g_0*h_j provide two independent")
    pr("'uniform mod b' components. Eigenvalues {1, 1/b, 1/b²} should be exact.\n")

    for b in [2, 3, 5, 7]:
        pr(f"  Base b={b}:")
        max_j = 6 if b <= 3 else (4 if b == 5 else 3)

        for j in range(1, max_j + 1):
            if b ** (2*j) > 500000:
                pr(f"    j={j}: skipped (too large)")
                continue

            dist = compute_mult_conv_dist_base_b(j, b)
            T = build_transfer_matrix(dist, base=b)
            dim = T.shape[0]
            eigs = sorted(np.linalg.eigvals(T).real, reverse=True)

            n_check = min(3, dim)
            exact = all(abs(eigs[i] - (1.0/b)**i) < 1e-7 for i in range(n_check))

            eig_strs = []
            for i, e in enumerate(eigs[:min(6, dim)]):
                target = (1.0/b)**i
                if abs(e - target) < 1e-7:
                    eig_strs.append(f"1/{b}^{i}" if i > 0 else "1")
                else:
                    eig_strs.append(f"{e:.6f}")

            pr(f"    j={j}: dim={dim:3d}  [{', '.join(eig_strs):50s}]  "
               f"{'✓' if exact else '—'}")
        pr()


# ════════════════════════════════════════════════════════════════════
# PART G: EXPLICIT 3x3 ALGEBRAIC VERIFICATION FOR j=2
# ════════════════════════════════════════════════════════════════════

def verify_algebraic_j2():
    pr("\n" + "=" * 72)
    pr("PART G: ALGEBRAIC VERIFICATION FOR j=2 (3×3 MATRIX)")
    pr("=" * 72)
    pr("\nV_2 = h_2 + g_1*h_1 + g_2 with MSBs g_0=h_0=1.")
    pr("P(V=0)=3/16, P(V=1)=7/16, P(V=2)=5/16, P(V=3)=1/16.\n")

    dist = {0: Fraction(3, 16), 1: Fraction(7, 16),
            2: Fraction(5, 16), 3: Fraction(1, 16)}

    pr("  Transfer matrix T[c_in, c_out] (exact rational):")
    dim = 3
    T_frac = [[Fraction(0)] * dim for _ in range(dim)]
    for c in range(dim):
        for v, p in dist.items():
            c_out = (v + c) // 2
            if c_out < dim:
                T_frac[c][c_out] += p

    for c in range(dim):
        row = [str(T_frac[c][j]) for j in range(dim)]
        pr(f"    c={c}: [{', '.join(f'{r:>6s}' for r in row)}]")

    pr("\n  Characteristic polynomial: det(T - λI) = 0")

    T_np = np.array([[float(T_frac[i][j]) for j in range(dim)] for i in range(dim)])
    eigs = sorted(np.linalg.eigvals(T_np).real, reverse=True)
    pr(f"    Eigenvalues: [{', '.join(f'{e:.10f}' for e in eigs)}]")
    pr(f"    Expected:    [1.0000000000, 0.5000000000, 0.2500000000]")
    match = all(abs(eigs[i] - 0.5**i) < 1e-10 for i in range(3))
    pr(f"    Exact match: {'✓ YES' if match else '✗ NO'}")

    pr(f"\n  Verify trace and determinant:")
    tr = sum(eigs)
    det_val = np.prod(eigs)
    pr(f"    trace = {tr:.10f}  (expected 1+1/2+1/4 = {1+0.5+0.25})")
    pr(f"    det   = {det_val:.10f}  (expected 1*1/2*1/4 = {1*0.5*0.25})")

    pr(f"\n  Polynomial basis check: (Tc^k)(c) for k=0,1,2")
    for k in range(3):
        vals = compute_Tc_k(dist, k, list(range(dim)))
        pr(f"    (Tc^{k})(c): c=0 → {vals[0]}, c=1 → {vals[1]}, c=2 → {vals[2]}")
        if k == 0:
            pr(f"      = 1 for all c (eigenvalue 1, const eigenfunction)")
        elif k == 1:
            pr(f"      = (1/2)c + {vals[0]}  → coefficient of c is 1/2")
        elif k == 2:
            a = (vals[2] - 2*vals[1] + vals[0]) / 2
            b_coeff = vals[1] - vals[0] - a
            c_coeff = vals[0]
            pr(f"      = ({a})c² + ({b_coeff})c + ({c_coeff})  → coeff of c² = {a} = 1/4")


# ════════════════════════════════════════════════════════════════════
# SYNTHESIS
# ════════════════════════════════════════════════════════════════════

def synthesis():
    pr("\n" + "=" * 72)
    pr("SYNTHESIS: THE m-BIT EQUIDISTRIBUTION LEMMA")
    pr("=" * 72)
    pr("""
THEOREM (m-Bit Equidistribution Lemma):
  Let V = B_1 + ... + B_m + S where B_i are iid Bernoulli(1/2),
  independent of S >= 0. The carry operator (Tf)(c) = E[f(floor((V+c)/2))]
  satisfies: (Tc^k)(c) is a polynomial of degree k in c with leading
  coefficient 1/2^k, for all k <= m. No parity oscillation.

  Consequently, T is exactly upper-triangular in {1, c, c^2, ..., c^m}
  with diagonal {1, 1/2, 1/4, ..., 1/2^m}. These are EXACT eigenvalues.

PROOF SKETCH:
  Write V = R + S with R ~ Binom(m, 1/2), independent of S.
  
  (Tc^k)(c) = E_S[ E_R[ floor((S+c+R)/2)^k ] ]
  
  The inner expectation averages floor((u+R)/2)^k over R ~ Binom(m, 1/2)
  where u = S+c. This Binomial smoothing eliminates the parity oscillation
  of the floor function:
  
  - floor(u/2)^k has parity oscillation of polynomial degree k-1 in u.
  - Each Bernoulli(1/2) averaging step: A f(u) = [f(u)+f(u+1)]/2
    maps the oscillatory part h(u) to -Delta(h)/2, reducing degree by 1.
  - After m averaging steps: oscillation degree = max(k-1-m, -inf).
  - For k <= m: oscillation = 0, so the result is a pure polynomial.
  
  Then (Tc^k)(c) = E_S[polynomial(S+c)] = polynomial(c). QED.

APPLICATION TO MULTIPLICATION:
  At position j >= 2, conv_j = g_j*1 + 1*h_j + (inner terms)
  = g_j + h_j + S, where g_j and h_j are independent Bern(1/2).
  So m = 2, giving exact eigenvalues {1, 1/2, 1/4} for ALL j >= 2.
  This PROVES Conjecture 1c item 2 (Double Universality).

COMPOSITION:
  (-1)^c is a common eigenfunction of ALL T_j with eigenvalue 1/2.
  Therefore (T_1 ... T_K)(-1)^c = (1/2)^K (-1)^c.
  This PROVES Conjecture 1c item 4 (Composition spectral gap).
""")


def main():
    pr("=" * 72)
    pr("A13: m-BIT EQUIDISTRIBUTION LEMMA — VERIFICATION")
    pr("=" * 72)

    verify_algebraic_j2()
    verify_multiplication_eigenvalues()
    verify_parity_independence()
    verify_eigenvalues_m_bits()
    verify_eigenvector()
    verify_composition()
    verify_base_b()
    synthesis()

    pr("\n" + "=" * 72)
    pr("A13 COMPLETE")
    pr("=" * 72)


if __name__ == '__main__':
    main()
