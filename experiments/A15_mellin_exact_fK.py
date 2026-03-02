"""
A15: Mellin transform of exact f_K(T) from B29.
Uses the exact (non-MC) isochrone profile to compute M[f_K](s)
and test pole structure vs zeta zeros.
"""
import numpy as np
from scipy import integrate
import math, time, sys

try:
    from mpmath import mp, mpf, zeta as mzeta
    mp.dps = 30
    HAS_MP = True
except ImportError:
    HAS_MP = False

pi = math.pi
C1E = 1 + 3 * math.log(3.0 / 4.0)
SOT = pi / 18 - C1E

ZETA_ZEROS = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351,
              37.5862, 40.9187, 43.3271, 48.0052, 49.7738]


def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def Kker(T):
    if abs(T) < 1e-8:
        return T - 2 * T**3 / 3
    return 2 * (T**2 - math.log(1 + T**2)) / T**3


def MK_series(s, N=500):
    """M[K](s) = 2 sum_{n>=2} (-1)^n / (n(2n+s-3))."""
    total = 0.0
    for n in range(2, N + 2):
        d = n * (2 * n + s - 3)
        if abs(d) < 1e-15:
            continue
        total += (-1)**n / d
    return 2 * total


# Import enumeration function (B29_transfer_operator.py in carry-arithmetic-B-zeta-approximation)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from B29_transfer_operator import enum_exact

NB = 128

pr("=" * 78)
pr("  A15: MELLIN TRANSFORM OF EXACT f_K(T)")
pr("=" * 78)

pr("\n--- Computing exact f_K for K=6..10 ---")
fK_data = {}
for K in range(6, 11):
    r = enum_exact(K, NB)
    fK_data[K] = r
    pr(f"  K={K}: nd={r['nd']}, c1o={r['c1o']:.8f}, So={r['So']:.10f}")


def mellin_fK(sigma, t, fK, Tc, dT):
    """M[f_K](sigma+it) from binned exact data."""
    re_val = 0.0
    im_val = 0.0
    for b in range(len(Tc)):
        if np.isnan(fK[b]) or abs(Tc[b]) < 1e-10:
            continue
        w = fK[b] * Tc[b]**(sigma - 1) * dT
        if abs(t) < 1e-10:
            re_val += w
        else:
            lt = t * math.log(Tc[b])
            re_val += w * math.cos(lt)
            im_val += w * math.sin(lt)
    return complex(re_val, im_val)


pr("\n--- PART A: M[f_K](s) for real s ---\n")
pr(f"  {'s':>6} {'K=6':>12} {'K=8':>12} {'K=10':>12} {'M[K](s)':>12}")
for s in [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
    vals = []
    dT = 1.0 / NB
    for K in [6, 8, 10]:
        r = fK_data[K]
        mf = mellin_fK(s, 0, r['fK'], r['Tc'], dT).real
        vals.append(mf)
    mk = MK_series(s)
    pr(f"  {s:6.1f} {vals[0]:12.8f} {vals[1]:12.8f} {vals[2]:12.8f} {mk:12.8f}")

pr(f"\n  M[K](1) = 2ln2-1 = {2*math.log(2)-1:.10f}")
pr(f"  Sigma_odd_target  = {SOT:.10f}")
pr(f"  If M[f](1)=const, then Sigma = const * M[K](1)")

pr("\n--- PART B: |M[f_K](1/2+it)| at zeta zeros ---\n")
pr(f"  {'t':>8} {'K=8':>12} {'K=10':>12} {'zeta?':>10}")
dT = 1.0 / NB

t_scan = list(np.linspace(0, 55, 56))
for z in ZETA_ZEROS:
    if z not in t_scan:
        t_scan.append(z)
t_scan.sort()

for t in t_scan:
    near = ""
    for z in ZETA_ZEROS:
        if abs(t - z) < 0.5:
            near = f"<- z={z:.1f}"
    mf8 = abs(mellin_fK(0.5, t, fK_data[8]['fK'], fK_data[8]['Tc'], dT))
    mf10 = abs(mellin_fK(0.5, t, fK_data[10]['fK'], fK_data[10]['Tc'], dT))
    if t % 5 < 0.1 or near:
        pr(f"  {t:8.2f} {mf8:12.6f} {mf10:12.6f} {near:>10}")

pr("\n--- PART C: Mellin-Parseval check ---\n")
K_test = 10
r = fK_data[K_test]
dT = 1.0 / NB
direct = sum(r['fK'][b] * Kker(r['Tc'][b]) * dT
             for b in range(NB) if not np.isnan(r['fK'][b]))

sigma_p = 1.5
T_MAX = 30
N_Q = 61
t_q = np.linspace(-T_MAX, T_MAX, N_Q)
dt_q = t_q[1] - t_q[0]
pars = 0.0
for t in t_q:
    mf = mellin_fK(sigma_p, t, r['fK'], r['Tc'], dT)
    mk_re = sum(Kker(r['Tc'][b]) * r['Tc'][b]**(sigma_p-1) *
                math.cos(-t * math.log(r['Tc'][b])) * dT
                for b in range(NB) if r['Tc'][b] > 1e-10)
    mk_im = sum(Kker(r['Tc'][b]) * r['Tc'][b]**(sigma_p-1) *
                math.sin(-t * math.log(r['Tc'][b])) * dT
                for b in range(NB) if r['Tc'][b] > 1e-10)
    mk = complex(mk_re, mk_im)
    pars += (mf * mk).real * dt_q
pars /= (2 * pi)

pr(f"  K={K_test}: direct int f*K = {direct:.10f}")
pr(f"  Parseval (sigma={sigma_p}) = {pars:.10f}")
pr(f"  Ratio = {pars/direct if abs(direct) > 1e-15 else 0:.6f}")

pr("\n--- PART D: Convergence of M[f_K](s) ---\n")
pr(f"  Test: does M[f_K](s) converge as K grows?")
for s_val in [1.0, 2.0]:
    pr(f"\n  s = {s_val}:")
    vals = []
    for K in sorted(fK_data.keys()):
        r = fK_data[K]
        mf = mellin_fK(s_val, 0, r['fK'], r['Tc'], dT).real
        vals.append((K, mf))
        pr(f"    K={K:2d}: M[f_K]({s_val}) = {mf:+.10f}")
    if len(vals) >= 3:
        d1 = vals[-1][1] - vals[-2][1]
        d2 = vals[-2][1] - vals[-3][1]
        if abs(d2) > 1e-12:
            rho = d1 / d2
            extrap = vals[-1][1] + d1 * rho / (1 - rho) if abs(1-rho) > 0.01 else vals[-1][1]
            pr(f"    rho = {rho:.4f}, M[f](s) extrap = {extrap:+.10f}")

pr("\n--- PART E: Product M[f]*M[K] structure ---\n")
pr(f"  If Sigma_odd = (1/2pi) int M[f]*M[K] dt (Parseval),")
pr(f"  then the integrand structure reveals how 2/9 arises.")
pr(f"\n  Integrand |M[f]*M[K]| along sigma=1.5:")
for t in [0, 1, 2, 5, 10, 14.13, 20, 25.01, 30]:
    mf = mellin_fK(1.5, t, fK_data[10]['fK'], fK_data[10]['Tc'], dT)
    mk_re = sum(Kker(fK_data[10]['Tc'][b]) * fK_data[10]['Tc'][b]**(0.5) *
                math.cos(-t * math.log(fK_data[10]['Tc'][b])) * dT
                for b in range(NB) if fK_data[10]['Tc'][b] > 1e-10)
    mk_im = sum(Kker(fK_data[10]['Tc'][b]) * fK_data[10]['Tc'][b]**(0.5) *
                math.sin(-t * math.log(fK_data[10]['Tc'][b])) * dT
                for b in range(NB) if fK_data[10]['Tc'][b] > 1e-10)
    mk = complex(mk_re, mk_im)
    prod = mf * mk
    pr(f"  t={t:7.2f}: |M[f]*M[K]| = {abs(prod):.8f}, Re = {prod.real:+.8f}")

pr("\n" + "=" * 78)
pr("  END A15")
pr("=" * 78)
