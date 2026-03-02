# Experiments — Paper A (Spectral Theory of Carries)

| Script | Description | Referenced in |
|--------|-------------|---------------|
| `A01_perfactor_identity.py` | Per-factor identity exact verification | Proposition 5 (§6) |
| `A02_constant_c_precision.py` | High-precision boundary coefficient computation | §6.1 (Boundary Coefficients) |
| `A03_constant_c_method.py` | Boundary coefficient fixed methodology | §6.1 (Boundary Coefficients) |
| `A04_c2_decomposition.py` | $c_2$ analytical decomposition | §6 (Correction Series) |
| `A05_markov_transition.py` | Carry Markov transition matrix | §2 (Carry Transition Operator) |
| `A06_diaconis_fulman.py` | Diaconis-Fulman spectrum verification | §2, Corollary 2.1 |
| `A07_rouche_rmax.py` | Rouché theorem analysis for $r_{\max}$ | Conjecture 10 (§8) |
| `A08_boundary_proof.py` | Carry boundary proof | Lemma 7 (§7.1) |
| `A09_gap3_closure.py` | Gap-3 $r_{\max}$ closure argument | §8 (Spectral Radius) |
| `A10_counterexample.py` | Counterexample analysis | §8.3 (Computational Evidence) |
| `A11_rouche_bound.py` | $r_{\max}$ Rouché bound | §8.1 (Rouché Path), Conjecture 10 |
| `A12_anticorrelation.py` | Anti-correlation spectral analysis | Conjecture 4 (§5) |
| `A13_equidistribution.py` | $m$-bit equidistribution lemma | Theorem 2 (§3), Corollaries 2.1–2.4, Proposition 3 |
| `A14_isochrone_mc.c` | Isochrone Monte Carlo profile | §4.2 |
| `A15_mellin_exact_fK.py` | Mellin transform of exact $f_K(T)$ | §4.2 |
| `A16_cascade_closed_forms.py` | Cascade closed-form expressions | §4.2 |
| `A17_precision_extrapolation.py` | Precision extrapolation for convergence rate | §4.2 |
| `A18_final_extrapolation.py` | Final extrapolation for multi-base $\rho(K)$ | §4.2 |
| `A19_prop3_exponent.py` | Proposition 3 exponent verification ($\alpha \approx k-3$) | Proposition 3 (§3.4) |
| `A20_constant_C.py` | Theorem 1 constant $C$ explicit computation | Theorem 1 (§2.2) |
| `A21_rouche_formal.py` | Rouché gap formalization on $\lvert z\rvert=b$ | Conjecture 10 (§8) |
| `A22_adversarial.py` | Adversarial carry profiles — $r_{\max}$ maximization | Conjecture 10 (§8) |
| `A23_boundary_transfer.py` | Boundary transfer theorem for anti-correlation | Conjecture 4 (§5) |

## Shared Utilities

`../src/carry_utils.py` — common functions for carry matrix construction.

## Requirements

Python >= 3.8, NumPy, SciPy, SymPy, mpmath. C compiler for A14.
