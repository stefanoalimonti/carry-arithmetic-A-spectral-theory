# carry-arithmetic-A-spectral-theory

**Spectral Theory of Carries in Positional Multiplication**

*Author: Stefano Alimonti* · [ORCID 0009-0009-1183-1698](https://orcid.org/0009-0009-1183-1698)

## Main Result

The $m$-Bit Equidistribution Lemma (Theorem 2): if the convolution sum at a given position contains $m$ independent uniform-mod-$b$ additive components, the carry transfer operator has exact eigenvalues $\{1, 1/b, \ldots, 1/b^m\}$. For multiplication at position $j \geq 2$, two free digit bits yield $m = 2$ and universal eigenvalues $\{1, 1/b, 1/b^2\}$; higher eigenvalues converge exponentially to $1/b^k$ as $j \to \infty$ (Proposition 3).

## Status

Complete. Ready for submission.

## Repository Structure

```
paper/carry_spectral_theory.md       The paper
experiments/
  A01_perfactor_identity.py          Per-factor identity (Proposition 5, §6)
  A02_constant_c_precision.py        Trace anomaly high-precision (§6.1)
  A03_constant_c_method.py           Trace anomaly methodology (§6.1)
  A04_c2_decomposition.py            c₂ analytical decomposition (§6)
  A05_markov_transition.py           Carry Markov transition matrix (§2)
  A06_diaconis_fulman.py             Diaconis-Fulman verification (Corollary 2.1)
  A07_rouche_rmax.py                 Rouché analysis for r_max (§8)
  A08_boundary_proof.py              Carry boundary structure (§7.1)
  A09_gap3_closure.py                Gap-3 r_max closure (Corollary 9)
  A10_counterexample.py              Counterexample search for r_max > b (Conjecture 10)
  A11_rouche_bound.py                Rouché bound on |z|=b (§8.1)
  A12_anticorrelation.py             Anti-correlation spectral (Conjecture 4)
  A13_equidistribution.py            m-Bit Equidistribution Lemma (Theorem 2, Corollaries 2.1–2.4, Proposition 3)
  A14_isochrone_mc.c                 Isochrone Monte Carlo profile (§4.2)
  A15_mellin_exact_fK.py             Mellin transform of exact f_K(T) (§4.2)
  A16_cascade_closed_forms.py        Cascade closed-form expressions (§4.2)
  A17_precision_extrapolation.py     Precision extrapolation for convergence rate (§4.2)
  A18_final_extrapolation.py         Final extrapolation for multi-base ρ(K) (§4.2)
  A19_prop3_exponent.py              Proposition 3 exponent verification (§3.4)
  A20_constant_C.py                  Theorem 1 constant C computation (§2.2)
  A21_rouche_formal.py               Rouché gap formalization on |z|=b (§8)
  A22_adversarial.py                 Adversarial carry profiles — r_max maximization (§8)
  A23_boundary_transfer.py           Boundary transfer theorem for anti-correlation (§5)
src/
  carry_utils.py                     Shared utility functions
```

## Reproduction

```bash
pip install numpy scipy sympy mpmath
python experiments/A01_perfactor_identity.py
# ... through A23
```

## Dependencies

- Python >= 3.8, NumPy, SciPy, SymPy, mpmath

## Companion Papers

| Label | Title | Repository |
|-------|-------|------------|
| [C] | Eigenvalue Statistics of Carry Companion Matrices: Markov-Driven GOE↔GUE Transition | [`carry-arithmetic-C-matrix-statistics`](https://github.com/stefanoalimonti/carry-arithmetic-C-matrix-statistics) |
| [D] | The Carry-Zero Entropy Bound | [`carry-arithmetic-D-factorization-limits`](https://github.com/stefanoalimonti/carry-arithmetic-D-factorization-limits) |
| [E] | The Trace Anomaly of Binary Multiplication | [`carry-arithmetic-E-trace-anomaly`](https://github.com/stefanoalimonti/carry-arithmetic-E-trace-anomaly) |
| [F] | Exact Covariance Structure | [`carry-arithmetic-F-covariance-structure`](https://github.com/stefanoalimonti/carry-arithmetic-F-covariance-structure) |
| [G] | The Angular Uniqueness of Base 2 | [`carry-arithmetic-G-angular-uniqueness`](https://github.com/stefanoalimonti/carry-arithmetic-G-angular-uniqueness) |
| [B] | Carry Polynomials and the Euler Product | [`carry-arithmetic-B-zeta-approximation`](https://github.com/stefanoalimonti/carry-arithmetic-B-zeta-approximation) |
| [P2] | The Sector Ratio in Binary Multiplication | [`carry-arithmetic-P2-sector-ratio`](https://github.com/stefanoalimonti/carry-arithmetic-P2-sector-ratio) |

### Citation

```bibtex
@article{alimonti2026spectral_theory,
  author  = {Alimonti, Stefano},
  title   = {Spectral Theory of Carries in Positional Multiplication},
  year    = {2026},
  note    = {Preprint},
  url     = {https://github.com/stefanoalimonti/carry-arithmetic-A-spectral-theory}
}
```

## License

Paper: CC BY 4.0. Code: MIT License.
