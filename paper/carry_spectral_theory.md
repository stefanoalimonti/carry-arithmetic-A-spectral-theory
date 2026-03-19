# Spectral Theory of Carries in Positional Multiplication

**Author:** Stefano Alimonti
**Affiliation:** Independent Researcher
**Date:** March 2026

---

## Abstract

Diaconis and Fulman [1] proved that carries in base-$b$ addition form a Markov chain with eigenvalues $1/b^k$. We extend this to multiplication. Our main result is the $m$-Bit Equidistribution Lemma: in the binary case proved here, if the convolution sum contains $m$ independent uniform-mod-$b$ additive components, the carry transfer operator has eigenvalues $\lbrace{}1, 1/b, \ldots, 1/b^m\rbrace{}$ exactly, via a Bernoulli-smoothing argument on the monomial basis; the general prime-base extension is given in proof-sketch form in §3.2. For multiplication at position $j \geq 2$, two free digit bits yield $m = 2$ and universal eigenvalues $\lbrace{}1, 1/b, 1/b^2\rbrace{}$; the remaining $j-1$ product terms provide exponential damping at rate $1/b$ per term, giving $\lambda_k^{(j)} = 1/b^k + O(b^{-j})$. Corollaries include: composition spectral gap $(1/b)^K$ for the $K$-position operator, total variation mixing in $O(\log c_{\max})$ steps, and a Gershgorin perturbation bound for general input distributions. All results are verified with exact rational arithmetic for bases 2, 3, 5, 7.

**Keywords:** carries, Markov chains, Diaconis-Fulman, positional multiplication, spectral gap, equidistribution

**MSC 2020:** 60J10, 11A63, 15A18, 60C05

---

## 1. Introduction

When two integers are multiplied by hand in base $b$, one writes out partial products and sums them column by column, propagating carries to the left. This is the same algorithm taught in primary school, and it is the process we analyze. The central object is the carry sequence $c_0, c_1, c_2, \ldots$ that propagates through the columns. We show that this sequence forms a Markov chain whose spectral properties are surprisingly rigid.

Diaconis and Fulman [1] proved that carries in base-$b$ *addition* form a Markov chain with eigenvalues $1/b^k$. We extend this to multiplication and develop the full spectral theory.

### 1.1 A Small Example: $13 \times 11$ in Base 2

To ground the discussion, we work through a complete example. Consider multiplying $13 \times 11 = 143$ in binary. The digits are:

$$13 = (1101)_2, \quad g_3 = 1,\; g_2 = 1,\; g_1 = 0,\; g_0 = 1$$
$$11 = (1011)_2, \quad h_3 = 1,\; h_2 = 0,\; h_1 = 1,\; h_0 = 1$$

The grade-school multiplication algorithm writes out one partial product for each digit of $h$, shifted by the appropriate number of positions, and sums the columns:

```
          1 1 0 1          (13)
        × 1 0 1 1          (11)
        ---------
          1 1 0 1          h₀ = 1:  13 × 1
        1 1 0 1 ·          h₁ = 1:  13 × 1, shifted left 1
      0 0 0 0 · ·          h₂ = 0:  13 × 0, shifted left 2
    1 1 0 1 · · ·          h₃ = 1:  13 × 1, shifted left 3
    -----------------
    1 0 0 0 1 1 1 1        (143)
```

Reading each column from bottom to top, the *column sum* at position $k$ is

$$\text{conv}_k = \sum_{i+j=k} g_i \, h_j$$

which is the discrete convolution of the digit sequences — exactly the same convolution that appears in polynomial multiplication. The carry at each position is then determined by a simple recursion: divide the column sum plus incoming carry by the base $b$; the remainder is the output digit, and the quotient is the outgoing carry:

$$c_{k+1} = \left\lfloor \frac{\text{conv}_k + c_k}{b} \right\rfloor, \qquad f_k = (\text{conv}_k + c_k) \bmod b$$

Here is the complete trace for our example ($b = 2$, $c_0 = 0$):

| Pos $k$ | Column sum $\text{conv}_k$ | Terms | Carry in $c_k$ | Total | Digit $f_k$ | Carry out $c_{k+1}$ |
|---|---|---|---|---|---|---|
| 0 | $g_0 h_0 = 1$ | 1 term | 0 | 1 | 1 | 0 |
| 1 | $g_1 h_0 + g_0 h_1 = 0 + 1 = 1$ | 2 terms | 0 | 1 | 1 | 0 |
| 2 | $g_2 h_0 + g_1 h_1 + g_0 h_2 = 1 + 0 + 0 = 1$ | 3 terms | 0 | 1 | 1 | 0 |
| 3 | $g_3 h_0 + g_2 h_1 + g_1 h_2 + g_0 h_3 = 1 + 1 + 0 + 1 = 3$ | 4 terms | 0 | **3** | 1 | **1** |
| 4 | $g_3 h_1 + g_2 h_2 + g_1 h_3 = 1 + 0 + 0 = 1$ | 3 terms | **1** | 2 | 0 | **1** |
| 5 | $g_3 h_2 + g_2 h_3 = 0 + 1 = 1$ | 2 terms | **1** | 2 | 0 | **1** |
| 6 | $g_3 h_3 = 1$ | 1 term | **1** | 2 | 0 | **1** |

The output digits $(f_0, \ldots, f_6, c_7) = (1,1,1,1,0,0,0,1)$ give $143 = (10001111)_2$.

**The carry chain as a Markov process.** The carry sequence $c_0 = 0, 0, 0, 0, 1, 1, 1, 1$ is the trajectory of a Markov chain: the next carry $c_{k+1}$ depends on the current carry $c_k$ and the random column sum $\text{conv}_k$, but not on earlier carries. The figure below traces the carry propagation through the multiplication:

```
Position:   0      1      2      3      4      5      6
           ┌──┐   ┌──┐   ┌──┐   ┌──┐   ┌──┐   ┌──┐   ┌──┐
conv_k:    │ 1│   │ 1│   │ 1│   │ 3│   │ 1│   │ 1│   │ 1│
           └──┘   └──┘   └──┘   └──┘   └──┘   └──┘   └──┘
c_k:  0 ───► 0 ───► 0 ───► 0 ───► 1 ───► 1 ───► 1 ───► 1
                                   ▲
                              first carry
           ┌──┐   ┌──┐   ┌──┐   ┌──┐   ┌──┐   ┌──┐   ┌──┐
f_k:       │ 1│   │ 1│   │ 1│   │ 1│   │ 0│   │ 0│   │ 0│
           └──┘   └──┘   └──┘   └──┘   └──┘   └──┘   └──┘
```

At position 3, the column sum reaches 3 — exceeding the base — and the first carry is born. Once present, it persists through subsequent positions. This persistence is governed by the *spectral gap* of the carry transfer operator: the second eigenvalue $\lambda_1 = 1/b$ controls how quickly the chain forgets its initial state.

**Why multiplication differs from addition.** In addition of $n+1$ numbers, each column sum is $\sum_{i=1}^{n+1} a_i^{(k)}$ — a sum of $n+1$ independent digits, all uniform modulo $b$. The reachable carry chain then has the classical Diaconis-Fulman eigenvalues $\lbrace{}1, 1/b, \ldots, 1/b^n\rbrace{}$ [1]. In multiplication, the column sum $\text{conv}_k = \sum_{i+j=k} g_i h_j$ has a different structure. In the local position-$j$ analysis of Part I — equivalently, for interior columns where both endpoint terms are present — there are two distinguished summands, $g_j h_0$ and $g_0 h_j$, that are independent uniform-mod-$b$ digits (since $g_0, h_0 \geq 1$ are coprime to $b$ for prime $b$). The remaining terms are products $g_i h_j$, which are *not* uniformly distributed. Our main result (Theorem 2, the $m$-Bit Equidistribution Lemma) shows that these two uniform components determine the first three eigenvalues $\lbrace{}1, 1/b, 1/b^2\rbrace{}$ exactly, while the product terms provide exponential damping of higher eigenvalues.

### 1.2 Prior Work on Multiplication Carries

Diaconis and Fulman [1, §6] briefly discuss the extension of their carry analysis to multiplication, observing that "carries in multiplication are more complex" and that the convolution structure replaces the simple column sum of addition. Holte [4] similarly restricts attention to carries in addition. Neither paper develops the spectral theory of the multiplication operator — the column sums in multiplication have fewer uniform components, the state space grows with position, and the transfer operator is position-dependent. The present paper addresses all three difficulties: the $m$-Bit Equidistribution Lemma (Theorem 2) handles the reduced uniformity, the Gershgorin perturbation bound (Theorem 1) controls non-uniform inputs, and the exponential convergence result (Proposition 3) shows that position-dependence vanishes rapidly.

The connection to *riffle shuffling* deserves mention. The eigenvalues $1/b^k$ of the Diaconis-Fulman carry chain are the same as those of the Gilbert-Shannon-Reeds riffle shuffle [1, §4]. The Bernoulli smoothing operator $\mathcal{A}$ used in our proof of the $m$-Bit Equidistribution Lemma (§3.2) is structurally analogous to the descent-smoothing operator that appears in the shuffle analysis: both average over a uniform-mod-$b$ shift and reduce the degree of oscillatory Fourier modes by one per application. For multiplication, the key difference is that only $m = 2$ such smoothings are available at each position, limiting the exact eigenvalues to $\lbrace{}1, 1/b, 1/b^2\rbrace{}$.

### 1.3 Structure of the Paper

The paper develops two complementary spectral viewpoints on the same algebraic object.

- **Part I** (§2–4): *Local dynamics.* The carry transition operator $T_j$ at each position $j$ is a finite Markov chain. We determine its exact eigenvalues ($m$-Bit Equidistribution Lemma) and mixing time.
- **Part II** (§5–6): *Boundary structure.* The carry sequence near the top of the product exhibits a boundary layer whose coefficients $\alpha_k$ encode the per-factor correction to the Euler product. The main proved result is an exact algebraic identity (Proposition 5) linking $\alpha_k$ to the companion matrix spectrum. The empirical observation driving this section — the $(b-1)/b$ anti-correlation law (Conjecture 4) — remains open.
- **Part III** (§7–8): *Global algebra.* The carry polynomial $C(x) = g(x)h(x) - f(x)$ and its quotient $Q(x) = C(x)/(x-b)$ encode the *entire* carry chain as a single algebraic object. The companion matrix $M$ of $Q(x)$ has the carry polynomial roots as eigenvalues. The local mixing rates from Part I constrain how quickly the carry sequence reaches its stationary regime, which in turn controls the root distribution of $Q(x)$. We prove $r_{\max} \leq 3$ for all-positive carry profiles (Corollary 9) via the Eneström-Kakeya theorem, and present evidence for the sharp bound $r_{\max} \leq b$ (Conjecture 10).

### Notation

$N = pq$ (semiprime), $g(x), h(x), f(x)$ (digit polynomials), $C(x) = g(x)h(x) - f(x)$ (carry polynomial), $Q(x) = C(x)/(x-b)$ (quotient), $M$ (companion matrix of $Q$), $D = \deg(Q)$.

---

# Part I: The Diaconis-Fulman Spectrum

## 2. The Carry Transition Operator

### 2.1 Setup and Finite State Space

Let $P(v)$ be any distribution on $\lbrace{}0, \ldots, v_{\max}\rbrace{}$ with $v_{\max} \geq b-1$. The carry operator on functions $\varphi: \mathbb{Z}_{\geq 0} \to \mathbb{R}$:

$$(T\varphi)(c) = \sum_v P(v) \cdot \varphi\!\left(\left\lfloor \frac{v+c}{b}\right\rfloor\right)$$

**Finite state space.** For bounded input ($V \leq v_{\max}$), the carry state space is finite. If $c' = \lfloor(V+c)/b\rfloor$ and $c \geq v_{\max}/(b-1)$, then $c' < c$: the operator is eventually contracting. The invariant state space is $\lbrace{}0, \ldots, N_c\rbrace{}$ with $N_c = \lfloor v_{\max}/(b-1) \rfloor$, and $T$ is an $(N_c+1) \times (N_c+1)$ transition matrix.  (We write $N_c$ for the maximum carry state to avoid collision with $N = pq$ used elsewhere for the semiprime.) All spectral results in this paper concern this finite-dimensional operator.

### 2.2 The Averaged Operator and Gershgorin Bound

**Theorem 1** (Gershgorin Spectral Bound). Let $T$ be the carry transition matrix on $\lbrace{}0, \ldots, N_c\rbrace{}$. Write $v + c = bq + r$ with $r = (v+c) \bmod b$. Then:

$$(Tc^k)(c) = \frac{1}{b^k}\sum_{j=0}^k \binom{k}{j} c^j \cdot \mu_{k-j}(c \bmod b)$$

where $\mu_{k-j}(c_0) = \mathbb{E}[(V - r)^{k-j} \mid c \equiv c_0 \pmod{b}]$ depends on $c$ only through $c \bmod b$. Define the averaged operator $\bar{T}$ by replacing $\mu_{k-j}(c_0)$ with its mean $\bar{\mu}_{k-j}$ over $c_0 \in \lbrace{}0, \ldots, b-1\rbrace{}$. Then $\bar{T}$ is exactly upper-triangular in the monomial basis with diagonal entries $1/b^k$. The perturbation satisfies $\|T - \bar{T}\|_{\ell^1} \leq C \cdot \rho$ where $\rho = \max_{m \neq 0}|\hat{P}(m)| < 1$ is the maximal non-trivial Fourier coefficient of $P$ modulo $b$. For each $k$ such that $\rho < (b-1)/(2b^{k+1})$, the Gershgorin disc around $1/b^k$ contains exactly one eigenvalue of $T$.

*Proof.* We derive the perturbation bound via the Fourier structure of $P$ modulo $b$.

**Claim** (Conditional Moment Bound). For all $s \geq 0$ and $c_0 \in \lbrace{}0, \ldots, b-1\rbrace{}$: $|\mu_s(c_0) - \bar{\mu}_s| \leq v_{\max}^s \cdot (b-1)\rho/b$.

*Proof of Claim.* For fixed $c_0$, write $V + c_0 = bQ + R$ with $R \in \lbrace{}0, \ldots, b-1\rbrace{}$, so $\mu_s(c_0) = \mathbb{E}[(bQ - c_0)^s]$. The residue $R \equiv V + c_0 \pmod{b}$ has distribution $P(R = r) = P(V \equiv r - c_0 \pmod{b})$. Expanding via the DFT of the residue distribution of $P$, with $\hat{P}_b(m) = \mathbb{E}[\omega^{mV}]$ and $\omega = e^{2\pi i/b}$:

$$P(R = r) = \frac{1}{b}\sum_{m=0}^{b-1}\hat{P}_b(m)\,\omega^{m(r - c_0)} = \frac{1}{b} + \frac{1}{b}\sum_{m=1}^{b-1}\hat{P}_b(m)\,\omega^{m(r - c_0)}$$

When $\rho = 0$ (uniform residues), $R$ is uniform regardless of $c_0$, making $\mu_s$ constant. For $\rho > 0$, each non-trivial Fourier mode ($m \neq 0$) perturbs the residue distribution by at most $|\hat{P}_b(m)|/b \leq \rho/b$ per residue class. The moment $\mu_s(c_0) = \mathbb{E}[(bQ - c_0)^s]$ is a linear functional of the distribution of $R$; since $|(bQ - c_0)^s| \leq v_{\max}^s$, the total variation bound

$$\|P_R - \mathrm{Unif}\|_1 \leq (b-1)\rho/b$$

gives $|\mu_s(c_0) - \bar{\mu}_s| \leq v_{\max}^s \cdot (b-1)\rho/b$. $\square$

The constant $C$ admits the explicit bound $C \leq (N_c+1) \cdot v_{\max}^{N_c}$: each row of $T - \bar{T}$ has at most $N_c+1$ nonzero entries, and each entry $(T - \bar{T})_{i,k}$ is a conditional moment difference $|\mu_{k-i}(c_0) - \bar{\mu}_{k-i}| \leq v_{\max}^{k-i} \cdot (b-1)\rho/b$ (Claim above). Summing over $i$ gives the row bound. In practice this bound is extremely loose: since $T$ and $\bar{T}$ are both stochastic matrices with row sums equal to 1, $\|T - \bar{T}\|_{\ell^1} \leq 2$ a priori, and systematic computation across families of non-uniform distributions (biased binomial, truncated geometric, skewed uniform) yields $C_{\text{actual}} \leq 2$ with slack factors exceeding $10^6$ for moderate $N_c$ (A20). $\square$

*Remark (isolation condition).* The condition $\rho < (b-1)/(2b^{k+1})$ ensures that the Gershgorin disc of radius $C\rho$ around the diagonal entry $1/b^k$ does not overlap the disc around $1/b^{k+1}$. More precisely, the disc radius in the monomial basis includes both the perturbation $\|T - \bar{T}\|$ and the off-diagonal entries of $\bar{T}$ itself. Since $\bar{T}$ is upper-triangular with entries decaying as $O(1/b^j)$, the off-diagonal contribution is bounded by $\sum_{j > k} |\bar{T}_{k,j}| \leq C'/b^k$ for a constant $C'$ depending on the moments of $P$. The stated condition absorbs this into the bound on $\rho$ and is sufficient but not sharp.

In particular, when $\rho$ is small (smooth input distribution), the eigenvalues of $T$ are close to $\lbrace{}1, 1/b, \ldots, 1/b^{N_c}\rbrace{}$. For the Diaconis-Fulman addition case ($V \sim \text{Binom}(n, 1/b)$), the m-Bit Equidistribution Lemma (Theorem 2) gives the exact spectrum without perturbation bounds.

---

## 3. The $m$-Bit Equidistribution Lemma

This section contains the paper's main result: a structural theorem that yields exact eigenvalues for carry operators whose input distribution has independent uniform components.

### 3.1 The Parity Lemma

**Proposition 1** (Parity Lemma for Multiplication). For the carry step $c' = \lfloor (V+c)/b \rfloor$ in base $b=2$, the eigenvalue $\lambda_1 = 1/2$ holds for ANY convolution distribution $P(V)$ provided $P(V \text{ odd}) = 1/2$. For multiplication, $\text{conv}_j$ always contains the independent term $g_j \cdot h_0 = g_j \sim \text{Bernoulli}(1/2)$, ensuring uniform parity.

*Proof.* This is the $m = 1$, $k = 1$ case of Theorem 2. Decompose $V = U + S$ where $U \sim \text{Bernoulli}(1/2)$ is independent of $S$ (so $V$ has one uniform-mod-2 component). By the Bernoulli smoothing argument (Step 2 of Theorem 2): $(Tc)(c) = \mathbb{E}_S[\mathbb{E}_U[\lfloor(S + c + U)/2\rfloor]]$. For fixed $w = S + c$: $\frac{1}{2}[\lfloor w/2 \rfloor + \lfloor(w+1)/2 \rfloor] = w/2$. Thus $(Tc)(c) = \mathbb{E}_S[(S+c)/2] = c/2 + \mathbb{E}[S]/2$, which is a degree-1 polynomial in $c$ with leading coefficient $1/2$. Since $T$ maps degree-1 polynomials to degree-1 polynomials with this leading coefficient, the operator is upper-triangular in the monomial basis $\lbrace{}1, c\rbrace{}$ with diagonal entries $\lbrace{}1, 1/2\rbrace{}$, giving eigenvalue $\lambda_1 = 1/2$. $\square$

### 3.2 Statement and Proof

**Theorem 2** ($m$-Bit Equidistribution Lemma). Let $b \geq 2$ and $V = U_1 + \cdots + U_m + S$ where $U_i$ are independent random variables each uniform modulo $b$ (i.e., $P(U_i \equiv r \pmod{b}) = 1/b$ for all $r$), independent of $S \geq 0$. Then the carry operator $(T\varphi)(c) = \mathbb{E}_V[\varphi(\lfloor(V+c)/b\rfloor)]$ satisfies: $(Tc^k)(c)$ is a polynomial of degree $k$ in $c$, with leading coefficient $1/b^k$ and no residue-class oscillation, for all $k \leq m$. Consequently, $T$ is exactly upper-triangular in the monomial basis $\lbrace{}1, c, c^2, \ldots, c^m\rbrace{}$ with diagonal entries $\lbrace{}1, 1/b, 1/b^2, \ldots, 1/b^m\rbrace{}$, and these are the first $m+1$ eigenvalues of $T$ on any state space of dimension $\geq m+1$.

Proof (base $b = 2$; general $b$ below). Write $V = R + S$ with $R = U_1 + \cdots + U_m \sim \text{Binom}(m, 1/2)$, independent of $S$. Then:

$$(Tc^k)(c) = \mathbb{E}_S\!\left[\mathbb{E}_R\!\left[\left\lfloor\frac{S + c + R}{2}\right\rfloor^{\!k}\right]\right]$$

Define $H_k^{(m)}(u) = \mathbb{E}_R[\lfloor(u+R)/2\rfloor^k]$ for integer $u \geq 0$. We show $H_k^{(m)}(u)$ is a polynomial in $u$ (independent of $u \bmod 2$) whenever $m \geq k$.

**Step 1 (Initial oscillation).** For $m = 0$: $H_k^{(0)}(u) = \lfloor u/2 \rfloor^k$. For even $u$: $(u/2)^k$; for odd $u$: $((u-1)/2)^k$. In the parity decomposition $H_k^{(0)}(u) = A_k(u) + (-1)^u B_k(u)$, the smooth part $A_k$ has degree $k$ with leading coefficient $1/2^k$, and the oscillatory part $B_k$ has degree $k-1$.

**Step 2 (Bernoulli smoothing).** The averaging operator

$$\mathcal{A}\varphi(u) = \frac{1}{2}[\varphi(u) + \varphi(u+1)]$$

acts on the parity decomposition as

$$\mathcal{A}[A(u) + (-1)^u B(u)] = A'(u) + (-1)^u B'(u), \qquad B'(u) = -\frac{1}{2}\Delta B(u).$$

Since the forward difference $\Delta$ reduces the degree of a polynomial by 1, $\deg(B') = \deg(B) - 1$.

**Step 3 (Induction).** After $m$ applications of $\mathcal{A}$: $\deg(B^{(m)}) = \max(k - 1 - m, -\infty)$. For $m \geq k$: $B^{(m)} \equiv 0$, so $H_k^{(m)}(u)$ is a pure polynomial in $u$ with leading term $u^k/2^k$.

Since $H_k^{(m)}$ is a polynomial for $k \leq m$: $(Tc^k)(c) = \mathbb{E}_S[H_k^{(m)}(S+c)] = c^k/2^k + (\text{lower degree})$. This establishes exact upper-triangularity with diagonal $1/2^k$ in the monomial basis. The eigenvalues of an upper-triangular matrix are its diagonal entries. $\square$

Generalization to base $b$ (proof sketch; see Open Problem 3 for the formal completion). For general $b$, replace the parity decomposition with the discrete Fourier decomposition modulo $b$: write $H_k^{(0)}(u) = \lfloor u/b \rfloor^k = A_k(u) + \sum_{m=1}^{b-1} \omega^{mu} B_k^{(m)}(u)$, where $\omega = e^{2\pi i/b}$ and each $B_k^{(m)}$ is a polynomial of degree $\leq k-1$. The averaging operator for a uniform-mod-$b$ component is

$$\mathcal{A}\varphi(u) = \frac{1}{b}\sum_{j=0}^{b-1} \varphi(u+j).$$

On the $m$-th Fourier mode ($m \neq 0$):

$$\frac{1}{b}\sum_{j=0}^{b-1} \omega^{m(u+j)} B(u+j) = \omega^{mu} \cdot \frac{1}{b}\sum_j \omega^{mj} B(u+j).$$

Since $\frac{1}{b}\sum_j \omega^{mj} = 0$ for $m \neq 0$, the leading term of $B$ vanishes and only the forward-difference terms survive, reducing $\deg(B)$ by 1. The induction proceeds identically: after $m$ uniform-mod-$b$ averagings, all oscillatory modes of degree $\leq k-1-m$ vanish. For $m \geq k$: no oscillation remains, and the leading coefficient is $1/b^k$. $\square$

**Corollary 2.1** (Diaconis-Fulman Spectrum for Addition). For base-$b$ addition of $n$ independent $\text{Uniform}\lbrace{}0, \ldots, b-1\rbrace{}$ digits, $V$ has $n$ independent uniform-mod-$b$ components and $v_{\max} = n(b-1)$, so the full finite operator acts on $\lbrace{}0, \ldots, n\rbrace{}$. Starting from $c_0 = 0$, however, the reachable carries are $\lbrace{}0, \ldots, n-1\rbrace{}$. Restricting to this reachable subspace recovers the classical Diaconis-Fulman spectrum $\lbrace{}1, 1/b, \ldots, 1/b^{n-1}\rbrace{}$ [1, Theorem 1].

### 3.3 Corollaries for Multiplication

**Corollary 2.2** (Double Universality). For base-2 multiplication at position $j \geq 2$: $\text{conv}_j = g_j + h_j + S$ where $g_j, h_j \sim \text{Bernoulli}(1/2)$ are independent of $S = \sum_{i=1}^{j-1} g_i h_{j-i}$. Thus $m = 2$, and the transfer operator $T_j$ has eigenvalues $\lbrace{}1, 1/2, 1/4\rbrace{}$ exactly. Verified for $j = 2, \ldots, 10$ with exact rational arithmetic (A13).

**Corollary 2.3** (Composition Spectral Gap). The position-dependent operators $T_1, \ldots, T_K$ are each upper-triangular in the monomial basis through degree 2 with diagonal entries $(1, 1/b, 1/b^2)$ (Corollary 2.2 for $b=2$; and, for general prime $b$, by the proof-sketch extension summarized in Corollary 2.4). The product of upper-triangular matrices is upper-triangular with the product of diagonal entries, so the composed operator $T_1 \cdots T_K$ has eigenvalues $(1, (1/b)^K, (1/b^2)^K)$ on the degree-$\leq 2$ polynomial subspace. The dominant sub-leading eigenvalue is $(1/b)^K$, governing the convergence rate of boundary observables. In particular, the trace anomaly satisfies $\alpha_1(K) = c_1 + P(K) \cdot (1/2)^K$ (base 2), where the limit $c_1$ and the polynomial prefactor $P(K)$ are characterized in [E]; direct enumeration to $K = 21$ yields $c_1 = 0.17453\ldots$, matching $\pi/18$ to 4.3 digits (conjectured but not proved).

**Corollary 2.4** (Base-$b$ Generalization; proof-sketch status). For prime base $b$, the multiplication convolution at position $j \geq 2$ contains the terms $g_j \cdot h_0$ and $g_0 \cdot h_j$, where $g_j, h_j \sim \text{Uniform}\lbrace{}0, \ldots, b-1\rbrace{}$ and $g_0, h_0 \in \lbrace{}1, \ldots, b-1\rbrace{}$ are coprime to $b$. These provide two independent uniform-mod-$b$ additive components, and the proof sketch above yields the same first three eigenvalues $\lbrace{}1, 1/b, 1/b^2\rbrace{}$. Verified exactly for bases 2, 3, 5, 7 (A13).

### 3.4 Exponential Convergence of Higher Eigenvalues

**Proposition 3** (Exponential Convergence). For prime base $b$, the $k$-th eigenvalue of the multiplication transfer operator $T_j$ satisfies:

$$\lambda_k^{(j)} = \frac{1}{b^k} + O\!\left(j^{k-3} \cdot b^{-(j-1)}\right) \quad \text{for all } k \geq 3, \; j \geq 2, \; \dim(T_j) > k$$

In particular, $\lambda_k^{(j)} \to 1/b^k$ exponentially fast as $j \to \infty$.  (The condition $\dim(T_j) > k$ is needed so that the $k$-th eigenvalue exists; since $\dim(T_j) = \lfloor v_{\max}(j)/(b-1) \rfloor + 1$ grows linearly with $j$, it is satisfied for all $j \geq k$ in practice.)

*Proof.* The convolution $\text{conv}_j = g_j h_0 + g_0 h_j + \sum_{i=1}^{j-1} g_i h_{j-i}$ decomposes into $j+1$ mutually independent terms (each uses a unique pair of digit bits). The first two terms provide two uniform-mod-$b$ components (Corollary 2.4). The remaining $n = j-1$ terms $X_i = g_i h_{j-i}$ are products of independent uniform digits.

**Damping factor.** For $g, h \sim \text{Uniform}\lbrace{}0, \ldots, b-1\rbrace{}$ independently and $b$ prime, the characteristic function of $X = gh$ at frequency $m/b$ ($m \not\equiv 0$) is:

$$\mathbb{E}[\omega^{mX}] = \frac{1}{b^2}\sum_{a=0}^{b-1}\sum_{c=0}^{b-1} \omega^{mac} = \frac{1}{b}$$

where $\omega = e^{2\pi i/b}$. (The inner sum over $c$ is $b$ when $a = 0$ and $0$ otherwise, since $ma \not\equiv 0 \pmod{b}$ for $a \neq 0$ when $b$ is prime.)

**Iterated damping.** The averaging operator $\mathcal{A}_{X_i}$ acts on the parity decomposition $\varphi(c) = A(c) + (-1)^c B(c)$ (base 2) as:

$$\mathcal{A}_{X_i}[(-1)^c B(c)] = (-1)^c\!\left[(1-2p)B(c) - p\,\Delta B(c)\right]$$

where $p = P(X_i = 1)$. For $X_i \sim \text{Bernoulli}(1/4)$: $|1 - 2p| = 1/2$. This damps the leading coefficient of $B$ by $1/2$ per application, while $\Delta$ reduces subleading terms by one degree. Iterating: after $n$ applications, the leading coefficient of $B^{(n)}$ is multiplied by $(1/2)^n$, and subleading terms satisfy $|b_l^{(n)}| \leq C \cdot n^{k-3-l} \cdot (1/2)^n$ for each degree $l$. To see this, represent $B^{(n)}(c) = \sum_{l=0}^{k-3} b_l^{(n)} c^l$ by its coefficient vector $\boldsymbol{b}^{(n)} \in \mathbb{R}^{k-2}$. Each smoothing step acts as $\boldsymbol{b}^{(n+1)} = L \cdot \boldsymbol{b}^{(n)}$, where $L$ is $(k-2) \times (k-2)$ lower-triangular: the diagonal entries are $L_{ll} = (1 - 2p) = 1/2$ (the damping factor on each coefficient's leading contribution), and the subdiagonal entries arise from the forward difference $\Delta$ mixing degree $l+1$ into degree $l$. Since $L$ is lower-triangular with constant diagonal $1/2$, its eigenvalues are $1/2$ with multiplicity $k-2$, and the Jordan form gives the polynomial prefactor $n^{k-3-l}$.

For general prime base $b$, the parity decomposition is replaced by the DFT modulo $b$: write $\varphi(c) = A(c) + \sum_{m=1}^{b-1}\omega^{mc}B^{(m)}(c)$. The averaging operator $\mathcal{A}_{X_i}$ acts on each mode $m \neq 0$ as

$$\mathcal{A}_{X_i}[\omega^{mc}B^{(m)}(c)] = \omega^{mc}[\mathbb{E}[\omega^{mX_i}]\,B^{(m)}(c) + \text{lower-order terms from } \Delta B^{(m)}].$$

Since $|\mathbb{E}[\omega^{mX_i}]| = 1/b$ (damping factor, proved above), the leading coefficient of each $B^{(m)}$ is damped by $1/b$ per application, and the iteration matrix $L$ has diagonal entries $1/b$ instead of $1/2$. The Jordan form argument is identical.

**Combining.** The two uniform-mod-$b$ components (Theorem 2) eliminate oscillatory degrees 0 through 1 for $k \leq 2$. For $k \geq 3$, after the two uniform smoothings the residual oscillation has degree $k-3$. The $n = j-1$ product terms then damp the $l$-th coefficient to $|b_l^{(n)}| \leq C \cdot n^{k-3-l} \cdot b^{-n}$. On the finite state space $\lbrace{}0, \ldots, N_j\rbrace{}$ with $N_j = O(j)$, expressing the oscillatory function in the monomial basis and summing the Gershgorin off-diagonal entries in row $k$ yields:

$$|\lambda_k^{(j)} - 1/b^k| \leq \sum_{l=0}^{k-3} |b_l^{(j-1)}| \cdot N_j^l \leq C' \cdot j^{k-3} \cdot b^{-(j-1)} \cdot \sum_{l=0}^{k-3} \left(\frac{j}{j-1}\right)^l = O(j^{k-3} \cdot b^{-(j-1)})$$

where the geometric sum $\sum_{l=0}^{k-3}(j/(j-1))^l$ is bounded by $k-2$ for all $j \geq 2$ (since $j/(j-1) \leq 2$ and $l \leq k-4$), contributing only a constant factor absorbed into $O(\cdot)$. $\square$

*Empirical verification (A13, A19):* $\lambda_3^{(3)} = 3/32 = 0.09375$ (error $0.031$ from $1/8$); $\lambda_3^{(5)} \approx 0.1195$ (error $0.006$); $\lambda_3^{(8)} \approx 0.1250$ (error $< 0.001$). The geometric rate $\sim (1/2)^j$ is confirmed. Fitting $|\lambda_k^{(j)} - 1/b^k| \sim C \cdot j^\alpha \cdot (1/b)^j$ for $k = 3, 4, 5$ yields $\alpha \approx k - 3$ (A19), confirming the tight exponent.

**Observation** ([G, experiment G08]). The base-3 multiplication transfer operator $T_j$ has universal eigenvalues $\lbrace{}1, 1/3, 1/9\rbrace{}$ for all $j = 1, \ldots, 6$, with $\lambda_3^{(j)} \to 1/27$ as $j \to \infty$.

---

## 4. Total Variation Mixing

**Theorem 3** (Mixing Bound). Let the convolution input $V$ satisfy the uniform-mod-$b$ property (Theorem 2 with $\rho = 0$), so that the remainder $r_t = (V_t + c_t) \bmod b$ is uniform on $\lbrace{}0, \ldots, b-1\rbrace{}$ at each step. Then the carry chain mixes in $O(\log(c_{\max})/\log b)$ steps:

$$\|P_t - \pi\|_{TV} \leq c_{\max}/b^t$$

*Proof via coupling.* Consider two copies of the carry chain started from arbitrary states $c_0$ and $c'_0$ with $c'_0 > c_0 \geq 0$, and couple them by feeding the *same* random convolution digit $v_t$ to both at each step. The carry update is $c_{t+1} = \lfloor(v_t + c_t)/b\rfloor$. For the gap $d_t = c'_t - c_t \geq 0$, write $r_t = (v_t + c_t) \bmod b$. Then:

$$d_{t+1} = \left\lfloor\frac{r_t + d_t}{b}\right\rfloor$$

since the integer parts cancel. For the multiplication carry chain, $r_t$ is uniform on $\lbrace{}0, \ldots, b-1\rbrace{}$ (this follows from the uniform-mod-$b$ property of the convolution, Theorem 2 with $\rho = 0$). Using the identity $\sum_{r=0}^{b-1}\lfloor(r+d)/b\rfloor = d$ for any non-negative integer $d$:

$$\mathbb{E}[d_{t+1} \mid \mathcal{F}_t] = \frac{1}{b}\sum_{r=0}^{b-1}\left\lfloor\frac{r + d_t}{b}\right\rfloor = \frac{d_t}{b}$$

where

$$\mathcal{F}_t = \sigma(v_0, \ldots, v_{t-1})$$

is the natural filtration. By the tower property: $\mathbb{E}[d_t] \leq d_0/b^t \leq c_{\max}/b^t$. By Markov's inequality, $\Pr(d_t \geq 1) \leq \mathbb{E}[d_t] \leq c_{\max}/b^t$. Once $d_t = 0$, the two chains occupy the same state and, receiving identical inputs, remain coupled for all subsequent steps. By the coupling inequality:

$$\|P_t - \pi\|_{TV} \leq \Pr(\text{not coupled by step } t) \leq \Pr(d_t \geq 1) \leq c_{\max}/b^t$$

This gives exact coalescence with probability 1 in finite time (since $\sum_t c_{\max}/b^t < \infty$), but the coalescence is *probabilistic*, not deterministic: the ceiling bound $d_{t+1} \leq \lceil d_t/b \rceil$ only guarantees $d_t \leq 1$ after $\lceil\log_b c_{\max}\rceil$ steps, and the final reduction from $d = 1$ to $d = 0$ requires $(v_t + c_t) \bmod b \neq b - 1$, which occurs with probability $(b-1)/b$ at each step. $\square$

*Remark.* For general input distributions with $\rho = \max_{m \neq 0}|\hat{P}(m)| > 0$, the expected gap satisfies $\mathbb{E}[d_{t+1} \mid d_t] = d_t/b + O(d_t \rho/b)$, and the TV bound becomes $\|P_t - \pi\|_{TV} \leq C(\rho) \cdot c_{\max}/b^t$ with a constant depending on $\rho$.

### 4.1 Computational Verification

| Base | $n$ | $\lambda_1$ | $\lambda_2$ | $\lambda_3$ | $\lambda_4$ |
|------|-----|-------------|-------------|-------------|-------------|
| 2 | 32 | 0.500000 | 0.250000 | 0.125000 | 0.062500 |
| 2 | 128 | 0.500000 | 0.250000 | 0.125000 | 0.062500 |
| 3 | 64 | 0.333333 | 0.111111 | 0.037037 | 0.012346 |
| 5 | 64 | 0.200000 | 0.040000 | 0.008000 | 0.001600 |
| 7 | 32 | 0.142857 | 0.020408 | 0.002915 | 0.000416 |

Agreement to 15+ digits for $n \geq 32$, all bases tested.

### 4.2 Convergence of Boundary Observables

The Diaconis-Fulman rate $\lambda_1 = 1/b$ governs *local* carry-value mixing. The trace anomaly $\alpha_1 := \mathbb{E}[c_{M-1}] - 1$ converges to $\pi/18$ (base 2) at rate $\rho(K) \to 1/2$ from above [E]. This rate is precisely the sub-dominant eigenvalue $1/b$, verified for bases 2, 3, 5, 7, 10 [experiments A14–A18]. Per-cascade-depth convergence rates satisfy $\rho_J \to 1$ for all $J \geq 2$: carry values mix exponentially, but per-cascade-depth boundary observables (which depend on bit-carry covariances from shared bits between adjacent convolutions) decay sub-exponentially.

**Connection to Fulman 2023 [8].** Fulman's property $K_a \cdot K_b = K_{ab}$ (multiplicative composition of carries transition matrices) constrains the spectral structure: the carries chain in base $b^k$ factors through $k$ applications of $K_b$, each with spectral gap $1-1/b$. For our multiplication carries, the composition $T_1 \cdots T_K$ inherits the dominant eigenvalue $(1/2)^K$ (Corollary 2.3).

---

# Part II: Anti-Correlation and the Correction Series

## 5. The $(b-1)/b$ Anti-Correlation Law

**Conjecture 4** (Anti-Correlation Law). For the multiplication carry chain in base $b$, the correction coefficients $\alpha_k = \langle c_{\text{top}-k}\rangle - \langle c_{\text{top}-k+1}\rangle$ satisfy:
1. Stationary mean: $\alpha_k \to (b-1)/(2b)$ for $k \geq 5$ (base 2).
2. One-step anti-correlation: the correction to $\alpha_k$ relative to the bulk value $(b-1)/(2b)$ decays with effective rate $(b-1)/b$, driven by the spectral gap $1 - 1/b$.

Part 1 (stationary mean $(b-1)/(2b)$) is verified across 3.17 million measurements in 18 bases; Part 2 (decay rate $(b-1)/b$) is confirmed with all bases within $1\sigma$ of the predicted rate; Z-score 32.4 against random control; Z-score > 127 across scaling tests 20–512 bit. A formal proof of both parts remains an open problem; the natural approach is a boundary transfer theorem extending the covariance induction of [F].

**Clarification (A12).** The *raw* carry correlation $\mathrm{Corr}(c_j, c_{j+1})$ at interior positions is **strongly positive** ($\approx 0.85$–$0.88$ for base 2, increasing with $K$), because the mean $\mathbb{E}[c_j] = (j-1)/4$ grows linearly and consecutive carries are co-monotone. The "anti-correlation" refers to the **top-boundary regime**: the coefficients $\alpha_k$ overshoot $(b-1)/(2b)$ and converge oscillatorily, with the oscillation damping at rate $(b-1)/b$.

*Partial analytical support.* The off-diagonal covariance

$$\mathrm{Cov}(c_j, g_i h_{j-i}) = 1/8$$

(universal, proved in [F, Lemma D]) gives

$$\mathrm{Cov}(c_j, \mathrm{conv}_j) = (j-1)/8$$

for odd $j$. Combined with the carry recurrence, the positive covariance contributes to the boundary layer structure. The spectral gap $1/b$ quantifies the decorrelation rate.

## 6. The Correction Series

**Proposition 5** (Algebraic Identity). For test prime $l$:

$$h(l) = \langle|\det(I - M/l)|\rangle \cdot (1 - 1/l) = 1 + \sum_{k=1}^{D-1} \frac{\alpha_k}{l^k}$$

where $\alpha_k = \langle c_{\text{top}-k}\rangle - \langle c_{\text{top}-k+1}\rangle$.

*Proof.* This is an algebraic consequence of the Carry Representation Theorem (CRT): $|\det(I - M/l)| = |C(l)|/|l-b|^D$ where $C(l) = g(l)h(l) - n(l)$, and $C(l)/(l-b) = Q(l)$ is the carry quotient polynomial with coefficients $c_{\text{top}-k}$. The identity follows by expanding $|Q(l)|/l^D$ as a power series in $1/l$. (Full derivation of the CRT in [B, §2]; the expansion step is elementary.) Verified to $\Delta < 10^{-15}$ across 155,000 semiprimes.

### 6.1 Boundary Coefficients

| Bit size | $\alpha_1$ | $\alpha_2$ | $\alpha_3$ | $\alpha_k$ ($k = 5, \ldots, 15$) |
|----------|------------|------------|------------|-----------------------------------|
| 16-bit | 0.164 | 0.188 | 0.234 | $0.250 \pm 0.005$ |
| 32-bit | 0.171 | 0.191 | 0.239 | $0.250 \pm 0.003$ |
| 64-bit | 0.178 | 0.187 | 0.240 | $0.250 \pm 0.002$ |

For $k \geq 5$: $\alpha_k \to 1/4$ (base 2), reflecting steady-state carry increment. The 2-position boundary layer is explained by the mixing time $t_{\text{mix}} \leq 2$ (Theorem 3) and the ULC constraint $c_{\text{top}} = 1$.

**High-precision verification (A02, A03).** Monte Carlo ($d = 48$, 500K samples) confirms:
- $\alpha_1 \approx 0.1746$ (consistent with $\pi/18 = 0.17453\ldots$; see [E] for the $\pi/18$ conjecture — 4.3 digits from direct enumeration at $K = 21$)
- $\alpha_5 \approx 0.257$ (overshoots $1/4$, oscillatory approach)
- $S_{\text{boundary}} = \sum_{k=1}^{\infty}(\alpha_k - 1/4) \approx -0.148$ (converged, no closed form found)

**Proved covariance structure (from [F]).** In the bulk (without MSB conditioning),

$$\mathbb{E}[\mathrm{carry}_j] = (j-1)/4 \quad \text{for } j \geq 1,$$

$$\mathrm{Cov}(\mathrm{carry}_j, \mathrm{conv}_j) = (j-1)/8 \quad \text{for odd } j \geq 3 \; (\text{verified to } j = 11),$$

and the induction step is

$$\mathrm{Cov}(c_{j+2}, v_{j+2}) - \mathrm{Cov}(c_j, v_j) = 1/4 \quad \text{for odd } j \geq 3.$$

For even $j$, the covariance exceeds $(j-1)/8$ by a correction $\epsilon_j > 0$ decaying exponentially: $\epsilon_2 = 1/16$, $\epsilon_4 = 1/32$, $\epsilon_6 = 1/128$, $\epsilon_8 = 1/256$, $\epsilon_{10} = 3/4096$, $\epsilon_{12} = 11/32768$.

**Remark (Towards a proof of Conjecture 4).** The bulk covariance results from [F] establish the interior carry statistics unconditionally. The missing ingredient is a *boundary transfer theorem*: showing that conditioning on $c_{\text{top}} = 1$ (the ULC constraint, proved in [B]) perturbs the interior covariance by exactly the oscillatory correction $\alpha_k - (b-1)/(2b) \sim (-1)^k \cdot ((b-1)/b)^k$. The spectral gap $1 - 1/b$ from Theorem 2 governs the exponential decorrelation rate; what remains is to compute the boundary amplitude. This requires extending the covariance induction of [F] from the unconditional bulk to the MSB-conditioned regime — a tractable but technically demanding computation that we leave to future work. Monte Carlo experiments at 24-bit and 32-bit confirm the interior mean converges to $(b-1)/(2b) = 0.250 \pm 0.003$ with boundary decay ratios consistent with $-(b-1)/b$ (A23).

---

# Part III: Spectral Radius Bounds

Parts I–II studied carries as a *stochastic process*: position-by-position dynamics, eigenvalues of the transfer operator, and the boundary layer near the top of the product. We now shift to the *algebraic* perspective. The Carry Representation Theorem (CRT, proved in [B]) establishes that $g(x)h(x) = f(x) + (x-b)Q(x)$, where the carry quotient polynomial $Q(x)$ has coefficients $q_k = c_{\text{top}-k}$ — the very same carry values studied in Part II. The companion matrix $M$ of $Q(x)$ has characteristic polynomial $\det(zI - M) = z^D Q(1/z)$, so its eigenvalues are the reciprocals of the roots of $Q$. The spectral radius $r_{\max} = \max|z_i|$ over the roots of $Q$ governs the convergence of the per-factor Euler correction series (Proposition 5): $r_{\max} < b$ is equivalent to absolute convergence of $\sum \alpha_k / l^k$ for all primes $l \geq b$. Thus the algebraic question "how large can the roots of $Q$ be?" is the same as the analytical question "does the correction series converge?"

**Concrete example.** $N = 49 = 7 \times 7$ has carry profile $c = (0, 1, 2, 2, 1)$ and $Q(x) = x + 2x^2 + 2x^3 + x^4$. The companion matrix of $Q$ has spectral radius $r_{\max} \approx 1.31 < b = 2$: the carry values studied as stochastic observables in Part II become the polynomial coefficients whose root distribution Part III bounds. The mixing results of Part I (Theorem 3, Proposition 3) constrain how fast the carries reach their stationary distribution, which in turn limits the growth of the $c_k$ and hence of $r_{\max}$.

## 7. The Carry Companion Matrix

### 7.1 Structural Lemmas

**Lemma 6** (Backward Recursion). $c_{\text{top}-1} = b \cdot c_{\text{top}} + f_{m-1} - \text{conv}_{m-1}$.

*Proof.* The carry recurrence at position $m$ reads

$$c_m = \lfloor(\mathrm{conv}_m + c_{m-1})/b\rfloor,$$

so

$$\mathrm{conv}_m + c_{m-1} = b \cdot c_m + f_m,$$

where $f_m$ is the $m$-th digit of the product. Setting $m = m{-}1$ (the penultimate position) and rearranging:

$$c_{\text{top}-1} = b \cdot c_{\text{top}} + f_{m-1} - \mathrm{conv}_{m-1}.$$

$\square$ (Verified with zero exceptions across 15,000 tests.)

**Lemma 7** (Tight Upper Bound). For base 2: $c_{\text{top}-1} \leq 3$.

*Proof.* From Lemma 6 with ULC ($c_{\text{top}} = 1$): $c_{\text{top}-1} = 2 + f_{m-1} - \text{conv}_{m-1} \leq 2 + 1 - 0 = 3$. $\square$

**Lemma 8** ($D = 2d$ case). When $D = 2d$ (i.e., both factors have the same number of digits): $c_{\text{top}-1} \leq 2$.

*Proof.* When $D = 2d$: $\text{conv}_{m-1} = g_{d-1}h_{d-1} = 1$ (both MSBs are 1), so $c_{\text{top}-1} = 1 + f_{m-1} \leq 2$. $\square$

The condition $D = 2d$ holds for approximately 61.4% of semiprimes in the range tested (A08). This proportion depends on the bit-length distribution of the factors and is not a structural constant.

### 7.2 The Spectral Radius Bound

**Corollary 9** (Eneström-Kakeya Bound). For base-2 semiprimes with all carry coefficients positive: $r_{\max} \leq 3$.

*Proof.* We use the carry identity $\text{conv}_k + c_k = b \cdot c_{k+1} + f_k$, giving $c_k = b \cdot c_{k+1} + f_k - \text{conv}_k$ for each position $k$. Since $f_k \in \lbrace{}0, 1\rbrace{}$ and $\text{conv}_k \geq 0$:

$$c_k \leq 2 c_{k+1} + 1 \quad \text{for all } k$$

When $c_{k+1} > 0$, this gives $c_k / c_{k+1} \leq 2 + 1/c_{k+1} \leq 3$, with equality only when $c_{k+1} = 1$ (i.e., at the top of the carry profile). The Eneström-Kakeya theorem [2, 3] states that for a polynomial $a_0 + a_1 z + \cdots + a_n z^n$ with all $a_i > 0$, all roots satisfy $|z| \leq \max_k(a_k / a_{k+1})$. Applied to $Q(x) = c_1 + c_2 x + \cdots + c_D x^{D-1}$ when all $c_k > 0$:

$$r_{\max} \leq \max_{k=1}^{D-1}\frac{c_k}{c_{k+1}} \leq 3 \qquad \square$$

*Remark (zero-coefficient profiles).* The Eneström-Kakeya theorem requires all coefficients to be strictly positive. Approximately 7% of semiprimes have $c_{\text{top}-1} = 0$ (§7.3), producing a zero coefficient in $Q$. For these:

- **Leading zeros** ($c_1 = \cdots = c_j = 0$, $c_{j+1} > 0$): $Q(x) = x^j \tilde{Q}(x)$, and the nonzero roots are those of $\tilde{Q}$, to which E-K applies if the remaining coefficients are all positive.
- **Interior zeros** ($c_k = 0$ with $c_{k-1}, c_{k+1} > 0$): E-K does not apply. By the Cauchy bound, $r_{\max} \leq 1 + \max_k c_k / c_{\text{top}}$. Since $c_{\text{top}} = 1$ (ULC), interior carries can be as large as $O(K)$ for $K$-bit factors, so this bound is weak. We rely on computational verification: among all 100,000+ semiprimes tested (A07, A11, A22), no zero-coefficient profile produced $r_{\max} > 2$.

**Proved vs. conjectured in Part III.** The E-K bound $r_{\max} \leq 3$ is proved for all-positive profiles (~93% of semiprimes). For zero-coefficient profiles, $r_{\max} \leq 3$ holds computationally but lacks a uniform analytical proof. The sharper bound $r_{\max} \leq b = 2$ is conjectured (Conjecture 10, §8).

### 7.3 Statistical Distribution

| $c_{\text{top}-1}$ | Frequency | $r_{\max}$ range |
|---------------------|-----------|------------------|
| 0 | ~7% | $[0.5, 1.0]$ |
| 1 | ~76% | $[0.8, 1.5]$ |
| 2 | ~17% | $[1.0, 2.0]$ |
| 3 | ~0.002% | $[1.5, 2.0]$ |

Median $r_{\max} \approx 1.19$, 97% have $r_{\max} < \sqrt{2}$.

## 8. The Tight Bound $r_{\max} \leq b$

**Conjecture 10** (Tight Spectral Radius Bound). For base-$b$ semiprimes: $r_{\max} \leq b$, with equality iff the carry profile is $[0, \ldots, 0, b, 1]$ (degenerate).

Verified for 100,000+ semiprimes across bit sizes 8–24 (A07, A11). The degenerate profile produces $C(x) = x + b$ with root $z = -b$, giving $r_{\max} = b$ exactly. Exhaustive search for $D = 3 \ldots 7$ and gradient ascent for $D \leq 15$ confirm no realizable profile achieves $r_{\max} \geq b$; the minimum Rouché gap $\min_\theta |Q(be^{i\theta})|$ grows exponentially as $\sim 2^{0.79D}$ (A21, A22).

### 8.1 The Rouché Path (A11)

On $|z| = b$, the carry polynomial $C(z) = g(z)h(z) - n(z)$ satisfies $C(b) = 0$ and, experimentally, $C(z) \neq 0$ for $z \neq b$ on $|z| = b$:

| Bits | Semiprimes | $\min_{\theta \neq 0} \lvert C(be^{i\theta}) \rvert$ | All $> 0$? |
|------|------------|------------------------------------------|------------|
| 8 | 476 | 0.39 | **Yes** |
| 12 | 498 | 6.5 | **Yes** |
| 16 | 500 | 493 | **Yes** |
| 20 | 500 | 951,615 | **Yes** |

The minimum grows as $\sim 2^{2K}$, suggesting the bound strengthens with factor size. If $C(z) \neq 0$ on $|z| = b, z \neq b$ could be proved, then $Q(z) = C(z)/(z-b)$ would have all roots inside $|z| < b$ by the argument principle.

### 8.2 Near-Extremal Structure (A11)

The extremal eigenvalue is always **real negative** ($z \approx -r_{\max}$). Near-extremal carry profiles are $[\ldots, 1, 0, 1, 0, 1]$: alternating carries maximizing the negative-real root. The ratio $r_{\max}/b$ decreases with base:

| Base | $\max(r_{\max}/b)$ | $\text{mean}(r_{\max}/b)$ |
|------|-------------------|--------------------------|
| 2 | 0.77 | 0.59 |
| 3 | 0.59 | 0.42 |
| 5 | 0.41 | 0.28 |
| 7 | 0.31 | 0.21 |
| 10 | 0.26 | 0.16 |

### 8.3 Computational Evidence

Extensive computational searches over all semiprimes $N = pq$ with $p, q < 10^4$ revealed no counterexample to $r_{\max} \leq 2$, consistent with Conjecture 10.

---

## 9. Conclusion

We have developed the spectral theory of carries in positional multiplication, extending the Diaconis-Fulman framework from addition to multiplication. The main contribution is the $m$-Bit Equidistribution Lemma (Theorem 2), which in the binary case yields exact eigenvalues $\lbrace{}1, 1/b, \ldots, 1/b^m\rbrace{}$ for carry operators with $m$ independent uniform-mod-$b$ components, via a Bernoulli-smoothing argument on the monomial basis; the general prime-base extension is given in proof-sketch form in §3.2. For multiplication, this gives universal eigenvalues $\lbrace{}1, 1/b, 1/b^2\rbrace{}$ at each position $j \geq 2$, with higher eigenvalues converging exponentially (Proposition 3). The algebraic perspective (Part III) connects the local mixing rates to the global root distribution of the carry quotient polynomial, yielding a proved spectral radius bound $r_{\max} \leq 3$ for all-positive carry profiles via the Eneström-Kakeya theorem. Two conjectures remain open: the $(b-1)/b$ anti-correlation law (Conjecture 4), for which a boundary transfer theorem extending the covariance induction of [F] is the natural path; and the tight bound $r_{\max} \leq b$ (Conjecture 10), for which the Rouché gap analysis provides strong numerical evidence; a formal proof via the Rouché gap remains open.

### Open Problems

1. **Anti-correlation law (Conjecture 4).** Prove that $\alpha_k \to (b-1)/(2b)$ for all $k$ sufficiently large and all bases $b$, with oscillatory corrections decaying at rate $(b-1)/b$ (the effective anti-correlation rate from Conjecture 4, part 2). The partial proof (boundary transfer from [F]) reduces this to showing that the per-position covariance induction extends from the central region to the boundary layer. The convergence rate is at most $(1/b)^k$ (Proposition 3).

2. **Tight spectral radius bound (Conjecture 10).** Prove $r_{\max} \leq b$ for all carry quotient polynomials. The Eneström-Kakeya theorem gives $r_{\max} \leq 3$ for all-positive profiles (Corollary 9, §7.2). Phase dispersion analysis (§7.2) provides strong numerical evidence for $r_{\max} \leq b$ but the Rouché gap at the critical radius remains open.

3. Base-$b$ generalization of Theorem 2. The $m$-Bit Equidistribution Lemma (Theorem 2) is stated and proved for base 2. The generalization to base $b$ (used in Corollary 2.1) follows by replacing Bernoulli(1/2) with Uniform-mod-$b$ in the DFT-based smoothing argument sketched in §3.2. A fully detailed write-up of the base-$b$ induction would complete the formal proof.

---

## 10. Reproducibility

All experiments (23 scripts + 1 shared utility module) are in `experiments/`. Key scripts by topic:

| Topic | Key scripts |
|-------|------------|
| Spectrum | `A05_markov_transition.py`, `A06_diaconis_fulman.py` |
| $m$-Bit Lemma | `A13_equidistribution.py` (Theorem 2, Corollaries 2.1–2.4, Proposition 3) |
| Boundary coefficients | `A02_constant_c_precision.py`, `A03_constant_c_method.py` |
| Anti-correlation | `A04_c2_decomposition.py`, `A01_perfactor_identity.py`, `A12_anticorrelation.py`, `A23_boundary_transfer.py` |
| Spectral radius | `A07_rouche_rmax.py`, `A09_gap3_closure.py`, `A10_counterexample.py`, `A08_boundary_proof.py`, `A11_rouche_bound.py`, `A21_rouche_formal.py`, `A22_adversarial.py` |
| Proposition 3 exponent | `A19_prop3_exponent.py` |
| Theorem 1 constant $C$ | `A20_constant_C.py` |

Requirements: Python 3.8+, NumPy, SciPy.

---

## References

1. P. Diaconis, J. Fulman, "Carries, Shuffling, and Symmetric Functions," *Adv. Appl. Math.* 43(2), 176–196, 2009.
2. G. Eneström, "Härledning af en allmän formel," *Öfversigt af K. Vet.-Akad. Förhandlingar* 50, 405–415, 1893.
3. S. Kakeya, "On the limits of the roots of an algebraic equation," *Tôhoku Math. J.* 2, 140–142, 1912.
4. J. Holte, "Carries, Combinatorics, and an Amazing Matrix," *Amer. Math. Monthly* 104(2), 138–149, 1997.
5. D. E. Knuth, *The Art of Computer Programming, Vol. 2*, 3rd ed., Addison-Wesley, 1997.
6. [E] Companion paper: "The Trace Anomaly of Binary Multiplication," this series. doi:10.5281/zenodo.18895604
7. [F] Companion paper: "Exact Covariance Structure of Binary Carry Chains," this series (`carry-arithmetic-F-covariance-structure` repo). doi:10.5281/zenodo.18895607
8. J. Fulman, "Carries and a Map on the Space of Rational Functions," *arXiv:2306.05529*, 2023.
9. [B] Companion paper: "Carry Polynomials and the Euler Product: An Approximation Framework," this series. (Contains the Unit Leading Carry theorem used in §6.). doi:10.5281/zenodo.18895597

---

*CC BY 4.0. Code: MIT License.*
