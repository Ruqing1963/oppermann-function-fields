# Oppermann intervals in F_q[t]: what is true, what is proved, what is out of reach

Independent project, 2026-10-12. Parts I–IV are not modified. Scripts in `code/`, outputs in `data/`.

## 1. The statement

For monic f of degree d the Legendre interval of Parts I–IV is I_f = {f² + s : deg s ≤ d} (q^{d+1} elements). It splits into the q sub-intervals

  O_{f,c} = {f² + c·f + s : deg s < d} = {f² + c·t^d + s : deg s < d},   c ∈ F_q,

each with q^d elements: the short intervals of degree 2d and length q^d = |f| next to f². For q = 2 these are exactly two halves (c = 0, c = 1), the analogue of Oppermann's (n², n²+n) and (n²+n, (n+1)²). **Function-field Oppermann:** every O_{f,c} contains an irreducible polynomial.

## 2. The statement is false in characteristic 2 (proved)

**Proposition A (d = 1, every even q).** O_{f,0} = {f² + s₀ : s₀ ∈ F_q} consists of squares, (f + √s₀)², so it contains no irreducible. There are q such empty sub-intervals.

**Proposition B (d = 2, q = 4^j).** Let f = t² + a₁t + a₀. Since f² = t⁴ + a₁²t² + a₀², the sub-interval O_{f,c} is the block {t⁴ + (a₁² + c)t² + αt + β : α, β ∈ F_q}. If k = log₂ q is even, the block with a₁² + c = 0, i.e. {t⁴ + αt + β}, contains no irreducible. This gives q² empty Oppermann sub-intervals for every q = 4, 16, 64, ….

*Proof.* Let P = t⁴ + e t² + αt + β with α ≠ 0, roots r₁, …, r₄ (sum 0). For each of the three pairings {ij|kl}, u = r_i + r_j = r_k + r_l satisfies (t² + ut + v)(t² + ut + w) = P, hence u² + v + w = e, u(v + w) = α, and u³ + eu + α = 0. So the three pairings correspond to the three roots of R(y) = y³ + ey + α. If P is irreducible, Frobenius acts on its roots as a 4-cycle, and a 4-cycle of S₄ fixes exactly one of the three pairings (its image in S₃ is a transposition). So R has exactly one root in F_q. For e = 0, R = y³ + α; when k is even, 3 | q − 1, so y³ = α has 0 or 3 solutions in F_q, never exactly one. If α = 0, P = (t + β^{1/4})⁴. ∎

For odd k ≥ 3 cubing is bijective on F_q and no d = 2 block is empty. The block counts are q = 8: 12 (×7), 28 (e = 0); q = 32: 240 (×31), 496 (e = 0). For q = 2 the empty block is e = 1. Script `char2_d2_blocks.py` (q ≤ 64).

**Sporadic case.** For q = 2, d = 4, the block {t⁸ + s : deg s ≤ 3} contains no irreducible (16 polynomials, checked). It gives the 8 empty sub-intervals f = t⁴ + a₂t² + a₁t + a₀ with c = a₂.

## 3. The census (`verify_pointwise_oppermann.py`, `data/oppermann_census.txt`)

Every Oppermann sub-interval for q^{2d} ≤ 2²⁶ was checked: a Ben-Or flag for every monic polynomial of degree 2d, summed in blocks. Empty sub-intervals occur only at

| q | d | empty / total | explanation |
|---|---|---|---|
| 2 | 1 | 2 / 4 | Prop. A |
| 2 | 2 | 4 / 8 | block t⁴ + t² + … |
| 2 | 4 | 8 / 32 | sporadic block t⁸ + (deg ≤ 3) |
| 4 | 1 | 4 / 16 | Prop. A |
| 4 | 2 | 16 / 64 | Prop. B |
| 8 | 1 | 8 / 64 | Prop. A |
| 16 | 1 | 16 / 256 | Prop. A |
| 16 | 2 | 256 / 4096 | Prop. B |

In every other computed case every sub-interval contains an irreducible:
- q = 2: d = 3, 5–13;
- q = 4: d = 3–6;
- q = 8: d = 2–4;
- q = 16: d = 3;
- **all odd q:** q = 3 (d ≤ 8), 5 (d ≤ 5), 7 (d ≤ 4), 9 (d ≤ 4), 11 (d ≤ 3), 13 (d ≤ 3).

The minimum count grows with d; for example q = 2, d = 13 has minimum 285 against mean 314. Summing over c recovers the Legendre counts of Part I, all positive.

## 4. Why the explicit formula cannot prove the statement (instruction 1)

The sub-intervals O_{f,c} are the classes of the Hayes group at level l = d (q^d characters). For a nontrivial character χ, L(u, χ) is a polynomial of degree ≤ d − 1 whose inverse roots have absolute value √q or 1 (Weil). The explicit formula therefore gives

  |Ψ(O_{f,c}) − q^d| ≤ (d − 1) q^d.

The main term q^d never beats the bound for d ≥ 2, whatever q is: Oppermann sits exactly on the square-root barrier, as it does over Z, where RH gives primes in [x, x + C√x log x] but not in [x, x + √x]. No degree threshold d₀ can come out of Weil II.

The rank-4 Witt sheaf and the sign law of Part IV do not enter here. They govern the *variance* of Ψ over Legendre intervals (the constant X), an average over f. A single interval is governed by the q^d Hayes characters above, not by the Gauss sums S″(λ). The variance controls how many intervals can deviate by a given amount. It gives no information about any individual interval, and the counterexamples of §2 show that no "measure-zero exceptional set" statement (instruction 2) can rule out empty intervals.

What is known in the other direction: for q → ∞ with the degree fixed, Galois-theoretic methods (Bank, Bary-Soroker and Rosenzweig, Duke Math. J. 2015) give asymptotics for primes in short intervals. Proposition B shows that in characteristic 2 any such result must exclude h = N/2 at N = 4 for q = 4^j. Whether their hypotheses already exclude it has to be checked against the paper.

## 5. The number-field transfer (instructions 4–5) is not a viable programme

Oppermann's conjecture over Z is open and beyond the reach of RH. The structures of Parts III–IV concern supersingular Artin–Schreier curves over F_{2^k}. Lifting such a Jacobian to characteristic 0 (Serre–Tate, Gross) produces an abelian scheme over a 2-adic ring, which carries no information about the distribution of rational primes in [n², (n+1)²]. Harman's sieve needs Type II bilinear bounds for sums over integers; no map takes Gauss sums over F_{q^m} to such bounds. I have not started this part. A formal "analogy framework" would be presented as a route to a theorem it cannot reach.

## 6. Realistic next steps

1. **Classify all empty Oppermann blocks in characteristic 2.** Are there finitely many for each q, and do they occur only for d ≤ 4? The resolvent argument suggests looking for blocks where the Frobenius cycle type forced by the shape of the polynomial (no odd coefficients above t^d) is incompatible with irreducibility.
2. **Odd characteristic.** Extend the census, for example q = 3 to d = 10 and q = 5 to d = 7, with the exact Hayes FFT of Part I at level l = d. State the verified range as a computer-assisted theorem, as Part I did for Legendre intervals.
3. **Large q.** Determine what the Galois-theoretic method gives for h = N/2 in odd characteristic, and whether Proposition B is the only obstruction in characteristic 2.
