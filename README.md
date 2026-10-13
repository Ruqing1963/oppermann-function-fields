# Oppermann intervals in F_q[t]

Independent of the Legendre-interval repositories, Parts I–IV.

For monic f of degree d, the Legendre interval {f² + s : deg s ≤ d} splits into q **Oppermann sub-intervals**

  O_{f,c} = {f² + c·f + s : deg s < d},  c ∈ F_q,

each of size q^d. For q = 2 these are the two halves, the analogue of (n², n²+n) and (n²+n, (n+1)²). Does every O_{f,c} contain an irreducible polynomial?

The preprint draft is `paper/oppermann_function_fields.tex`.

## Results

**Characteristic 2: no, and the failures are explained completely (proved).**
- **(A) d = 1:** the c = 0 piece consists of squares.
- **(B) d = 2:** O_{f,c} = {t⁴ + e t² + αt + β}, with e = a₁² + c. Its exact number of irreducibles is:

  | | e = 0 | e ≠ 0 |
  |---|---|---|
  | k even | 0 | q²/4 |
  | k odd | q(q−1)/2 | q(q−2)/4 |

  The proof lets Frobenius act on the roots of the additive polynomial t⁴ + e t² + αt through AGL₂(F₂) ≅ S₄. The criterion is Tr(e³/α²) ≠ Tr(1). Consequently there are empty sub-intervals exactly for q = 2 and q = 4^j.
- **(C) q = 2, d = 4:** the block {t⁸ + s : deg s ≤ 3} contains no irreducible. This follows from the Stickelberger–Swan sign of Frobenius: when 8 | deg, the sign is decided by the tail alone.

**Odd characteristic, d = 2: proved for every q.** π(O_{f,c}) = (q² − 1 − (q − η(−1))·η(2e))/4 > 0. The proof uses the ternary trace form Tr(x²) on ker Tr ⊂ F_{q⁴}.

**Characteristic 2, d = 3, p₃ = 0 (the classes BBR does not cover): proved for every q.** All q classes are equivalent under translation. Each contains q(q²−1)/6 irreducibles (k odd) or q(q−1)²/6 (k even), by a trace computation in F_q ⊂ F_{q³} ⊂ F_{q⁶}.

**Large q:** for every d ≥ 3, all Oppermann sub-intervals are nonempty once q ≥ q₀(d). This follows from Bank–Bary-Soroker–Rosenzweig together with the d = 3, p₃ = 0 theorem; q₀(d) is not explicit.

**Odd characteristic: yes, in the whole computed range (computer-assisted).**
- Every Hayes class of level d is an Oppermann sub-interval.
- Exact counts come from the Hayes explicit formula, using the FFT algorithm of Part I at level l = d.
- Every sub-interval contains irreducibles for q = 3, d ≤ 16; q = 5, d ≤ 10; q = 7, d ≤ 9; q = 9, d ≤ 8; q = 11, d ≤ 7; q = 13, d ≤ 6.

**Characteristic 2, computed range:** there are no empty Oppermann sub-intervals beyond (A)–(C) for q = 2, d ≤ 26; q = 4, d ≤ 13; q = 8, d ≤ 7; q = 16, d ≤ 5.

**Limits of the method:** counted by exact level, the Weil bound gives |Ψ − q^d| ≤ (d − 1)q^d − Σ_{j<d} q^j.
- For d ≥ 3 this is the square-root barrier, so it proves nothing pointwise.
- For d = 2 it gives Ψ ≥ q. In the empty classes of (B) this is attained exactly, by fourth powers.

## Files

| file | purpose |
|---|---|
| `code/hayes_oppermann.py` | exact counts in all Oppermann sub-intervals (Hayes explicit formula + FFT) → `data/hayes_oppermann_census.txt` |
| `code/verify_pointwise_oppermann.py` | brute-force Ben-Or census, q^{2d} ≤ 2²⁶ → `data/oppermann_census.txt` (cross-check) |
| `code/char2_d2_blocks.py`, `code/char2_d2_theorem.py` | d = 2 in characteristic 2: block counts and the exact formula, q ≤ 64 |
| `code/swan_parity.py` | zero blocks {t^{2d} + s} over F₂: factorisation, disc mod 8, Stickelberger check → `data/swan_parity.txt` |
| `code/make_figure.py` | `figures/oppermann_min_ratio.pdf` |
| `code/hayes_fft_fq.py`, `code/fqlib.py`, `code/small_q_gf.py` | finite-field and Hayes-group routines copied from Part I (unchanged) |
| `notes/REPORT.md` | working notes |

Run with:

```
cd code
python hayes_oppermann.py 3 1 16 <tmpdir>      # q d_min d_max [memmap dir]
python char2_d2_theorem.py 6
python swan_parity.py 10
python make_figure.py
cd ../paper && pdflatex oppermann_function_fields.tex
```
