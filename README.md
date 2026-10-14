# Oppermann intervals in F_q[t]

**Galois obstructions, exact counts for small degree, and the Oppermann subgroup in characteristic 2**

Ruqing Chen. DOI: [10.5281/zenodo.23353949](https://doi.org/10.5281/zenodo.23353949)

Paper: `paper/oppermann_function_fields.tex` / `.pdf` (23 pages).

For monic f of degree d, the Legendre interval {f² + s : deg s ≤ d} splits into q **Oppermann sub-intervals**

  O_{f,c} = {f² + c·f + s : deg s < d},  c ∈ F_q,

each of size q^d. For q = 2 these are the two halves, the analogue of (n², n²+n) and (n²+n, (n+1)²). Does every O_{f,c} contain an irreducible polynomial?

## Main results

**Characteristic 2: no. The failures are explained by three Galois obstructions (§3).**
- **(A) d = 1:** the c = 0 piece consists of squares.
- **(B) d = 2:** Frobenius acts on the roots of an additive polynomial through AGL₂(F₂) ≅ S₄. The exact counts are q²/4, 0, q(q−2)/4 and q(q−1)/2. Empty sub-intervals occur exactly for q = 2 and q = 4^j.
- **(C) q = 2, d = 4:** the Stickelberger–Swan sign of Frobenius obstructs.
- **Dead zone (§3.4):** for r ≥ 3, AGL_r(F₂) contains no 2^r-cycle. Hence no affine polynomial of degree 2^r ≥ 8 is irreducible, and the mechanism of (B) does not extend to d ≥ 4.

**Exact counts for small degree**
- **Odd q, d = 2 (§4):** a ternary trace form on F_{q⁴} gives π = (q² − 1 − (q − η(−1))η(2e))/4 > 0.
- **Characteristic 2, d = 3 (§5):** a count in the field tower F_q ⊂ F_{q³} ⊂ F_{q⁶}, together with a cubic Gauss sum, gives every class. Let c be the t³-coefficient and ε = (−1)^{k/2}.
  - k odd: π = q(q² − 1)/6 for every class.
  - k even, c = 0: π = q(q − 1)²/6.
  - k even, c a cube: π = q(q² + 1 − 2ε√q)/6.
  - k even, c a non-cube: π = q(q² + 1 + ε√q)/6.

**Large q and computation (§6)**
- **Large q:** for fixed d ≥ 3, Bank–Bary-Soroker–Rosenzweig together with the d = 3 count gives nonemptiness once q is large.
- **Odd q:** exact Hayes-FFT counts find every sub-interval nonempty for q = 3, d ≤ 16; q = 5, d ≤ 10; q = 7, d ≤ 9; q = 9, d ≤ 8; q = 11, d ≤ 7; q = 13, d ≤ 6.
- **Characteristic 2:** the same computation finds no empty sub-intervals beyond (A)–(C) for q = 2, d ≤ 26; q = 4, d ≤ 13; q = 8, d ≤ 7; q = 16, d ≤ 5.
- **Weil bound:** it sits at the square-root barrier. At d = 2 every level-2 L-function in characteristic 2 is rigid, with γ⁴ = (−1)^k q².

**The Oppermann subgroup in characteristic 2 (§2.1, §7)**
- **Structure:** H_d = ker(P ↦ S₁, S₃, …), so G_d/H_d ≅ F_q^{⌊d/2⌋}. The annihilator's explicit-formula coefficients are Artin–Schreier sums.
- **Mean drift:** E_H[Ψ] − q^d is 0 for d = 2, 3 and for d = 4 with k even. It is −2q² + 2q for d = 4 with k odd, and q² − q for d = 5.
- **Variance:** it is a coherent sum over annihilator cosets.
  - d = 2: Var_H/q² = q − 1, against Var_G/q² = 1 − 1/q.
  - d = 3: Var_H/q³ = 1 − 1/q (k odd) or 3(1 − 1/q) (k even).
  - For d ≥ 4 there is a numerical census, compared with the Keating–Rudnick value d − 1.

## Files

| script (`code/`) | output (`data/`, `figures/`, `paper/`) | paper |
|---|---|---|
| `hayes_oppermann.py` | `hayes_oppermann_census.txt` | Thms 1.2, 1.3 |
| `verify_pointwise_oppermann.py`, `compare_bruteforce.py` | `oppermann_census.txt`, `compare_bruteforce.txt` | cross-checks |
| `char2_d2_blocks.py`, `char2_d2_theorem.py` | `char2_d2_blocks.txt`, `char2_d2_theorem.txt` | Thm 1.1(B) |
| `swan_parity.py` | `swan_parity.txt` | Thm 1.1(C), Table 1 |
| `verify_additive_slices.py` | `additive_slices.txt` | Thm 3.6, Cor 3.7 |
| `odd_d2_formula.py` | `odd_d2_formula.txt` | Thm 1.4 |
| `char2_d3_p3zero.py`, `char2_d3_all.py` | `char2_d3_p3zero.txt`, `char2_d3_all.txt` | Thm 5.2 |
| `cubic_gauss_m.py`, `make_d3_table.py` | `cubic_gauss_m.txt`, `paper/tab_d3.tex` | Thm 5.3, Cor 5.4, Table 2 |
| `char2_d2_rigidity.py` | `char2_d2_rigidity.txt` | Prop 6.1 |
| `subgroup_lib.py`, `verify_subgroup_mean.py` | `verify_subgroup_mean.txt` | Thms 2.4, 2.5, 7.3 |
| `subgroup_variance_census.py`, `verify_variance23.py` | `subgroup_variance.{txt,json}`, `verify_variance23.txt` | Thm 7.5 |
| `make_figure.py`, `make_figure2.py`, `make_tables.py` | `figures/oppermann_min_ratio.pdf`, `figures/oppermann_char2_rigidity.pdf`, `paper/census_tables.tex` | Figs 1, 2 |
| `make_subgroup_figures.py`, `make_subgroup_tables.py` | `figures/subgroup_mean_gauss.pdf`, `figures/subgroup_variance.pdf`, `paper/tab_subgroup_*.tex` | Figs 3, 4 |
| `hayes_fft_fq.py`, `fqlib.py`, `small_q_gf.py` | — | from Part I. The only change is extra moduli in `fqlib.py` for q = 32, …, 256 |

`notes/REPORT.md` contains the working notes of the first round.

Run with:

```
cd code
python hayes_oppermann.py 3 1 16 <tmpdir>      # q d_min d_max [memmap dir]
python verify_subgroup_mean.py > ../data/verify_subgroup_mean.txt
python subgroup_variance_census.py <tmpdir>
python cubic_gauss_m.py 14 > ../data/cubic_gauss_m.txt
python make_figure.py; python make_figure2.py; python make_subgroup_figures.py
python make_tables.py; python make_d3_table.py; python make_subgroup_tables.py
cd ../paper && pdflatex oppermann_function_fields.tex
```
