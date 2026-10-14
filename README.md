# Oppermann intervals in F_q[t]

**Galois obstructions, exact counts for small degree, and the Oppermann subgroup in characteristic 2**

Ruqing Chen. DOI: [10.5281/zenodo.23353949](https://doi.org/10.5281/zenodo.23353949)

Paper: `paper/oppermann_function_fields.tex` / `.pdf` (22 pages).

For monic f of degree d, the Legendre interval {fÂ² + s : deg s â‰¤ d} splits into q **Oppermann sub-intervals**

  O_{f,c} = {fÂ² + cÂ·f + s : deg s < d},  c âˆˆ F_q,

each of size q^d. For q = 2 these are the two halves, the analogue of (nÂ², nÂ²+n) and (nÂ²+n, (n+1)Â²). Does every O_{f,c} contain an irreducible polynomial?

## Main results

**Characteristic 2: no. The failures are explained by three Galois obstructions (Â§3).**
- **(A) d = 1:** the c = 0 piece consists of squares.
- **(B) d = 2:** Frobenius acts on the roots of an additive polynomial through AGLâ‚‚(Fâ‚‚) â‰… Sâ‚„. The exact counts are qÂ²/4, 0, q(qâˆ’2)/4 and q(qâˆ’1)/2. Empty sub-intervals occur exactly for q = 2 and q = 4^j.
- **(C) q = 2, d = 4:** the Stickelbergerâ€“Swan sign of Frobenius obstructs.
- **Dead zone (Â§3.4):** for r â‰¥ 3, AGL_r(Fâ‚‚) contains no 2^r-cycle. Hence no affine polynomial of degree 2^r â‰¥ 8 is irreducible, and the mechanism of (B) does not extend to d â‰¥ 4.

**Exact counts for small degree**
- **Odd q, d = 2 (Â§4):** a ternary trace form on F_{qâ´} gives Ï€ = (qÂ² âˆ’ 1 âˆ’ (q âˆ’ Î·(âˆ’1))Î·(2e))/4 > 0.
- **Characteristic 2, d = 3 (Â§5):** a count in the field tower F_q âŠ‚ F_{qÂ³} âŠ‚ F_{qâ¶}, together with a cubic Gauss sum, gives every class. Let c be the tÂ³-coefficient and Îµ = (âˆ’1)^{k/2}.
  - k odd: Ï€ = q(qÂ² âˆ’ 1)/6 for every class.
  - k even, c = 0: Ï€ = q(q âˆ’ 1)Â²/6.
  - k even, c a cube: Ï€ = q(qÂ² + 1 âˆ’ 2Îµâˆšq)/6.
  - k even, c a non-cube: Ï€ = q(qÂ² + 1 + Îµâˆšq)/6.

**Large q and computation (Â§6)**
- **Large q:** for fixed d â‰¥ 3, Bankâ€“Bary-Sorokerâ€“Rosenzweig together with the d = 3 count gives nonemptiness once q is large.
- **Odd q:** exact Hayes-FFT counts find every sub-interval nonempty for q = 3, d â‰¤ 16; q = 5, d â‰¤ 10; q = 7, d â‰¤ 9; q = 9, d â‰¤ 8; q = 11, d â‰¤ 7; q = 13, d â‰¤ 6.
- **Characteristic 2:** the same computation finds no empty sub-intervals beyond (A)â€“(C) for q = 2, d â‰¤ 26; q = 4, d â‰¤ 13; q = 8, d â‰¤ 7; q = 16, d â‰¤ 5.
- **Weil bound:** it sits at the square-root barrier. At d = 2 every level-2 L-function in characteristic 2 is rigid, with Î³â´ = (âˆ’1)^k qÂ².

**The Oppermann subgroup in characteristic 2 (Â§2.1, Â§7)**
- **Structure:** H_d = ker(P â†¦ Sâ‚, Sâ‚ƒ, â€¦), so G_d/H_d â‰… F_q^{âŒŠd/2âŒ‹}. The annihilator's explicit-formula coefficients are Artinâ€“Schreier sums.
- **Mean drift:** E_H[Î¨] âˆ’ q^d is 0 for d = 2, 3 and for d = 4 with k even. It is âˆ’2qÂ² + 2q for d = 4 with k odd, and qÂ² âˆ’ q for d = 5.
- **Variance:** it is a coherent sum over annihilator cosets.
  - d = 2: Var_H/qÂ² = q âˆ’ 1, against Var_G/qÂ² = 1 âˆ’ 1/q.
  - d = 3: Var_H/qÂ³ = 1 âˆ’ 1/q (k odd) or 3(1 âˆ’ 1/q) (k even).
  - For d â‰¥ 4 there is a numerical census, compared with the Keatingâ€“Rudnick value d âˆ’ 1.

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
| `hayes_fft_fq.py`, `fqlib.py`, `small_q_gf.py` | â€” | from Part I. The only change is extra moduli in `fqlib.py` for q = 32, â€¦, 256 |

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
