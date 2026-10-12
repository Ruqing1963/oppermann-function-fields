# Oppermann intervals in F_q[t]

Independent of the Legendre-interval repositories, Parts I–IV.

The question: does every Oppermann sub-interval O_{f,c} = {f² + c·f + s : deg s < d} contain an irreducible polynomial? Here f is monic of degree d over F_q, c ∈ F_q, and for q = 2 the two values of c give the two Oppermann halves.

## Findings so far

The answer is **no in characteristic 2**:
- d = 1, every even q: the c = 0 piece consists of squares.
- d = 2, q = 4^j: there are q² empty sub-intervals, by a resolvent-cubic argument.
- q = 2, d = 4: one sporadic empty block.

The answer is **yes in every other computed case**: all q^{2d} ≤ 2²⁶, including q = 3, 5, 7, 9, 11, 13 for every degree in range.

The Weil/explicit-formula bound is exactly at the square-root barrier, |Ψ − q^d| ≤ (d − 1)q^d, so it cannot prove the statement for any d ≥ 2. See `notes/REPORT.md`.

## Files

- `code/verify_pointwise_oppermann.py` → `data/oppermann_census.txt`: an exhaustive census of irreducibles in every sub-interval, with numba and a Ben-Or irreducibility test.
- `code/char2_d2_blocks.py`: the d = 2 blocks {t⁴ + e t² + αt + β} in characteristic 2.
- `code/fqlib.py`, `code/small_q_gf.py`: finite-field tables, copied from Part I.

Run with:

```
cd code
python verify_pointwise_oppermann.py 26 > ../data/oppermann_census.txt
python char2_d2_blocks.py 6
```
