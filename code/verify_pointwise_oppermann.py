"""
verify_pointwise_oppermann.py -- exhaustive census of irreducibles in the Oppermann sub-intervals of F_q[t].

For monic f of degree d the Legendre interval I_f = {f^2 + s : deg s <= d} (Parts I-IV) splits into the q
sub-intervals
        O_{f,c} = {f^2 + c f + s : deg s < d} = {f^2 + c t^d + s : deg s < d},    c in F_q,
each with q^d elements (the analogue of (n^2, n^2+n) and (n^2+n, (n+1)^2); for q = 2 there are exactly two
halves, c = 0 and c = 1).  The function-field Oppermann statement is: every O_{f,c} contains an irreducible.

Method: a Ben-Or irreducibility flag for every monic polynomial of degree 2d (q^{2d} of them, numba, parallel),
reshaped so that each row is one block {prefix + s : deg s < d}; the row sums are the irreducible counts of all
short intervals of degree 2d and length q^d at once.  The Oppermann blocks are the rows whose prefix (coefficients
of t^d .. t^{2d-1}) is that of f^2 + c t^d.  Cross-check: summing over c gives the Legendre counts N_irr(f) of
Part I.  Output: one line per (q, d) with the minimum count, the number of empty Oppermann sub-intervals, and the
list of empty ones (f, c) when there are few.

Usage: python verify_pointwise_oppermann.py [max_log2_size=24]
"""
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fqlib  # noqa: E402

FIELDS = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16]


@njit(parallel=True, cache=True)
def irreducible_flags(n, q, A, M, NEG, INV):
    total = q ** n
    out = np.zeros(total, np.uint8)
    nchunks = 256
    per = (total + nchunks - 1) // nchunks
    for c in prange(nchunks):
        W = np.zeros((8, 2 * n + 2), np.int64)
        P = np.zeros(n + 1, np.int64)
        lo = c * per
        hi = min(total, lo + per)
        for idx in range(lo, hi):
            r = idx
            for i in range(n):
                P[i] = r % q
                r //= q
            P[n] = 1
            if P[0] == 0:        # divisible by t
                continue
            if fqlib._is_irreducible(P, n, q, A, M, NEG, INV, W):
                out[idx] = 1
    return out


def prefix_index(coeffs_top, q):
    """coeffs_top[i] = coefficient of t^{d+i}, i = 0..d-1 -> integer index (digit i = coefficient of t^{d+i})"""
    v = 0
    for i in reversed(range(len(coeffs_top))):
        v = v * q + int(coeffs_top[i])
    return v


def run(q, d):
    A, M, NEG, INV = fqlib.field_tables(q)
    n = 2 * d
    t0 = time.time()
    flags = irreducible_flags(n, q, A, M, NEG, INV)
    # index = low + q^d * prefix, where low encodes t^0..t^{d-1} and prefix encodes t^d..t^{2d-1}
    counts = flags.reshape(q ** d, q ** d).sum(axis=1, dtype=np.int64)   # counts[prefix]
    total_irr = int(flags.sum())
    del flags
    rows = []
    empty = []
    legendre = []
    for fidx in range(q ** d):
        f = [0] * (d + 1)
        r = fidx
        for i in range(d):
            f[i] = r % q
            r //= q
        f[d] = 1
        sq = fqlib.py_mul(f, f, A, M)                     # length 2d+1, low first
        leg = 0
        for c in range(q):
            top = list(sq[d:2 * d])
            top[0] = int(A[top[0], c])                     # + c t^d
            cnt = int(counts[prefix_index(top, q)])
            rows.append(cnt)
            leg += cnt
            if cnt == 0:
                empty.append((tuple(f), c))
        legendre.append(leg)
    rows = np.array(rows)
    return dict(q=q, d=d, n_sub=len(rows), min=int(rows.min()), mean=float(rows.mean()), n_empty=len(empty),
                empty=empty, leg_min=int(min(legendre)), total_irr=total_irr, secs=time.time() - t0,
                expected=q ** d / n)


def fmt_poly(f, q):
    terms = []
    for i in reversed(range(len(f))):
        c = f[i]
        if c == 0:
            continue
        mon = "1" if i == 0 else ("t" if i == 1 else f"t^{i}")
        terms.append(mon if (c == 1 and i > 0) else (f"{c}" if i == 0 else f"{c}{mon}"))
    return "+".join(terms)


def main():
    max_log2 = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    print(f"{'q':>3} {'d':>3} {'#O_f,c':>8} {'min':>5} {'mean':>9} {'q^d/2d':>9} {'#empty':>7} {'Legendre min':>12} {'secs':>7}")
    for q in FIELDS:
        d = 1
        while (2 * d) * np.log2(q) <= max_log2:
            res = run(q, d)
            print(f"{q:>3} {d:>3} {res['n_sub']:>8} {res['min']:>5} {res['mean']:>9.2f} {res['expected']:>9.2f} "
                  f"{res['n_empty']:>7} {res['leg_min']:>12} {res['secs']:>7.1f}")
            if 0 < res["n_empty"] <= 12:
                for f, c in res["empty"]:
                    print(f"        empty: f = {fmt_poly(f, q)}, c = {c}")
            sys.stdout.flush()
            d += 1


if __name__ == "__main__":
    main()
