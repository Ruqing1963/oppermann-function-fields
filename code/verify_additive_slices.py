"""
verify_additive_slices.py -- Theorem 6.1 and Corollary 6.2 (additive slices, AGL_r(F_2)).

(a) For r = 1..4: the orders of all elements of AGL_r(F_2) acting on F_2^r; whether a 2^r-cycle exists; the largest
    2-power order, compared with 2^{ceil(log2(r+1))}.
(b) For degree 2^r = 8 (q = 2, 4, 8, 16) and 2^r = 16 (q = 2, 4): every affine polynomial
    t^{2^r} + c_{r-1} t^{2^{r-1}} + ... + c_0 t + beta over F_q is reducible.  For comparison, degree 4 (r = 2) and
    degree 2 (r = 1) affine polynomials do contain irreducibles.
Usage: python verify_additive_slices.py > ../data/additive_slices.txt
"""
import itertools
import math
import os
import sys

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fqlib  # noqa: E402


def agl_orders(r):
    """multiset of (order, is a 2^r-cycle) over AGL_r(F_2); matrices as tuples of column bitmasks."""
    n = 1 << r
    cols = list(range(1, n))
    count_cycles = 0
    max2 = 1
    total = 0
    for colset in itertools.permutations(range(1, n), r):
        # invertibility: columns span F_2^r
        span = {0}
        for c in colset:
            span |= {s ^ c for s in span}
        if len(span) != n:
            continue

        def F(w):
            out = 0
            for i in range(r):
                if w >> i & 1:
                    out ^= colset[i]
            return out
        Fw = [F(w) for w in range(n)]
        for v in range(n):
            total += 1
            perm = [Fw[w] ^ v for w in range(n)]
            # order and cycle type
            seen = [False] * n
            order = 1
            longest = 0
            for s in range(n):
                if seen[s]:
                    continue
                L, x = 0, s
                while not seen[x]:
                    seen[x] = True
                    x = perm[x]
                    L += 1
                order = order * L // math.gcd(order, L)
                longest = max(longest, L)
            if longest == n:
                count_cycles += 1
            if order & (order - 1) == 0:
                max2 = max(max2, order)
    del cols
    return total, count_cycles, max2


@njit(cache=True)
def _count_irreducible_affine(r, q, A, M, NEG, INV):
    n = 1 << r
    W = np.zeros((8, 2 * n + 2), np.int64)
    P = np.zeros(n + 1, np.int64)
    cnt = 0
    tot = 0
    ncoef = r + 1  # c_0..c_{r-1}, beta
    for code in range(q ** ncoef):
        for i in range(n + 1):
            P[i] = 0
        P[n] = 1
        x = code
        for i in range(r):
            P[1 << i] = x % q
            x //= q
        P[0] = x % q
        tot += 1
        if P[0] == 0:
            continue
        if fqlib._is_irreducible(P, n, q, A, M, NEG, INV, W):
            cnt += 1
    return tot, cnt


def main():
    print("(a) AGL_r(F_2)")
    for r in range(1, 5):
        total, cyc, max2 = agl_orders(r)
        bound = 1 << math.ceil(math.log2(r + 1))
        print(f"  r = {r}: |AGL| = {total}, 2^r-cycles: {cyc}, largest 2-power order {max2}, "
              f"bound 2^ceil(log2(r+1)) = {bound}")
        sys.stdout.flush()
    print("(b) affine polynomials t^{2^r} + sum c_i t^{2^i} + beta over F_q: irreducible / total")
    for r, qs in ((1, (2, 4, 8, 16)), (2, (2, 4, 8, 16)), (3, (2, 4, 8, 16)), (4, (2, 4))):
        for q in qs:
            A, M, NEG, INV = fqlib.field_tables(q)
            tot, cnt = _count_irreducible_affine(r, q, A, M, NEG, INV)
            print(f"  degree {1 << r:2d}, q = {q:2d}: {cnt} / {tot}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
