"""
char2_d2_blocks.py -- which blocks B_e = {t^4 + e t^2 + a t + b : a, b in F_q} contain no irreducible, q = 2^k.

For d = 2 in characteristic 2 the Oppermann sub-interval O_{f,c} (f = t^2 + a1 t + a0) is the block B_e with
e = a1^2 + c (f^2 has no odd terms).  The census found q^2 empty sub-intervals for q = 2, 4, 16 and none for q = 8;
this script lists, for q = 2..256, the values e for which B_e is empty and the irreducible counts of all blocks.
Usage: python char2_d2_blocks.py [kmax=8]
"""
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from small_q_gf import make_field  # noqa: E402

IRRED = {1: [1, 1], 2: [1, 1, 1], 3: [1, 1, 0, 1], 4: [1, 1, 0, 0, 1], 5: [1, 0, 1, 0, 0, 1],
         6: [1, 1, 0, 0, 0, 0, 1], 7: [1, 1, 0, 0, 0, 0, 0, 1], 8: [1, 0, 1, 1, 1, 0, 0, 0, 1]}


def tables(k):
    q = 1 << k
    _, add, mul, neg, inv = make_field(2, k, IRRED[k])
    return q, np.array(mul, np.int64)


def trace(x, k, M):
    s, y = 0, x
    for _ in range(k):
        s ^= y
        y = int(M[y, y])
    return s


def roots_in_field(coeffs, q, M):
    """number of roots in F_q of the polynomial with coefficients (low first)"""
    cnt = 0
    for x in range(q):
        v, p = 0, 1
        for c in coeffs:
            if c:
                v ^= int(M[c, p])
            p = int(M[p, x])
        if v == 0:
            cnt += 1
    return cnt


def is_irreducible_quartic(e, a, b, q, M):
    """t^4 + e t^2 + a t + b over F_q (char 2): no root in F_q and no factorisation into two quadratics."""
    if b == 0:
        return False
    if roots_in_field([b, a, e, 0, 1], q, M):
        return False
    # (t^2 + u t + v)(t^2 + u t + w) since the t^3 coefficient is 0: u^2 + v + w = e, u(v + w) = a, v w = b
    for u in range(q):
        if u == 0:
            if a != 0:
                continue
            # v + w = e, v w = b: roots of z^2 + e z + b in F_q
            if roots_in_field([b, e, 1], q, M):
                return False
        else:
            # v + w = a / u, v w = b, and u^2 + a/u = e
            uinv = next(x for x in range(1, q) if M[u, x] == 1)
            s = int(M[a, uinv])
            if (int(M[u, u]) ^ s) != e:
                continue
            if roots_in_field([b, s, 1], q, M):
                return False
    return True


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for k in range(1, kmax + 1):
        q, M = tables(k)
        counts = {}
        for e in range(q):
            counts[e] = sum(1 for a in range(q) for b in range(q) if is_irreducible_quartic(e, a, b, q, M))
        empty = [e for e in counts if counts[e] == 0]
        dist = Counter(counts.values())
        print(f"q = {q:3d} (k = {k}): empty blocks e = {empty}  (Tr(e) = {[trace(e, k, M) for e in empty]});  "
              f"count distribution {dict(sorted(dist.items()))};  mean {sum(counts.values()) / q:.2f} = (q^2 - 1)/4 ?"
              f" {(q * q - 1) / 4:.2f}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
