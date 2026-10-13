"""
odd_d2_formula.py -- brute-force check of the exact d = 2 count in odd characteristic.

Theorem (paper, Thm 3.x).  q odd, B_e = {t^4 + e t^2 + a t + b : a, b in F_q}.  Every level-2 Hayes class in degree 4
(= every d = 2 Oppermann sub-interval) is carried bijectively onto some B_e by t -> t - p_3/4, and
    pi(B_0) = (q^2 - 1)/4,
    pi(B_e) = (q^2 - q*eta(2e) - 1 + eta(-2e))/4      (e != 0),
eta the quadratic character of F_q.  In particular pi >= (q - 2)(q + 1)/4 > 0.

Brute force: P is irreducible iff gcd(t^{q^2} - t, P) = 1 (no factor of degree 1 or 2).  Checked for every e and
every (a, b) for q in {3, 5, 7, 9, 11, 13, 25, 27, 49, 81}.  Also checks the reduction: for q <= 13, every class
(p3, p2) has the count of B_e with e = p2 - 3 p3^2 / 8.
Usage: python odd_d2_formula.py
"""
import os
import sys

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from small_q_gf import make_field  # noqa: E402

# (p, k, monic irreducible modulus of degree k over F_p, low first)
FIELDS = {3: (3, 1, [0, 1]), 5: (5, 1, [0, 1]), 7: (7, 1, [0, 1]), 9: (3, 2, [1, 0, 1]), 11: (11, 1, [0, 1]),
          13: (13, 1, [0, 1]), 25: (5, 2, [2, 0, 1]), 27: (3, 3, [1, 2, 0, 1]), 49: (7, 2, [1, 0, 1]),
          81: (3, 4, [2, 1, 0, 0, 1])}


def tables(q):
    p, k, mod = FIELDS[q]
    _, add, mul, neg, inv = make_field(p, k, mod)
    A, M, N, I = (np.array(x, np.int64) for x in (add, mul, neg, inv))
    assert all(M[a, I[a]] == 1 for a in range(1, q)), "modulus not irreducible"
    return A, M, N, I


@njit(cache=True)
def _mulmod(x, y, P, A, M, N):
    """x, y: length-4 arrays (deg < 4); P monic quartic (length 5). returns x*y mod P."""
    r = np.zeros(7, np.int64)
    for i in range(4):
        if x[i]:
            for j in range(4):
                if y[j]:
                    r[i + j] = A[r[i + j], M[x[i], y[j]]]
    for m in range(6, 3, -1):
        c = r[m]
        if c:
            for i in range(5):
                r[m - 4 + i] = A[r[m - 4 + i], N[M[c, P[i]]]]
    return r[:4].copy()


@njit(cache=True)
def _irreducible(P, q, A, M, N, I):
    # t^(q^2) mod P by repeated q-th powering
    x = np.zeros(4, np.int64)
    x[1] = 1
    for _ in range(2):
        # x <- x^q
        res = np.zeros(4, np.int64)
        res[0] = 1
        base = x.copy()
        e = q
        while e:
            if e & 1:
                res = _mulmod(res, base, P, A, M, N)
            base = _mulmod(base, base, P, A, M, N)
            e >>= 1
        x = res
    # g = x - t
    g = x.copy()
    g[1] = A[g[1], N[1]]
    # gcd(P, g)
    a = np.zeros(5, np.int64)
    a[:] = P
    b = np.zeros(5, np.int64)
    b[:4] = g
    da = 4
    db = 3
    while db >= 0 and b[db] == 0:
        db -= 1
    while db >= 0:
        # a <- a mod b
        inv_lead = I[b[db]]
        while da >= db and da >= 0:
            c = M[a[da], inv_lead]
            if c:
                for i in range(db + 1):
                    a[da - db + i] = A[a[da - db + i], N[M[c, b[i]]]]
            da -= 1
            while da >= 0 and a[da] == 0:
                da -= 1
        a, b = b, a
        da, db = db, da
    return da == 0          # gcd is a nonzero constant


@njit(cache=True)
def block_count(e, q, A, M, N, I):
    P = np.zeros(5, np.int64)
    P[4] = 1
    P[2] = e
    cnt = 0
    for a in range(q):
        for b in range(q):
            P[1] = a
            P[0] = b
            if _irreducible(P, q, A, M, N, I):
                cnt += 1
    return cnt


def main():
    for q in FIELDS:
        A, M, N, I = tables(q)
        squares = {int(M[x, x]) for x in range(1, q)}
        eta = lambda x: 0 if x == 0 else (1 if x in squares else -1)  # noqa: E731
        two = int(A[1, 1])
        bad = 0
        counts = {}
        for e in range(q):
            c = int(block_count(e, q, A, M, N, I))
            counts[e] = c
            if e == 0:
                pred = (q * q - 1) // 4
            else:
                x = int(M[two, e])
                pred = (q * q - q * eta(x) - 1 + eta(int(N[x]))) // 4
            bad += c != pred
        vals = sorted(set(counts.values()))
        line = f"q = {q:3d}: blocks B_e checked {q}, formula violations {bad}, values {vals}, min {min(vals)}"
        if q <= 13:
            # reduction check: class (p3, p2) -> e = p2 - 3 p3^2 / 8, count by translating every member
            inv8 = int(I[int(A[int(A[int(A[1, 1]), int(A[1, 1])]), int(A[int(A[1, 1]), int(A[1, 1])])])])
            three = int(A[two, 1])
            red_bad = 0
            for p3 in range(q):
                for p2 in range(q):
                    e = int(A[p2, N[M[M[three, M[p3, p3]], inv8]]])
                    cnt = 0
                    P = np.zeros(5, np.int64)
                    P[4], P[3], P[2] = 1, p3, p2
                    for a in range(q):
                        for b in range(q):
                            P[1], P[0] = a, b
                            cnt += _irreducible(P, q, A, M, N, I)
                    red_bad += cnt != counts[e]
            line += f", reduction violations {red_bad}"
        print(line)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
