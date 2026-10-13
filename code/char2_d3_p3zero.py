"""
char2_d3_p3zero.py -- the d = 3 Oppermann classes with p_3 = 0 in characteristic 2:
    C_e = {t^6 + e t^4 + s2 t^2 + s1 t + s0 : s2, s1, s0 in F_q},   q = 2^k.
These are the classes not covered by Bank--Bary-Soroker--Rosenzweig (every member has P' = s1, of degree <= 0).
Translation t -> t + a maps C_e onto C_{e + a^2} ((t+a)^6 = t^6 + a^2 t^4 + a^4 t^2 + a^6), so pi(C_e) does not depend on e;
the script checks this for q <= 16 and compares pi(C_0) with the conjectured exact value
    q(q^2 - 1)/6 (k odd),   q(q - 1)^2/6 (k even).
Usage: python char2_d3_p3zero.py [kmax=6]
"""
import os
import sys

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from small_q_gf import make_field  # noqa: E402

IRRED = {1: [1, 1], 2: [1, 1, 1], 3: [1, 1, 0, 1], 4: [1, 1, 0, 0, 1], 5: [1, 0, 1, 0, 0, 1],
         6: [1, 1, 0, 0, 0, 0, 1], 7: [1, 1, 0, 0, 0, 0, 0, 1], 8: [1, 0, 1, 1, 1, 0, 0, 0, 1]}
N = 6


def tables(k):
    _, add, mul, neg, inv = make_field(2, k, IRRED[k])
    return np.array(mul, np.int64), np.array(inv, np.int64)


@njit(cache=True)
def _mulmod(x, y, P, M):
    r = np.zeros(2 * N - 1, np.int64)
    for i in range(N):
        if x[i]:
            for j in range(N):
                if y[j]:
                    r[i + j] ^= M[x[i], y[j]]
    for m in range(2 * N - 2, N - 1, -1):
        c = r[m]
        if c:
            for i in range(N + 1):
                r[m - N + i] ^= M[c, P[i]]
    return r[:N].copy()


@njit(cache=True)
def _powq(x, q, P, M):
    res = np.zeros(N, np.int64)
    res[0] = 1
    base = x.copy()
    e = q
    while e:
        if e & 1:
            res = _mulmod(res, base, P, M)
        base = _mulmod(base, base, P, M)
        e >>= 1
    return res


@njit(cache=True)
def _gcd_deg(P, g, M, INV):
    a = np.zeros(N + 1, np.int64)
    a[:] = P
    b = np.zeros(N + 1, np.int64)
    b[:N] = g
    da = N
    db = N - 1
    while db >= 0 and b[db] == 0:
        db -= 1
    while db >= 0:
        il = INV[b[db]]
        while da >= db:
            c = M[a[da], il]
            if c:
                for i in range(db + 1):
                    a[da - db + i] ^= M[c, b[i]]
            da -= 1
            while da >= 0 and a[da] == 0:
                da -= 1
        tmp = a
        a = b
        b = tmp
        t2 = da
        da = db
        db = t2
    return da


@njit(cache=True)
def _type(P, q, M, INV):
    """number of distinct-degree factors: returns (deg of gcd with t^{q^i}-t for i=1,2,3); irreducible iff all 0
    and P squarefree-free of these."""
    x = np.zeros(N, np.int64)
    x[1] = 1
    out = np.zeros(3, np.int64)
    for i in range(3):
        x = _powq(x, q, P, M)
        g = x.copy()
        g[1] ^= 1
        out[i] = _gcd_deg(P, g, M, INV)
    return out


@njit(cache=True)
def count_block(e, q, M, INV):
    P = np.zeros(N + 1, np.int64)
    P[6] = 1
    P[4] = e
    irr = 0
    for s2 in range(q):
        for s1 in range(q):
            for s0 in range(q):
                P[2] = s2
                P[1] = s1
                P[0] = s0
                o = _type(P, q, M, INV)
                if o[0] == 0 and o[1] == 0 and o[2] == 0:
                    irr += 1
    return irr


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    print(f"{'q':>4} {'pi(C_0)':>9} {'pi(C_1)':>9} {'q^3/6':>10} {'(C_0-q^3/6)/q^2':>16} {'(C_1-q^3/6)/q^2':>16}"
          f" {'all e equal':>12} {'formula':>8}")
    for k in range(1, kmax + 1):
        q = 1 << k
        M, INV = tables(k)
        c0 = count_block(0, q, M, INV)
        c1 = count_block(1, q, M, INV) if q <= 64 else c0
        eq = True
        if q <= 16:
            eq = all(count_block(e, q, M, INV) == c1 for e in range(2, q))
        print(f"{q:>4} {c0:>9} {c1:>9} {q ** 3 / 6:>10.1f} {(c0 - q ** 3 / 6) / q ** 2:>16.4f} "
              f"{(c1 - q ** 3 / 6) / q ** 2:>16.4f} {str(eq) if q <= 16 else '-':>12} "
              f"{str(c0 == (q * (q * q - 1) // 6 if k % 2 else q * (q - 1) ** 2 // 6)):>8}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
