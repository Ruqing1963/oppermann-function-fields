"""
char2_d3_all.py -- all d = 3 Oppermann classes in characteristic 2 (paper, Theorem 6.4).

C_{e,c} = {t^6 + e t^4 + c t^3 + s2 t^2 + s1 t + s0},  q = 2^k.  Claimed:
    k odd :  pi(C_{e,c}) = q(q^2 - 1)/6                          for all e, c;
    k even:  pi(C_{e,c}) = q(q - 1)^2/6 + q m(c),
             m(c) = #{nu in K^* non-cube : Tr(c/nu) = 1}   (m(0) = 0,  (q+2)/6 <= m(c) <= q/2 for c != 0).
Exact counts from the Hayes FFT at level 3 (hayes_oppermann.interval_psi), all q^2 classes, q = 2..256.
Usage: python char2_d3_all.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fqlib  # noqa: E402
import hayes_fft_fq as H  # noqa: E402
from hayes_oppermann import interval_psi, prime_power_E  # noqa: E402


def main():
    for k in range(1, 9):
        q = 1 << k
        Psi, ctx = interval_psi(q, 3, 6)
        E = prime_power_E(q, 3, 6, ctx)
        cnt = (np.rint(Psi).astype(np.int64) - E) // 6
        A, M, NEG, INV = fqlib.field_tables(q)
        tr = np.zeros(q, np.int64)
        for x in range(q):
            s, y = 0, x
            for _ in range(k):
                s ^= y
                y = int(M[y, y])
            tr[x] = s
        cubes = {int(M[M[x, x], x]) for x in range(1, q)}
        noncubes = [x for x in range(1, q) if x not in cubes]
        bad = 0
        mvals = set()
        for c in range(q):
            if k % 2:
                pred = q * (q * q - 1) // 6
            else:
                m = sum(1 for nu in noncubes if tr[int(M[c, INV[nu]])] == 1)
                mvals.add(m)
                if c:
                    assert (q + 2) / 6 <= m <= q / 2
                pred = q * (q - 1) ** 2 // 6 + q * m
            for e in range(q):
                bad += int(cnt[H.class_index([1, 0, e, c], ctx)]) != pred
        print(f"q = {q:3d}: {q * q} classes, violations {bad}" + ("" if k % 2 else f", m(c) values {sorted(mvals)}"))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
