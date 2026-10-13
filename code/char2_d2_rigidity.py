"""
char2_d2_rigidity.py -- level-2 Hayes characters in degree 4, characteristic 2 (paper, Proposition 6.2).

Checks, for q = 2, 4, 8, 16:
  (1) E(B_e) = q for every e (proper prime powers in t^4 + e t^2 + a t + b);
  (2) every character of exact level 2 has L(u, chi) = 1 - gamma u with gamma^4 = (-1)^k q^2;
  (3) Psi(x) = q^2 - (-1)^k (q^2 [x = 1] - q [p3(x) = 0]) for every level-2 class x;
  (4) every class with p3 != 0 contains exactly q^2/4 irreducibles.
Usage: python char2_d2_rigidity.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hayes_fft_fq as H  # noqa: E402
from hayes_oppermann import interval_psi, prime_power_E  # noqa: E402


def main():
    for k in (1, 2, 3, 4):
        q = 1 << k
        sgn = (-1) ** k
        l, n = 2, 4
        Psi, ctx = interval_psi(q, l, n)
        E = prime_power_E(q, l, n, ctx)
        cnt = (np.rint(Psi).astype(np.int64) - E) // n
        # class index of the prefix (p3, p2): series [1, p3, p2]
        bad1 = bad3 = bad4 = 0
        for p3 in range(q):
            for p2 in range(q):
                x = H.class_index([1, p3, p2], ctx)
                pred = q * q - sgn * ((q * q if (p3 == 0 and p2 == 0) else 0) - (q if p3 == 0 else 0))
                bad3 += abs(Psi[x] - pred) > 1e-6
                if p3 == 0:
                    bad1 += E[x] != q
                else:
                    bad4 += cnt[x] != q * q // 4
        # L-functions: c1(chi) = FFT of the class distribution of monic linears t + a  (series [1, a, 0])
        N = q ** l
        ind = H._indices_degree(1, l, q, ctx["p"], ctx["k"], ctx["A"], ctx["M"], ctx["NEG"], ctx["FROB"], ctx["FINV"],
                                ctx["axis_of"], ctx["orders"], ctx["strides"])
        c1 = np.fft.fftn(np.bincount(ind, minlength=N).astype(complex).reshape(ctx["shape"])).reshape(N)
        # characters of exact level <= 1 have c1 = 0 (or q for the trivial one); the others should be exact level 2
        g4 = [(-c1[i]) ** 4 for i in range(1, N) if abs(c1[i]) > 1e-9]
        n_nonzero = len(g4)
        bad2 = sum(abs(g - sgn * q * q) > 1e-6 for g in g4)
        print(f"q = {q:2d}: E(B_e)=q violations {bad1}; characters with deg L = 1: {n_nonzero} (expected q^2-q = "
              f"{q * q - q}); gamma^4 != (-1)^k q^2: {bad2}; Psi formula violations {bad3}; "
              f"p3!=0 classes with pi != q^2/4: {bad4}")


if __name__ == "__main__":
    main()
