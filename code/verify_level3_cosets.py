"""
verify_level3_cosets.py -- Proposition 5.6: coset sums of psi_6 over the annihilator of H_3 (characteristic 2).

For lambda_{a,b}(e, c) = (-1)^{Tr(a e + b c)} on H_3 = {1 + e u^2 + c u^3} and
A(lambda) = sum of psi_6(chi) over the q characters chi of G_3 restricting to lambda:
    A(lambda_{a,b}) = 0                      (a != 0),
    A(lambda_{0,b}) = q^4                    (b != 0, k odd; or k even and b a cube),
    A(lambda_{0,b}) = -2 q^4                 (k even, b a non-cube).
A(lambda) = q * sum_{(e,c) in H_3} Psi(e,c) lambda(e,c) is computed from the exact class values.
Consequence checked with the L-functions (Hayes FFT): for k even at least 2q(q-1)/3 characters of exact level 3
have psi_6 = -2q^3, i.e. both inverse zeros satisfy gamma^6 = q^3.
Usage: python verify_level3_cosets.py > ../data/level3_cosets.txt
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fqlib  # noqa: E402
import hayes_fft_fq as H  # noqa: E402
from subgroup_lib import class_data, trace_table  # noqa: E402


def main():
    for k in range(1, 7):
        q = 1 << k
        A_, M, NEG, INV = fqlib.field_tables(q)
        tr = trace_table(q, A_, M)
        cubes = {int(M[M[x, x], x]) for x in range(1, q)}
        Psi, pi, E, ctx, _ = class_data(q, 3)
        val = {}
        for e in range(q):
            for c in range(q):
                val[(e, c)] = int(Psi[H.class_index([1, 0, e, c], ctx)])
        bad = 0
        for a in range(q):
            for b in range(q):
                if a == 0 and b == 0:
                    continue
                s = 0
                for (e, c), v in val.items():
                    s += v * (1 - 2 * int(tr[int(A_[M[a, e], M[b, c]])]))
                Aval = q * s
                if a != 0:
                    pred = 0
                elif k % 2 == 1 or b in cubes:
                    pred = q ** 4
                else:
                    pred = -2 * q ** 4
                bad += Aval != pred
        # L-function check: characters with psi_6 = -2 q^3
        extra = ""
        if k % 2 == 0:
            l, N = 3, q ** 3
            C = []
            for kdeg in range(l):
                ind = H._indices_degree(kdeg, l, q, ctx["p"], ctx["k"], ctx["A"], ctx["M"], ctx["NEG"], ctx["FROB"],
                                        ctx["FINV"], ctx["axis_of"], ctx["orders"], ctx["strides"])
                C.append(np.fft.fftn(np.bincount(ind, minlength=N).astype(complex).reshape(ctx["shape"])).reshape(N))
            c1, c2 = C[1], C[2]
            p = [None, -c1, c1 * c1 - 2 * c2]
            for n in range(3, 7):
                p.append(-c1 * p[n - 1] - c2 * p[n - 2])
            psi6 = -p[6]
            n_ext = int(np.sum(np.abs(psi6[1:] + 2 * q ** 3) < 1e-6 * q ** 3))
            extra = f";  characters with psi_6 = -2q^3: {n_ext} (>= 2q(q-1)/3 = {2 * q * (q - 1) // 3})"
        print(f"q = {q:2d}: coset-sum violations {bad} of {q * q - 1}{extra}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
