"""
verify_subgroup_mean.py -- Theorems 2.1, 2.2, 3.1, 3.2 (Oppermann subgroup H_d, characteristic 2).

For q in {2, 4, 8, 16} and d in {2, 3, 4, 5}:
  Thm 2.1  (a) e_j = 0 for all odd j < d  <=>  S_j = 0 for all odd j < d   (every class);
           (b) H_d is closed under the group law of G_d (checked on discrete logarithms).
  Thm 2.2  psi_{2d}(chi_beta) = sum_{x in F_{q^{2d}}} (-1)^{Tr(g_beta(x))}  for every beta in F_2^{floor(d/2)}
           (left side from the exact class values Psi, right side by brute force in F_{2^{2dk}}; when q^{2d} <= 2^16).
  Thm 3.1/3.2  E_{H_d}[Psi] and E_{H_d}[pi] against the closed forms
           d = 2: q^2,                    q(q-1)/4
           d = 3: q^3,                    q(q^2-1)/6
           d = 4: q^4 (k even),           (q^4-q^2)/8           ;  q^4-2q^2+2q (k odd),  (q^4-3q^2+2q)/8
           d = 5: q^5+q^2-q,              (q^5-q)/10
Usage: python verify_subgroup_mean.py  > ../data/verify_subgroup_mean.txt
"""
import os
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fqlib  # noqa: E402
from subgroup_lib import (GF2N, class_data, power_sums, prefix_of_index, subgroup_indices,  # noqa: E402
                          trace_table)


def closed_form(q, d):
    k = q.bit_length() - 1
    if d == 2:
        return Fraction(q ** 2), Fraction(q * (q - 1), 4)
    if d == 3:
        return Fraction(q ** 3), Fraction(q * (q * q - 1), 6)
    if d == 4:
        if k % 2 == 0:
            return Fraction(q ** 4), Fraction(q ** 4 - q ** 2, 8)
        return Fraction(q ** 4 - 2 * q ** 2 + 2 * q), Fraction(q ** 4 - 3 * q ** 2 + 2 * q, 8)
    if d == 5:
        return Fraction(q ** 5 + q ** 2 - q), Fraction(q ** 5 - q, 10)
    raise ValueError


def main():
    print(f"{'q':>3} {'d':>2} {'|H|':>6} {'E_H[Psi]':>14} {'formula':>14} {'E_H[pi]':>14} {'formula':>14} "
          f"{'thm2.1a':>8} {'thm2.1b':>8} {'thm2.2 (#beta)':>15}")
    for q in (2, 4, 8, 16):
        A, M, NEG, INV = fqlib.field_tables(q)
        tr = trace_table(q, A, M)
        k = q.bit_length() - 1
        for d in (2, 3, 4, 5):
            Psi, pi, E, ctx, resid = class_data(q, d)
            Hidx = subgroup_indices(q, d, ctx)
            mPsi = Fraction(int(Psi[Hidx].sum()), len(Hidx))
            mpi = Fraction(int(pi[Hidx].sum()), len(Hidx))
            fPsi, fpi = closed_form(q, d)
            # Theorem 2.1 (a)
            P = prefix_of_index(q, d, ctx)
            odd = [j for j in range(1, d) if j % 2]
            inH = np.zeros(q ** d, bool)
            inH[Hidx] = True
            ok_a = True
            Sodd = np.zeros((q ** d, len(odd)), np.int64)
            for idx in range(q ** d):
                e = [int(v) for v in P[idx]]
                S = power_sums(e, q, A, M)
                Sodd[idx] = [S[j - 1] for j in odd]
                ok_a &= (all(e[j - 1] == 0 for j in odd) == all(S[j - 1] == 0 for j in odd)) and \
                    (all(e[j - 1] == 0 for j in odd) == bool(inH[idx]))
            # Theorem 2.1 (b): closure on discrete logs
            strides, orders = ctx["strides"], ctx["orders"]
            def exps(i):
                return tuple((i // strides[a]) % orders[a] for a in range(len(strides)))
            Hset = {exps(int(i)) for i in Hidx}
            ok_b = all(tuple((x[a] + y[a]) % orders[a] for a in range(len(strides))) in Hset
                       for x in Hset for y in Hset) if len(Hset) <= 512 else None
            # Theorem 2.2
            n_beta = "-"
            Nbits = 2 * d * k
            if Nbits <= 16:
                F = GF2N(Nbits)
                good = 0
                for bmask in range(1, 1 << len(odd)):
                    beta = [odd[i] for i in range(len(odd)) if bmask >> i & 1]
                    sign = np.array([1 - 2 * tr[np.bitwise_xor.reduce(Sodd[idx, [odd.index(j) for j in beta]])]
                                     for idx in range(q ** d)])
                    lhs = int((Psi * sign).sum())
                    rhs = 0
                    for x in range(1 << Nbits):
                        s = 0
                        for j in beta:
                            s ^= F.pow(x, j)
                        rhs += 1 - 2 * F.trace(s)
                    good += lhs == rhs
                n_beta = f"{good}/{(1 << len(odd)) - 1}"
            print(f"{q:>3} {d:>2} {len(Hidx):>6} {str(mPsi):>14} {str(fPsi):>14} {str(mpi):>14} {str(fpi):>14} "
                  f"{str(ok_a):>8} {str(ok_b):>8} {n_beta:>15}"
                  + ("" if (mPsi == fPsi and mpi == fpi) else "   <-- MISMATCH"))
            sys.stdout.flush()


if __name__ == "__main__":
    main()
