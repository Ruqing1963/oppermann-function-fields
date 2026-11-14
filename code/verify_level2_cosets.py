"""
verify_level2_cosets.py -- the level-2 part of the coset sums A(lambda) on the Oppermann subgroup (characteristic 2).

For a != 0 let lambda_a(x) = (-1)^{Tr(a e_2(x))} on H_d.  The q characters chi of G_d of level <= 2 that restrict to
lambda_a have, by the level-2 rigidity gamma^4 = (-1)^k q^2,
    psi_{2d}(chi) = -(-1)^{km} q^d                 (d = 2m),
    psi_{2d}(chi) = (-1)^{km} q^{d-1} psi_2(chi)   (d = 2m+1),   sum over the q characters of psi_2 = 0,
so their total is  P(a) = -(-1)^{km} q^{d+1}  (d even),  0  (d odd).
From the exact class values:  P(a) = q * sum_{x in G_d, e_1(x) = 0} Psi(x) (-1)^{Tr(a e_2(x))},  and the full coset sum
is  A(lambda_a) = q^{floor(d/2)} * sum_{x in H_d} Psi(x) lambda_a(x).  The script checks P(a) for every a != 0 and
reports the share of the level-2 parts in the subgroup variance,
    share = sum_a P(a)^2 / sum_{lambda != 1} A(lambda)^2 .
Usage: python verify_level2_cosets.py > ../data/level2_cosets.txt
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fqlib  # noqa: E402
from subgroup_lib import class_data, prefix_of_index, subgroup_indices, trace_table  # noqa: E402


def main():
    print(f"{'q':>2} {'d':>2} {'P(a), all a != 0':>18} {'predicted':>12} {'ok':>5} {'A(lambda_a): min .. max':>28} "
          f"{'level-2 share of Var_H':>23}")
    for q, dmax in ((2, 10), (4, 6), (8, 4), (16, 4)):
        k = q.bit_length() - 1
        A_, M, NEG, INV = fqlib.field_tables(q)
        tr = trace_table(q, A_, M)
        for d in range(2, dmax + 1):
            Psi, pi, E, ctx, _ = class_data(q, d)
            P = prefix_of_index(q, d, ctx)
            e1, e2 = P[:, 0], P[:, 1]
            Hidx = subgroup_indices(q, d, ctx)
            m = d // 2
            pred = -((-1) ** (k * m)) * q ** (d + 1) if d % 2 == 0 else 0
            vals, Avals = [], []
            sel = e1 == 0
            for a in range(1, q):
                sign = 1 - 2 * tr[M[a, e2]]
                vals.append(q * int(np.sum(Psi[sel] * sign[sel])))
                Avals.append(q ** (d // 2) * int(np.sum(Psi[Hidx] * sign[Hidx])))
            ok = all(v == pred for v in vals)
            # sum over nontrivial lambda of A(lambda)^2 = q^{2d} Var_H(Psi)  (Proposition 7.4)
            x = Psi[Hidx].astype(np.float64)
            varH = float(np.mean(x * x) - np.mean(x) ** 2)
            share = sum(v * v for v in vals) / (q ** (2 * d) * varH) if varH > 0 else float("nan")
            print(f"{q:>2} {d:>2} {str(sorted(set(vals))):>18} {pred:>12} {str(ok):>5} "
                  f"{str(min(Avals)) + ' .. ' + str(max(Avals)):>28} {share:>23.6f}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
