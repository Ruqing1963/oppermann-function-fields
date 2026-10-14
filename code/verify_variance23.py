"""
verify_variance23.py -- exact check of Theorem 5.2 (q = 2, 4, ..., 64):
  d = 2:  Var_H(Psi)/q^2 = Var_H(4 pi)/q^2 = q - 1,   Var_G(Psi)/q^2 = 1 - 1/q;
  d = 3:  Var_H(Psi)/q^3 = 1 - 1/q (k odd), 3(1 - 1/q) (k even);
          Var_H(6 pi)/q^3 = 0 (k odd), 6(1 - 1/q) (k even).
Usage: python verify_variance23.py > ../data/verify_variance23.txt
"""
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from subgroup_lib import class_data, subgroup_indices  # noqa: E402


def var(vals):
    vals = [int(v) for v in vals]
    m = Fraction(sum(vals), len(vals))
    return sum((Fraction(v) - m) ** 2 for v in vals) / len(vals)


def main():
    for k in range(1, 7):
        q = 1 << k
        P2, pi2, _, c2, _ = class_data(q, 2)
        H2 = subgroup_indices(q, 2, c2)
        ok2 = (var(P2[H2]) / q ** 2 == q - 1 and var(4 * pi2[H2]) / q ** 2 == q - 1
               and var(P2) / q ** 2 == 1 - Fraction(1, q))
        P3, pi3, _, c3, _ = class_data(q, 3)
        H3 = subgroup_indices(q, 3, c3)
        f = 1 if k % 2 else 3
        ok3 = (var(P3[H3]) / q ** 3 == f * (1 - Fraction(1, q))
               and var(6 * pi3[H3]) / q ** 3 == (0 if k % 2 else 6 * (1 - Fraction(1, q))))
        print(f"q = {q:2d}: d = 2 {ok2},  d = 3 {ok3}   "
              f"(VarH Psi/q^3 = {var(P3[H3]) / q ** 3}, VarH 6pi/q^3 = {var(6 * pi3[H3]) / q ** 3})")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
