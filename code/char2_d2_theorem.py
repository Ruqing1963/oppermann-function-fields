"""
char2_d2_theorem.py -- brute-force check of the exact d = 2 count in characteristic 2.

Theorem (notes, LaTeX Thm 3.3).  q = 2^k, B_e = {t^4 + e t^2 + a t + b : a, b in F_q}.
  (i)  for a != 0, the cubic y^3 + e y + a has exactly one root in F_q  iff  Tr(e^3 / a^2) != Tr(1);
  (ii) for such (e, a), exactly q/2 values of b make t^4 + e t^2 + a t + b irreducible; for the other a != 0, none;
  (iii) hence #irreducibles in B_e = q^2/4 (k even, e != 0), 0 (k even, e = 0), q(q-2)/4 (k odd, e != 0),
        q(q-1)/2 (k odd, e = 0).
Usage: python char2_d2_theorem.py [kmax=7]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from char2_d2_blocks import is_irreducible_quartic, roots_in_field, tables, trace  # noqa: E402


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for k in range(1, kmax + 1):
        q, M = tables(k)
        inv = [0] * q
        for x in range(1, q):
            inv[x] = next(y for y in range(1, q) if M[x, y] == 1)
        tr1 = trace(1, k, M)
        bad_i = bad_ii = bad_iii = 0
        for e in range(q):
            e3 = int(M[M[e, e], e])
            total = 0
            for a in range(1, q):
                one_root = roots_in_field([a, e, 0, 1], q, M) == 1
                crit = trace(int(M[e3, M[inv[a], inv[a]]]), k, M) != tr1
                bad_i += one_root != crit
                nb = sum(1 for b in range(q) if is_irreducible_quartic(e, a, b, q, M))
                bad_ii += nb != (q // 2 if one_root else 0)
                total += nb
            if k % 2 == 0:
                pred = 0 if e == 0 else q * q // 4
            else:
                pred = q * (q - 1) // 2 if e == 0 else q * (q - 2) // 4
            bad_iii += total != pred
        print(f"q = {q:4d}: violations of (i) {bad_i}, (ii) {bad_ii}, (iii) {bad_iii}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
