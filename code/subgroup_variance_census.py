"""
subgroup_variance_census.py -- means and variances of Psi and 2d*pi on the Oppermann subgroup H_d and on the whole
group G_d (Section 5).  Exact class values from the Hayes FFT.

For each (q, d):  |H_d|, E_H[Psi]/q^d, E_H[2d pi]/q^d, Var_H(Psi)/q^d, Var_H(2d pi)/q^d,
                  Var_G(Psi)/q^d, Var_G(2d pi)/q^d, and the Keating--Rudnick value d - 1.
Characteristic 2: q = 2 (d <= 20), 4 (d <= 10), 8 (d <= 7), 16 (d <= 5).  Odd q (global only, H_d = G_d):
q = 3 (d <= 12), 5 (d <= 8), 7 (d <= 7).
Exact checks (Theorem 5.2): Var_H(Psi)/q^2 = q - 1 at d = 2; Var_H(Psi)/q^3 = 1 - 1/q (k odd), 3(1 - 1/q) (k even)
at d = 3; Var_G(Psi)/q^2 = 1 - 1/q at d = 2.
Output: ../data/subgroup_variance.txt and ../data/subgroup_variance.json
Usage: python subgroup_variance_census.py [work_dir]
"""
import json
import os
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from subgroup_lib import class_data, subgroup_indices  # noqa: E402

RANGES = {2: 20, 4: 10, 8: 7, 16: 5, 3: 12, 5: 8, 7: 7}
HERE = os.path.dirname(os.path.abspath(__file__))


def mv(v):
    v = v.astype(np.float64)
    m = v.mean()
    return m, ((v - m) ** 2).mean()


def exact_var(v):
    vals = [int(x) for x in v]
    m = Fraction(sum(vals), len(vals))
    return sum((Fraction(x) - m) ** 2 for x in vals) / len(vals)


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else None
    rows = []
    out = open(os.path.join(HERE, "..", "data", "subgroup_variance.txt"), "w", encoding="utf-8")
    hdr = (f"{'q':>3} {'d':>3} {'|H|':>7} {'E_H Psi/q^d':>12} {'E_H 2dpi/q^d':>13} {'VarH Psi':>10} {'VarH 2dpi':>10} "
           f"{'VarG Psi':>10} {'VarG 2dpi':>10} {'KR d-1':>7} {'exact check':>12}")
    print(hdr)
    out.write(hdr + "\n")
    for q, dmax in RANGES.items():
        for d in range(2, dmax + 1):
            Psi, pi, E, ctx, resid = class_data(q, d, work)
            qd = float(q) ** d
            n = 2 * d
            mG, vG = mv(Psi)
            _, vGpi = mv(n * pi)
            if q % 2 == 0:
                Hidx = subgroup_indices(q, d, ctx)
            else:
                Hidx = np.arange(q ** d)
            mH, vH = mv(Psi[Hidx])
            mHpi, vHpi = mv(n * pi[Hidx])
            check = "-"
            k = q.bit_length() - 1
            if q % 2 == 0 and d == 2:
                check = str(exact_var(Psi[Hidx]) == Fraction(q * q * (q - 1))
                            and exact_var(Psi) == Fraction(q * (q - 1)))
            if q % 2 == 0 and d == 3:
                pred = Fraction(q * q * (q - 1)) * (1 if k % 2 else 3)
                check = str(exact_var(Psi[Hidx]) == pred)
            r = dict(q=q, d=d, H=len(Hidx), meanH=mH / qd, meanHpi=mHpi / qd, varH=vH / qd, varHpi=vHpi / qd,
                     varG=vG / qd, varGpi=vGpi / qd, KR=d - 1, check=check, resid=resid)
            rows.append(r)
            line = (f"{q:>3} {d:>3} {len(Hidx):>7} {r['meanH']:>12.6f} {r['meanHpi']:>13.6f} {r['varH']:>10.4f} "
                    f"{r['varHpi']:>10.4f} {r['varG']:>10.4f} {r['varGpi']:>10.4f} {d - 1:>7} {check:>12}")
            print(line)
            out.write(line + "\n")
            sys.stdout.flush()
            out.flush()
    out.close()
    with open(os.path.join(HERE, "..", "data", "subgroup_variance.json"), "w") as fh:
        json.dump(rows, fh, indent=1)


if __name__ == "__main__":
    main()
