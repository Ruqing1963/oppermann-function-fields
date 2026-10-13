"""
swan_parity.py -- the sign obstruction (Stickelberger-Swan) for the zero block Z_d = {t^{2d} + s : deg s < d} over F_2.

In characteristic 2, f^2 has only even-degree terms, so for every d the block Z_d (all coefficients of t^d..t^{2d-1}
equal to 0) is an Oppermann sub-interval: f = t^d + sum_{i<d/2} a_i t^i ... with c chosen to cancel t^d (see notes).

Stickelberger-Swan: if P in Z[t] is monic of degree n and P mod 2 is squarefree with r irreducible factors, then
r = n (mod 2) iff disc(P) = 1 (mod 8).  An irreducible P of even degree n has r = 1, so it needs disc(P) = 5 (mod 8).

For P = t^n + g with 8 | n:  disc(P) = (-1)^{n(n-1)/2} Res(P, P') and P' = n t^{n-1} + g' = g' (mod 8), so
disc(P) = Res(P, g') (mod 8) -- the leading term drops out.  For n = 8 (d = 4) and deg g <= 3 this forces
disc = 1 (mod 8) on every squarefree member: the block contains no irreducible.

This script, for d = 1..dmax: lists every member of Z_d, its factorisation pattern over F_2, the integer-lift
discriminant mod 8 (lift with coefficients 0/1), checks Stickelberger on every squarefree member, and counts the
members whose sign permits irreducibility.
Usage: python swan_parity.py [dmax=8]
"""
import sys

from sympy import Poly, discriminant, factor_list, symbols

t = symbols("t")


def members(d):
    for s in range(1 << d):
        coeffs = [(s >> i) & 1 for i in range(d)]
        yield coeffs, t ** (2 * d) + sum(c * t ** i for i, c in enumerate(coeffs))


def pattern(P):
    _, fl = factor_list(P, modulus=2)
    return sorted(((Poly(f, t).degree(), e) for f, e in fl), reverse=True)


def main():
    dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for d in range(1, dmax + 1):
        n = 2 * d
        n_sqfree = n_sign_ok = n_irr = n_sign_noroot = 0
        bad = 0
        lines = []
        for coeffs, P in members(d):
            pat = pattern(P)
            sqfree = all(e == 1 for _, e in pat)
            r = len(pat)
            irr = sqfree and r == 1
            n_irr += irr
            if sqfree:
                n_sqfree += 1
                D = int(discriminant(P, t)) % 8
                if (D == 1) != (r % 2 == n % 2):
                    bad += 1
                n_sign_ok += (D == 5)
                noroot = coeffs[0] == 1 and (1 + sum(coeffs)) % 2 == 1      # P(0) = P(1) = 1
                n_sign_noroot += (D == 5 and noroot)
            else:
                D = None
            if d <= 4:
                lines.append(f"    {str(P):40s} factors (deg,mult) {pat}  disc mod 8 = {D}")
        print(f"d = {d:2d}  n = {n:2d}  |Z_d| = {1 << d:5d}  squarefree {n_sqfree:5d}  "
              f"disc = 5 mod 8: {n_sign_ok:5d}  (no root in F_2: {n_sign_noroot:5d})  irreducible {n_irr:5d}  Stickelberger violations {bad}")
        for ln in lines:
            print(ln)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
