"""
cubic_gauss_m.py -- Theorem 4.1: the two values of m(c) (d = 3, characteristic 2, k even).

m(c) = #{nu in F_q^* non-cube : Tr(c/nu) = 1}.  Claimed, with eps = (-1)^{k/2}:
    m(c) = (q - eps sqrt q)/3      if c is a nonzero cube,
    m(c) = (2q + eps sqrt q)/6     if c is a non-cube,
from the cubic Gauss sum G(omega) = sum_y omega(y)(-1)^{Tr y} = -eps sqrt q.
For q = 4, 16, ..., 2^16 the script computes G(omega) and m(c) for every c (q <= 1024) or for c = 1 and one
non-cube c (larger q), and compares.  Also prints the resulting class counts pi(C_{e,c}).
Usage: python cubic_gauss_m.py [kmax=16]
"""
import cmath
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from subgroup_lib import GF2N  # noqa: E402


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    print(f"{'q':>6} {'eps':>4} {'G(omega)':>22} {'-eps*sqrt q':>12} {'m(cube)':>8} {'pred':>6} {'m(non-cube)':>11} "
          f"{'pred':>6} {'all c':>6} {'pi(C_e,0)':>12} {'pi(C_e,cube)':>13} {'pi(C_e,non)':>12}")
    for k in range(2, kmax + 1, 2):
        q = 1 << k
        F = GF2N(k)
        sq = math.isqrt(q)
        eps = (-1) ** (k // 2)
        # generator of F_q^*: element of order q - 1
        fac = [p for p in range(2, q) if (q - 1) % p == 0 and all(p % r for r in range(2, int(p ** 0.5) + 1))]
        g = next(x for x in range(2, q) if all(F.pow(x, (q - 1) // p) != 1 for p in fac))
        log = [None] * q
        y = 1
        for i in range(q - 1):
            log[y] = i
            y = F.mul(y, g)
        tr = [F.trace(x) for x in range(q)]
        zeta = cmath.exp(2j * cmath.pi / 3)
        G = sum(zeta ** (log[x] % 3) * (1 - 2 * tr[x]) for x in range(1, q))
        inv = [0] * q
        for x in range(1, q):
            inv[x] = F.pow(g, (q - 1 - log[x]) % (q - 1))

        def m(c):
            return sum(1 for nu in range(1, q) if log[nu] % 3 and tr[F.mul(c, inv[nu])] == 1)

        pred_cube = (q - eps * sq) // 3
        pred_non = (2 * q + eps * sq) // 6
        noncube = g
        m_cube, m_non = m(1), m(noncube)
        allc = "-"
        if q <= 1024:
            allc = str(all(m(c) == (pred_cube if log[c] % 3 == 0 else pred_non) for c in range(1, q)))
        base = q * (q - 1) ** 2 // 6
        print(f"{q:>6} {eps:>4} {G.real:>11.4f}{G.imag:+.4f}i {-eps * sq:>12} {m_cube:>8} {pred_cube:>6} "
              f"{m_non:>11} {pred_non:>6} {allc:>6} {base:>12} {base + q * pred_cube:>13} {base + q * pred_non:>12}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
