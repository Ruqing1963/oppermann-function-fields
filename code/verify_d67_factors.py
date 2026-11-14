"""
verify_d67_factors.py -- the polynomial facts used in the evaluation of the mean drift for d = 6, 7.

g_{A,B}(X) = X^10 + A X^6 + B X^5 + B^2  over K = F_q, q = 2^k, B != 0.

Symbolic part (SymPy, coefficients mod 2):
  (a) d = 7: if g = M * C with M = X^7 + m5 X^5 + ... + m0 (no X^6 term) and C a monic cubic, then
      M = X^7 + s X^5 + s t X^4 + (t^4 + s t^2) X^3 + s t^3 X^2 + s t^4 X + s t^5 and M(t) = 0;
  (b) d = 6: if g = M * Q with M = X^6 + m4 X^4 + ... + m0 (no X^5 term) and Q a monic quartic, then with
      m0 = h^2, m1 = B = p:  M = (X^2 + (p^2/h^3) X + p/h)(X^4 + (p^2/h^3) X^3 + h^3/p);
  (c) every cubic X^3 + c1 X + c0 with c0 c1 != 0 divides g_{A,B} for B = c0^3/c1^2, A = c0^4/c1^4 + c0^2/c1 + c1^2;
  (d) g_{B+1,B} = (X^4 + X^3 + B)((X^2 + X)^3 + B).

Census part: for q = 2^k, k = 1..6, all pairs ((A, B), C) with B != 0 and C = X^3 + X^2 + c1 X + c0 an irreducible
cubic over F_q dividing g_{A,B}.  Expected: none for k odd; for k even exactly one for every non-cube B, with
A = B + 1 (the factor of (X^2+X)^3 + B with trace 1), i.e. 2(q-1)/3 pairs.  Every pair found is confirmed by an
explicit polynomial division.  Also counts I = #{c : V^3 + V + c irreducible} = (q+1)/3 (k odd), (q-1)/3 (k even).

Usage: python verify_d67_factors.py [kmax]     (default 6; output also in ../data/verify_d67_factors.txt)
"""
import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = open(os.path.join(HERE, "..", "data", "verify_d67_factors.txt"), "w", encoding="utf8")
IRRED = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011}


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")
    OUT.flush()


def symbolic():
    X, s, t, p, h, c0, c1, B = sp.symbols("X s t p h c0 c1 B")
    zero = lambda ex, *g: sp.Poly(sp.expand(ex), *g, modulus=2).is_zero
    say("=== symbolic identities (mod 2) ===")
    # (a) the coefficient equations of g = M * C for d = 7, solved in the text; here: the solution satisfies them
    M7 = X**7 + s*X**5 + s*t*X**4 + (t**4 + s*t**2)*X**3 + s*t**3*X**2 + s*t**4*X + s*t**5
    C3 = X**3 + s*X + s*t                                   # c2 = 0, c1 = m5, c0 = m4
    Bv = s*t**3
    Av = (t**4 + s*t**2) + s**2                             # m3 + m5^2
    say("  (a) M7 * (X^3 + m5 X + m4) = g_{A,B} with B = m2, A = m3 + m5^2: %s;  M7(t) = 0: %s" %
        (zero(M7*C3 - (X**10 + Av*X**6 + Bv*X**5 + Bv**2), X, s, t), zero(M7.subs(X, t), s, t)))
    # (b) d = 6; multiply by powers of h and p to clear denominators
    al, be, ga = p**2/h**3, p/h, h**3/p
    M6 = X**6 + (p/h + p**4/h**6)*X**4 + (p**3/h**4)*X**3 + (h**3/p)*X**2 + p*X + h**2
    Q4 = X**4 + (p/h + p**4/h**6)*X**2 + (p**3/h**4)*X + p**2/h**2     # q3 = 0, q2 = m4, q1 = m3, q0 = p^2/m0
    A6 = p**2/h**2 + h**3/p + (p/h + p**4/h**6)**2                       # q0 + m2 + m4^2
    e1 = sp.together(sp.expand(M6 - (X**2 + al*X + be)*(X**4 + al*X**3 + ga)))
    e2 = sp.together(sp.expand(M6*Q4 - (X**10 + A6*X**6 + p*X**5 + p**2)))
    say("  (b) M6 = (X^2 + (p^2/h^3) X + p/h)(X^4 + (p^2/h^3) X^3 + h^3/p): %s;  M6 * Q4 = g_{A,B}: %s" %
        (zero(sp.numer(e1), X, p, h), zero(sp.numer(e2), X, p, h)))
    # (c)
    Bc = c0**3/c1**2
    Ac = c0**4/c1**4 + c0**2/c1 + c1**2
    D = X**7 + c1*X**5 + c0*X**4 + (Ac + c1**2)*X**3 + Bc*X**2 + (c0**2 + c1*Ac + c1**3)*X + (c0*Ac + c0*c1**2 + c1*Bc)
    e3 = sp.together(sp.expand((X**3 + c1*X + c0)*D - (X**10 + Ac*X**6 + Bc*X**5 + Bc**2)))
    say("  (c) (X^3 + c1 X + c0) divides g_{A,B}, B = c0^3/c1^2, A = c0^4/c1^4 + c0^2/c1 + c1^2: %s" % zero(sp.numer(e3), X, c0, c1))
    # (d)
    say("  (d) g_{B+1,B} = (X^4 + X^3 + B)((X^2 + X)^3 + B): %s" %
        zero((X**4 + X**3 + B)*((X**2 + X)**3 + B) - (X**10 + (B + 1)*X**6 + B*X**5 + B**2), X, B))


def resultant_proof():
    """The identities of the proof of Theorem 7.11 (no irreducible cubic X^3+X^2+c1X+c0 divides g_{A,B} outside the
    family c0 = c1^2, A = B+1, B = c1^4+c1^3), all in F_2[c1, c0] (and F_4[c1, c0] for (j))."""
    X, c0, c1, om, u = sp.symbols("X c0 c1 om u")
    P = lambda ex, *g: sp.Poly(sp.expand(ex), *g, modulus=2)
    zero = lambda ex, *g: P(ex, *g).is_zero
    C = X**3 + X**2 + c1*X + c0
    red = lambda ex: sp.rem(P(ex, X, c1, c0), P(C, X, c1, c0)).as_expr()
    say("=== Theorem 7.11: elimination for cubic factors X^3 + X^2 + c1 X + c0 of g_{A,B} (mod 2) ===")
    w5, w6, w10 = red(X**5), red(X**6), red(X**10)
    e5 = (c0 + 1)*X**2 + (c1**2 + c1 + c0)*X + (c1 + 1)*c0
    e6 = (c1**2 + c1 + 1)*X**2 + (c1 + c0)*X + (c0**2 + c0)
    e10 = (c1**4 + c1**2 + c1*c0**2 + c1 + 1)*X**2 + (c0 + 1)**2*(c1 + c0)*X + (c1**2*c0**2 + c0**3 + c0**2 + c0)
    say("  (e) X^5, X^6, X^10 mod C as in the text: %s %s %s" % (zero(w5 - e5, X, c1, c0), zero(w6 - e6, X, c1, c0), zero(w10 - e10, X, c1, c0)))
    a2, b2, r2 = c1**2 + c1 + 1, c0 + 1, c1**4 + c1**2 + c1*c0**2 + c1 + 1
    a1, b1, r1 = c1 + c0, c1**2 + c1 + c0, (c0 + 1)**2*(c1 + c0)
    a0, b0, r0 = c0**2 + c0, (c1 + 1)*c0, c1**2*c0**2 + c0**3 + c0**2 + c0
    Delta = c1**4 + c1**2*c0 + c0**2
    NB = (c1 + c0)*(c1**2 + c1*c0 + c0)**2
    NA = (c1**6 + c1**5 + c1**4 + c0*c1**4 + c0*c1**2 + c0**2*(c1**3 + c1**2 + c1 + 1) + c0**3 + c0**4)
    say("  (f) det = Delta: %s;  Cramer: Delta*B = N_B: %s, Delta*A = N_A: %s" % (
        zero(a2*b1 - b2*a1 - Delta, c1, c0), zero(a2*r1 - a1*r2 - NB, c1, c0), zero(r2*b1 - r1*b2 - NA, c1, c0)))
    # Delta = 0: N_B = 0, so c1 = c0 (then C(1) = 0) or c1^2 + c1 c0 + c0 = 0; then (c1+1)^2 Delta = c1^4 (c1^2+c1+1),
    # so c1^2 + c1 + 1 = 0, c0 = 1 and C(c1) = c1^3 + 1 = 0
    say("  (g) C(1) = c1 + c0: %s;  with c0 = c1^2/(c1+1): (c1+1)^2 Delta = c1^4 (c1^2+c1+1): %s;  C(c1) at c0 = 1 equals c1^3 + 1: %s" % (
        zero(C.subs(X, 1) - (c1 + c0), c1, c0),
        zero((c1 + 1)**2*c1**4 + c1**2*c1**2*(c1 + 1) + c1**4 - c1**4*(c1**2 + c1 + 1), c1),
        zero(C.subs({X: c1, c0: 1}) - (c1**3 + 1), c1)))
    E = NB**2 + b0*NB*Delta + a0*NA*Delta + r0*Delta**2
    F = (c1**2 + c0)*(c1**4 + c1**2*c0 + c0**3)*((c1**2 + c0)**2 + c0**2*(c1**2 + c0) + c0**4)
    say("  (h) E := Delta^2 (B^2 + (c1+1)c0 B + (c0^2+c0) A + c1^2c0^2+c0^3+c0^2+c0) at A = N_A/Delta, B = N_B/Delta equals "
        "(c1^2+c0)(c1^4+c1^2c0+c0^3)((c1^2+c0)^2+c0^2(c1^2+c0)+c0^4): %s" % zero(E - F, c1, c0))
    say("  (i) c1^3 C((c0+c1)/c1) = c1^4 + c1^2 c0 + c0^3: %s" % zero(sp.together(c1**3*C.subs(X, (c0 + c1)/c1)) - (c1**4 + c1**2*c0 + c0**3), c1, c0))
    # (j) over F_4: om^2 + om + 1 = 0, u^2 = c0;  C(om^2 u)^2 = om c0 (c1^2 + om^2 c0^2 + c0)
    lhs = sp.expand(C.subs(X, om**2*u)**2).subs(c0, u**2)
    rhs = sp.expand(om*u**2*(c1**2 + om**2*u**4 + u**2))
    dj = sp.rem(P(lhs - rhs, om, u, c1), P(om**2 + om + 1, om, u, c1))
    pr = sp.rem(P((c1**2 + om**2*c0**2 + c0)*(c1**2 + om*c0**2 + c0) - ((c1**2 + c0)**2 + c0**2*(c1**2 + c0) + c0**4), om, c1, c0), P(om**2 + om + 1, om, c1, c0))
    say("  (j) over F_4: third factor = prod_om (c1^2 + om^2 c0^2 + c0): %s;  C(om^2 sqrt(c0))^2 = om c0 (c1^2 + om^2 c0^2 + c0): %s" % (pr.is_zero, dj.is_zero))
    # (k) the surviving branch c0 = c1^2
    sub = {c0: c1**2}
    Bk = c1**4 + c1**3
    say("  (k) c0 = c1^2: Delta = c1^4: %s;  N_B = c1^4 (c1^4+c1^3): %s;  N_A = c1^4 (c1^4+c1^3+1): %s;  (X^2+X)^3 = c1^4+c1^3 mod C: %s" % (
        zero(Delta.subs(sub) - c1**4, c1), zero(NB.subs(sub) - c1**4*Bk, c1), zero(NA.subs(sub) - c1**4*(Bk + 1), c1),
        zero(sp.rem(P(((X**2 + X)**3 - Bk), X, c1), P(C.subs(sub), X, c1)).as_expr(), X, c1)))


def tables(k):
    q = 1 << k
    mod = IRRED[k]

    def mul(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            a <<= 1
            if a & q:
                a ^= mod
            b >>= 1
        return r
    for g in range(1, q):
        x, n = 1, 0
        while True:
            x = mul(x, g)
            n += 1
            if x == 1:
                break
        if n == q - 1:
            break
    exp = np.zeros(2 * q, np.int64)
    log = np.zeros(q, np.int64)
    x = 1
    for i in range(max(q - 1, 1)):
        exp[i] = x
        log[x] = i
        x = mul(x, g)
    exp[q - 1:2 * q - 2] = exp[:q - 1]

    def vmul(a, b):
        a = np.asarray(a)
        b = np.asarray(b)
        return np.where((a == 0) | (b == 0), 0, exp[(log[a] + log[b]) % max(q - 1, 1)])
    return q, mul, vmul


def divides(mul, C, g):
    g = g[:]
    for i in range(len(g) - 1, 2, -1):
        c = g[i]
        if c:
            for j in range(4):
                g[i - 3 + j] ^= mul(c, C[j])
    return not any(g[:3])


def census(k):
    q, mul, vmul = tables(k)
    xs = np.arange(q)
    cube = lambda v: vmul(v, vmul(v, v))
    # I = number of c with V^3 + V + c irreducible
    vals = set(int(v) for v in (cube(xs) ^ xs))
    I = q - len(vals)
    c0, c1 = np.meshgrid(xs, xs, indexing="ij")
    irr = np.ones((q, q), bool)
    for x in range(q):
        irr &= ((mul(mul(x, x), x) ^ mul(x, x)) ^ vmul(c1, x) ^ c0) != 0
    c0sq, c1sq = vmul(c0, c0), vmul(c1, c1)
    a2 = c1sq ^ c1 ^ 1
    b2 = c0 ^ 1
    r2 = vmul(c0sq, c1) ^ vmul(c1sq, c1sq) ^ c1sq ^ c1 ^ 1
    a1 = c0 ^ c1
    b1 = c0 ^ c1sq ^ c1
    r1 = vmul(c0sq, c0) ^ vmul(c0sq, c1) ^ c0 ^ c1
    a0 = c0sq ^ c0
    b0 = vmul(c0, c1) ^ c0
    r0 = vmul(c0sq, c0) ^ vmul(c0sq, c1sq) ^ c0sq ^ c0
    cubes = set(mul(mul(x, x), x) for x in range(1, q))
    found = []
    for A in range(q):
        for B in range(1, q):
            ok = irr & ((vmul(a2, A) ^ vmul(b2, B) ^ r2) == 0) & ((vmul(a1, A) ^ vmul(b1, B) ^ r1) == 0) \
                & ((vmul(a0, A) ^ mul(B, B) ^ vmul(b0, B) ^ r0) == 0)
            for i, j in zip(*np.nonzero(ok)):
                found.append((A, B, int(i), int(j)))
    confirmed = all(divides(mul, [f[2], f[3], 1, 1], [mul(f[1], f[1]), 0, 0, 0, 0, f[1], f[0], 0, 0, 0, 1]) for f in found)
    family = all(f[0] == f[1] ^ 1 and f[1] not in cubes for f in found)
    perB = len(set(f[1] for f in found)) == len(found)
    expect = 0 if k % 2 else 2 * (q - 1) // 3
    say("q = %3d: I = %d (expected %d);  pairs ((A,B), irreducible X^3+X^2+c1X+c0 | g_{A,B}): %d (expected %d)%s" %
        (q, I, (q + 1) // 3 if k % 2 else (q - 1) // 3, len(found), expect,
         "" if not found else ";  divisions confirmed: %s;  all with A = B+1, B a non-cube, one per B: %s" % (confirmed, family and perB)))
    return len(found) == expect and (not found or (confirmed and family and perB)) and I == ((q + 1) // 3 if k % 2 else (q - 1) // 3)


if __name__ == "__main__":
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    symbolic()
    resultant_proof()
    say("=== census of irreducible cubic factors with X^2-coefficient 1 ===")
    res = [census(k) for k in range(1, kmax + 1)]
    say("all as expected: %s" % all(res))
