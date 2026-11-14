"""
verify_mean_d67.py -- the mean of Psi over the Oppermann subgroup H_d for d = 6, 7 in characteristic 2.

K = F_q, q = 2^k, E = F_{q^d}, T = Tr_{E/K}.  By the proof of the d = 6, 7 theorem of the paper,
    N_d = q^d + q^{d-2} * sum_{w in E^x, T(w) = 0}  sum_{(eps, a, nu) in ker_w} (-1)^{eps (1 + dk) + Tr_{E/F_2}(nu w^5)},
    ker_w = {(eps, a, nu) in F_2 x K x K : eps + a w^6 + nu w^5 + nu^2 w^10 = 0}          (w^2 = u),
and E_{H_d}[Psi] = N_d / q^{ceil(d/2)}.  The script enumerates ker_w for every w and reports the signed sums of the
branches   I (eps=0, nu=0; only the zero triple),  II (eps=1, nu=0),  III (eps=0, nu!=0),  IV (eps=1, nu!=0),
each split by the degree of w over K, together with the number of kernel elements (unsigned), and compares
N_d with the Hayes FFT where that is available (q^d <= 2^21).

Usage: python verify_mean_d67.py            (q = 2, 4, 8 for d = 6, 7 and q = 16 for d = 6; output also in
                                             ../data/verify_mean_d67.txt)
"""
import os
import sys

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from subgroup_lib import GF2N, class_data, subgroup_indices  # noqa: E402

OUT = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "verify_mean_d67.txt"), "w", encoding="utf8")


def say(s=""):
    print(s, flush=True)
    OUT.write(s + "\n")
    OUT.flush()


@njit(cache=True)
def mul(a, b, poly, N):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> N & 1:
            a ^= poly
    return r


@njit(cache=True)
def power(a, e, poly, N):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a, poly, N)
        a = mul(a, a, poly, N)
        e >>= 1
    return r


@njit(cache=True)
def frob(a, times, poly, N):
    for _ in range(times):
        a = mul(a, a, poly, N)
    return a


@njit(cache=True)
def branches(k, d, poly, N):
    q = 1 << k
    size = 1 << N
    Klist = np.zeros(q, np.int64)
    n = 0
    for x in range(size):
        if frob(x, k, poly, N) == x:
            Klist[n] = x
            n += 1
    signed = np.zeros((4, d + 1), np.int64)      # branch, degree of w over K
    count = np.zeros((4, d + 1), np.int64)
    nadm = 0
    for w in range(1, size):
        t = 0
        y = w
        deg = 0
        for i in range(d):
            t ^= y
            y = frob(y, k, poly, N)
            if deg == 0 and y == w:
                deg = i + 1
        if t != 0:
            continue
        nadm += 1
        w2 = mul(w, w, poly, N)
        w4 = mul(w2, w2, poly, N)
        w5 = mul(w4, w, poly, N)
        w6 = mul(w4, w2, poly, N)
        w10 = mul(w5, w5, poly, N)
        inv6 = power(w6, size - 2, poly, N)
        for eps in range(2):
            for j in range(q):
                nu = Klist[j]
                nw5 = mul(nu, w5, poly, N)
                z = eps ^ nw5 ^ mul(mul(nu, nu, poly, N), w10, poly, N)
                lam = mul(z, inv6, poly, N)
                if frob(lam, k, poly, N) != lam:
                    continue
                tr = 0
                y = nw5
                for i in range(N):
                    tr ^= y
                    y = mul(y, y, poly, N)
                s = 1 - 2 * ((eps * (1 + d * k) + tr) & 1)
                br = (0 if eps == 0 else 1) + (0 if nu == 0 else 2)
                signed[br, deg] += s
                count[br, deg] += 1
    return signed, count, nadm


def main():
    for k, d in ((1, 6), (2, 6), (3, 6), (4, 6), (1, 7), (2, 7), (3, 7)):
        q = 1 << k
        N = k * d
        poly = GF2N(N).poly
        signed, count, nadm = branches(k, d, poly, N)
        assert nadm == q ** (d - 1) - 1
        S = [int(signed[b].sum()) for b in range(4)]
        dual = q ** d + q ** (d - 2) * sum(S)
        mean = dual // q ** ((d + 1) // 2)
        assert mean * q ** ((d + 1) // 2) == dual
        line = "q=%2d d=%d: N_d(dual) = %d, E_H[Psi] - q^d = %d" % (q, d, dual, mean - q ** d)
        if q ** d <= 1 << 21:
            Psi, pi, E, ctx, _ = class_data(q, d)
            Nfft = int(Psi[subgroup_indices(q, d, ctx)].sum())
            line += ";  N_d(FFT) = %d  %s" % (Nfft, "OK" if Nfft == dual else "MISMATCH")
        else:
            line += ";  (no FFT at this size)"
        say(line)
        say("   S_I = %d (= q^(d-1) - 1: %s), S_II = %d, S_III = %d, S_IV = %d" % (S[0], S[0] == q ** (d - 1) - 1, S[1], S[2], S[3]))
        for b, name in ((1, "II"), (2, "III"), (3, "IV")):
            parts = ", ".join("deg %d: %+d (%d terms)" % (e, int(signed[b, e]), int(count[b, e]))
                              for e in range(1, d + 1) if count[b, e])
            say("   branch %-3s by degree of w over K: %s" % (name, parts if parts else "empty"))
        if d == 6:
            a6 = (q - 1) if k % 2 else 3 * (q - 1)
            say("   predicted S_II = -a_6 = %d: %s;  S_III(deg 1) = (q-1)^2 = %d: %s" %
                (-a6, S[1] == -a6, (q - 1) ** 2, int(signed[2, 1]) == (q - 1) ** 2))
            if k % 2:
                say("   k odd: S_III = 2q(q-1) = %d: %s;  S_IV = -2q(q-1): %s;  drift = q^3 - q^2 = %d: %s" %
                    (2 * q * (q - 1), S[2] == 2 * q * (q - 1), S[3] == -2 * q * (q - 1), q ** 3 - q ** 2,
                     mean - q ** d == q ** 3 - q ** 2))
        else:
            say("   predicted: all kernels trivial, drift = q^3 - q = %d: %s" %
                (q ** 3 - q, count[1:].sum() == 0 and mean - q ** d == q ** 3 - q))


if __name__ == "__main__":
    main()
