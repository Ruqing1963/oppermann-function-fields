"""
subgroup_lib.py -- shared routines for the Oppermann-subgroup paper.

  class_data(q, d)        exact Psi(x) and pi(x) for every Hayes class x of level d in degree 2d (Hayes FFT)
  subgroup_indices(q, d)  flat indices of H_d (classes with e_j = 0 for odd j < d), q even
  prefix_of_index(q, d)   the prefix (e_1, ..., e_d) of every flat class index
  power_sums(e, q, tabs)  S_1..S_d from e_1..e_d by Newton's identities in characteristic 2
  GF2N                    arithmetic in F_{2^N} with integers (carry-less multiplication), for brute-force sums
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fqlib  # noqa: E402
import hayes_fft_fq as H  # noqa: E402
from hayes_oppermann import interval_psi, oppermann_classes_char2, prime_power_E  # noqa: E402


def class_data(q, d, work_dir=None):
    n, l = 2 * d, d
    Psi, ctx = interval_psi(q, l, n, work_dir)
    N = Psi.shape[0]
    R = np.empty(N, np.int64)
    resid = 0.0
    for s0 in range(0, N, 1 << 22):
        x = np.ascontiguousarray(Psi[s0:s0 + (1 << 22)])
        r = np.rint(x)
        resid = max(resid, float(np.max(np.abs(x - r))))
        R[s0:s0 + len(r)] = r.astype(np.int64)
    del Psi
    E = prime_power_E(q, l, n, ctx)
    num = R - E
    assert np.all(num % n == 0)
    assert int(R.sum()) == q ** n
    return R, num // n, E, ctx, resid


def subgroup_indices(q, d, ctx):
    return oppermann_classes_char2(q, d, ctx)


def prefix_of_index(q, d, ctx):
    """array P[idx, i] = e_{i+1} for every class index (prefix = series coefficients 1..d)."""
    N = q ** d
    P = np.zeros((N, d), np.int64)
    for code in range(N):
        series = [1] + [(code // q ** i) % q for i in range(d)]
        P[H.class_index(series, ctx), :] = series[1:]
    return P


def power_sums(e, q, A, M):
    """S_1..S_d (list) from e_1..e_d over F_q, q even: S_j = sum_{i<j} e_i S_{j-i} + j e_j."""
    d = len(e)
    S = [0] * (d + 1)
    for j in range(1, d + 1):
        s = e[j - 1] if j % 2 else 0
        for i in range(1, j):
            s = int(A[s, M[e[i - 1], S[j - i]]])
        S[j] = s
    return S[1:]


def trace_table(q, A, M):
    k = q.bit_length() - 1
    tr = np.zeros(q, np.int64)
    for x in range(q):
        s, y = 0, x
        for _ in range(k):
            s = int(A[s, y])
            y = int(M[y, y])
        tr[x] = s
    return tr


class GF2N:
    """F_{2^N} = F_2[x]/(f), f found by search (Rabin's test)."""

    def __init__(self, N, poly=None):
        self.N = N
        self.poly = poly if poly is not None else self._find_irreducible(N)

    @staticmethod
    def _clmul(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            a <<= 1
            b >>= 1
        return r

    @staticmethod
    def _polymod(a, f):
        df = f.bit_length() - 1
        while a and a.bit_length() - 1 >= df:
            a ^= f << (a.bit_length() - 1 - df)
        return a

    @classmethod
    def _gcd(cls, a, b):
        while b:
            a, b = b, cls._polymod(a, b)
        return a

    @classmethod
    def _find_irreducible(cls, N):
        primes = [p for p in range(2, N + 1) if N % p == 0 and all(p % r for r in range(2, p))]
        for f in range((1 << N) | 1, 1 << (N + 1), 2):
            def xpow(e):  # x^(2^e) mod f
                y = 2
                for _ in range(e):
                    y = cls._polymod(cls._clmul(y, y), f)
                return y
            if xpow(N) != 2:
                continue
            if all(cls._gcd(f, xpow(N // p) ^ 2) == 1 for p in primes):
                return f
        raise ValueError

    def mul(self, a, b):
        return self._polymod(self._clmul(a, b), self.poly)

    def pow(self, a, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            e >>= 1
        return r

    def trace(self, a):
        s, y = 0, a
        for _ in range(self.N):
            s ^= y
            y = self.mul(y, y)
        return s  # 0 or 1
