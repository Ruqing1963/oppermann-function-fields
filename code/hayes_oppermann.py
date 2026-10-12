"""
hayes_oppermann.py -- exact irreducible counts in all Oppermann sub-intervals via the Hayes explicit formula
(the FFT algorithm of Part I, hayes_fft_fq.py, run at level l = d instead of l = d - 1).

The block {P + s : deg s < d} (P monic of degree n = 2d) is the Hayes class of level l = d determined by the
coefficients of t^{2d-1}, ..., t^d.  For every class x:
    Psi(x) = sum_{F in x, deg F = n} Lambda(F) = q^{-l} sum_chi conj(chi(x)) psi_n(chi),
with psi_n(chi) obtained from the L-polynomials L(u, chi) (degree <= l - 1) by the Newton recursion; the
L-coefficients are the FFTs of the class distributions of the monic polynomials of degree < l.  Then
    #irreducibles in x = (Psi(x) - E(x)) / n,     E(x) = sum_{k >= 2, k | n} sum_{deg pi = n/k, pi^k in x} n/k,
where the prime powers are enumerated directly (Ben-Or flags in degree <= d) and placed by the homomorphism
class(pi^k) = k * class(pi) in the Hayes group.  Psi is rounded; the largest rounding residual is reported.

Oppermann sub-intervals: O_{f,c} = class of f^2 + c t^d.
  * q odd: f -> top d coefficients of f^2 is a bijection (triangular with diagonal 2a_i), so every class is an
    Oppermann sub-interval: the minimum is taken over all q^d classes.
  * q even: f^2 has only even-degree terms; the Oppermann classes are those whose coefficients at odd degrees
    j in (d, 2d-1] vanish (q^{1 + #even j in (d, 2d-1]} classes); enumerated explicitly.

Usage: python hayes_oppermann.py q d_min d_max [work_dir]
Output: one line per (q, d): #Oppermann classes, min, mean, min / (q^d / n), #empty, max rounding residual, and
(when q is even) also the minimum over all classes.
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fqlib  # noqa: E402
import hayes_fft_fq as H  # noqa: E402
from hayes_fft_fq import _dlog  # noqa: E402
from verify_pointwise_oppermann import irreducible_flags  # noqa: E402


def interval_psi(q, l, n, work_dir=None):
    """Psi(x) for every Hayes class x of level l, degree n (flat index, float64)."""
    import scipy.fft as sfft
    ctx = H.setup(q, l)
    N = q ** l
    shape = ctx["shape"]
    big = work_dir is not None and N * l * 16 > 600 * 2 ** 20
    C = []
    for kdeg in range(0, l):
        ind = H._indices_degree(kdeg, l, q, ctx["p"], ctx["k"], ctx["A"], ctx["M"], ctx["NEG"],
                                ctx["FROB"], ctx["FINV"], ctx["axis_of"], ctx["orders"], ctx["strides"])
        arr = np.bincount(ind, minlength=N).astype(np.complex128)
        del ind
        F = sfft.fftn(arr.reshape(shape), overwrite_x=True, workers=-1).reshape(N)
        del arr
        if big:
            mm = np.memmap(os.path.join(work_dir, f"c{q}_{l}_{kdeg}.dat"), dtype=np.complex128, mode="w+",
                           shape=(N,))
            mm[:] = F
            del F
            C.append(mm)
        else:
            C.append(F)
    top = np.empty(N, np.complex128)
    chunk = 1 << 20
    for s0 in range(0, N, chunk):
        s1 = min(N, s0 + chunk)
        c = [np.asarray(C[t][s0:s1]) for t in range(l)]
        psi = [None] * (n + 1)
        for m in range(1, n + 1):
            acc = m * c[m] if m < l else np.zeros(s1 - s0, np.complex128)
            for kk in range(max(1, m - l + 1), m):
                acc = acc - psi[kk] * c[m - kk]
            psi[m] = acc
        top[s0:s1] = psi[n]
    top[0] = float(q) ** n
    del C
    if big:
        for kdeg in range(l):
            try:
                os.remove(os.path.join(work_dir, f"c{q}_{l}_{kdeg}.dat"))
            except OSError:
                pass
    Psi = sfft.ifftn(top.reshape(shape), overwrite_x=True, workers=-1).reshape(N).real
    return Psi, ctx


@njit(cache=True)
def _power_corrections(flags, m, kpow, l, q, p, k, A, M, NEG, FROB, FINV, axis_of, orders, strides, out):
    """out[class(pi^kpow)] += m for every monic irreducible pi of degree m flagged in `flags` (index = low digits)."""
    series = np.zeros(l + 1, np.int64)
    P = np.zeros(m + 1, np.int64)
    nax = len(strides)
    for idx in range(q ** m):
        if flags[idx] == 0:
            continue
        r = idx
        for i in range(m):
            P[i] = r % q
            r //= q
        P[m] = 1
        for i in range(l + 1):
            series[i] = 0
        series[0] = 1
        for i in range(1, min(m, l) + 1):
            series[i] = P[m - i]          # coefficient of t^{m-i}
        flat = _dlog(series, l, p, k, A, M, NEG, FROB, FINV, axis_of, orders, strides)
        nf = 0
        for a in range(nax):
            e = (flat // strides[a]) % orders[a]
            nf += ((e * kpow) % orders[a]) * strides[a]
        out[nf] += m


def prime_power_E(q, l, n, ctx):
    A, M, NEG, INV = fqlib.field_tables(q)
    E = np.zeros(q ** l, np.int64)
    for kpow in range(2, n + 1):
        if n % kpow:
            continue
        m = n // kpow
        flags = irreducible_flags(m, q, A, M, NEG, INV)
        _power_corrections(flags, m, kpow, l, q, ctx["p"], ctx["k"], ctx["A"], ctx["M"], ctx["NEG"],
                           ctx["FROB"], ctx["FINV"], ctx["axis_of"], ctx["orders"], ctx["strides"], E)
        if m == 1:
            E[0] += 1                     # pi = t (skipped by the flags, which require P(0) != 0): t^n is in class 1
    return E


def oppermann_classes_char2(q, d, ctx):
    """flat indices of the Oppermann classes f^2 + c t^d in characteristic 2."""
    free = [d] + [j for j in range(d + 1, 2 * d) if j % 2 == 0]
    out = set()
    nfree = len(free)
    for code in range(q ** nfree):
        series = [0] * (d + 1)
        series[0] = 1
        r = code
        for j in free:
            series[2 * d - j] = r % q      # series index i <-> coefficient of t^{2d-i}
            r //= q
        out.add(H.class_index(series, ctx))
    return np.array(sorted(out), np.int64)


def run(q, d, work_dir=None):
    t0 = time.time()
    n, l = 2 * d, d
    Psi, ctx = interval_psi(q, l, n, work_dir)
    N = Psi.shape[0]
    cnt = np.empty(N, np.int64)
    resid = 0.0
    for s0 in range(0, N, 1 << 22):                    # chunked to keep the peak memory low
        x = np.ascontiguousarray(Psi[s0:s0 + (1 << 22)])
        r = np.rint(x)
        resid = max(resid, float(np.max(np.abs(x - r))))
        cnt[s0:s0 + len(r)] = r.astype(np.int64)
    del Psi
    E = prime_power_E(q, l, n, ctx)
    cnt -= E
    assert np.all(cnt % n == 0), "non-integral prime count"
    cnt //= n
    assert int(cnt.sum()) * n + int(E.sum()) == q ** n
    if q % 2:
        sel = cnt
    else:
        sel = cnt[oppermann_classes_char2(q, d, ctx)]
    res = dict(q=q, d=d, n_classes=len(sel), min=int(sel.min()), mean=float(sel.mean()),
               ratio=float(sel.min()) / (q ** d / n), n_empty=int((sel == 0).sum()), resid=resid,
               all_min=int(cnt.min()), secs=time.time() - t0)
    return res, cnt, ctx


def main():
    q = int(sys.argv[1])
    d0, d1 = int(sys.argv[2]), int(sys.argv[3])
    work_dir = sys.argv[4] if len(sys.argv) > 4 else None
    print(f"{'q':>3} {'d':>3} {'#Opp classes':>13} {'min':>9} {'mean':>12} {'min/(q^d/2d)':>13} {'#empty':>7} "
          f"{'all-class min':>13} {'resid':>9} {'secs':>8}")
    for d in range(d0, d1 + 1):
        r, _, _ = run(q, d, work_dir)
        print(f"{r['q']:>3} {r['d']:>3} {r['n_classes']:>13} {r['min']:>9} {r['mean']:>12.2f} {r['ratio']:>13.4f} "
              f"{r['n_empty']:>7} {r['all_min']:>13} {r['resid']:>9.2e} {r['secs']:>8.1f}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
