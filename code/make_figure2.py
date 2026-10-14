"""
make_figure2.py -- Figure 2: characteristic-2 rigidity at d = 2, 3.

Left:  normalised counts pi(x) / (q^d / 2d) against k = log2 q (q = 2..256), exact values from the Hayes FFT
       (hayes_oppermann.interval_psi) with the formulas of Theorem 1.1(B) and Theorem 6.4 as lines:
         d = 2, B_0;  d = 2, B_e (e != 0);  d = 3, C_e (p_3 = 0);
       and the min-max range over the d = 3 Oppermann classes with p_3 != 0.
Right: normalised inverse zeros gamma / sqrt(q) of L(u, chi) for the characters of exact level 2 and 3 of G_3
       (degree 6): q = 8 (characteristic 2; marker area ~ multiplicity) against q = 7 (odd characteristic).
Output: ../figures/oppermann_char2_rigidity.pdf; the distinct characteristic-2 angles are printed.
"""
import os
import sys
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hayes_fft_fq as H  # noqa: E402
from hayes_oppermann import interval_psi, prime_power_E  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BLUE, ORANGE, AQUA, YELLOW, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#898781"
KMAX = 8


def counts(q, d):
    n, l = 2 * d, d
    Psi, ctx = interval_psi(q, l, n)
    E = prime_power_E(q, l, n, ctx)
    cnt = (np.rint(Psi).astype(np.int64) - E) // n
    return cnt, ctx


def left_data():
    ks = list(range(1, KMAX + 1))
    b0, be, c0, lo, hi = [], [], [], [], []
    for k in ks:
        q = 1 << k
        cnt2, ctx2 = counts(q, 2)
        norm2 = q * q / 4
        b0.append(cnt2[H.class_index([1, 0, 0], ctx2)] / norm2)
        be.append(cnt2[H.class_index([1, 0, 1], ctx2)] / norm2)
        cnt3, ctx3 = counts(q, 3)
        norm3 = q ** 3 / 6
        c0.append(cnt3[H.class_index([1, 0, 0, 0], ctx3)] / norm3)
        vals = [cnt3[H.class_index([1, 0, x, y], ctx3)] for x in range(q) for y in range(1, q)]
        lo.append(min(vals) / norm3)
        hi.append(max(vals) / norm3)
        print(f"k = {k}: B0 {b0[-1]:.4f}  Be {be[-1]:.4f}  C0 {c0[-1]:.4f}  p3!=0 range [{lo[-1]:.4f}, {hi[-1]:.4f}]")
        sys.stdout.flush()
    return ks, b0, be, c0, lo, hi


def theory(ks):
    kk = np.array(ks, float)
    odd = (kk % 2 == 1)
    b0 = np.where(odd, 2 * (1 - 2.0 ** -kk), 0.0)
    be = np.where(odd, 1 - 2.0 ** (1 - kk), 1.0)
    c0 = np.where(odd, 1 - 4.0 ** -kk, (1 - 2.0 ** -kk) ** 2)
    return b0, be, c0


def zeros(q, l=3):
    """normalised inverse zeros of L(u, chi) = 1 + c1 u + c2 u^2 for all characters of G_l with deg L >= 1."""
    ctx = H.setup(q, l)
    N = q ** l
    C = []
    for kdeg in range(l):
        ind = H._indices_degree(kdeg, l, q, ctx["p"], ctx["k"], ctx["A"], ctx["M"], ctx["NEG"], ctx["FROB"],
                                ctx["FINV"], ctx["axis_of"], ctx["orders"], ctx["strides"])
        C.append(np.fft.fftn(np.bincount(ind, minlength=N).astype(complex).reshape(ctx["shape"])).reshape(N))
    out = []
    for i in range(1, N):
        c1, c2 = C[1][i], C[2][i]
        if abs(c2) > 1e-9:
            out.extend(np.roots([1, c1, c2]))
        elif abs(c1) > 1e-9:
            out.append(-c1)
    return np.array(out) / np.sqrt(q)


def main():
    ks, b0, be, c0, lo, hi = left_data()
    tb0, tbe, tc0 = theory(ks)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [1.35, 1]})

    ax.axhline(1.0, color=GREY, linewidth=0.8, linestyle=":")
    ax.vlines(ks, lo, hi, color=YELLOW, linewidth=6, alpha=0.45, label=r"$d=3$, $p_3\neq0$ (min–max)")
    ax.plot(ks, tb0, color=BLUE, linewidth=1.4, linestyle="--")
    ax.plot(ks, b0, "o", color=BLUE, markersize=6, label=r"$d=2$, $B_0$")
    ax.plot(ks, tbe, color=ORANGE, linewidth=1.4, linestyle="--")
    ax.plot(ks, be, "s", color=ORANGE, markersize=5, label=r"$d=2$, $B_e$ ($e\neq0$)")
    ax.plot(ks, tc0, color=AQUA, linewidth=1.4)
    ax.plot(ks, c0, "D", color=AQUA, markersize=5, label=r"$d=3$, $C_e$ ($p_3=0$)")
    ax.set_xticks(ks)
    ax.set_xticklabels([f"{k}\n{1 << k}" for k in ks], fontsize=8)
    ax.set_xlabel(r"$k=\log_2 q$  /  $q$")
    ax.set_ylabel(r"$\pi(x)\,/\,(q^d/2d)$")
    ax.set_title("characteristic 2: parity oscillation", fontsize=10)
    ax.set_ylim(-0.08, 2.75)
    ax.grid(alpha=0.25, linewidth=0.5)
    ax.legend(fontsize=7.5, frameon=False, loc="upper center", ncol=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    th = np.linspace(0, 2 * np.pi, 400)
    bx.plot(np.cos(th), np.sin(th), color=GREY, linewidth=0.8)
    z7 = zeros(7)
    bx.scatter(z7.real, z7.imag, s=4, color=GREY, alpha=0.35, linewidths=0, label=r"$q=7$ (odd)")
    z8 = zeros(8)
    key = Counter((round(z.real, 6), round(z.imag, 6)) for z in z8)
    pts = np.array(list(key.keys()))
    mult = np.array(list(key.values()))
    bx.scatter(pts[:, 0], pts[:, 1], s=12 + 1.2 * mult, color=BLUE, alpha=0.85, linewidths=0,
               label=r"$q=8$ (char. 2)")
    ang = sorted({round((np.degrees(np.angle(complex(*p))) + 360) % 360, 3) for p in pts})
    print("q = 8: distinct normalised zeros", len(pts), "at angles (deg):", ang, " |z| =",
          sorted({round(abs(complex(*p)), 6) for p in pts}))
    print("q = 7: zeros", len(z7), " |z| range", round(min(abs(z7)), 6), round(max(abs(z7)), 6))
    bx.set_aspect("equal")
    bx.set_xlim(-1.35, 1.35)
    bx.set_ylim(-1.35, 1.35)
    bx.set_title(r"inverse zeros $\gamma/\sqrt{q}$, levels 2 and 3, $d=3$", fontsize=10)
    bx.legend(fontsize=7.5, frameon=False, loc="upper left", handletextpad=0.2, borderaxespad=0.0)
    bx.grid(alpha=0.2, linewidth=0.5)
    for s in ("top", "right"):
        bx.spines[s].set_visible(False)

    fig.tight_layout()
    out = os.path.join(HERE, "..", "figures", "oppermann_char2_rigidity.pdf")
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
