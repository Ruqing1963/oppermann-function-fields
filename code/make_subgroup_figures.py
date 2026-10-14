"""
make_subgroup_figures.py -- Figures 3 and 4 (Section 7; right panel of Figure 3: Theorem 5.3).

Figure 1 (../figures/subgroup_mean_gauss.pdf)
  left : mean drift on the Oppermann subgroup, (E_H[Psi] - q^d) / q^{d/2}, against d (q = 2, 4, 8, 16);
         the closed forms of Theorem 7.3 at d = 4, 5 are drawn as open circles.
  right: the cubic Gauss-sum sign behind m(c): (3 m(cube) - q)/sqrt q = -eps and (6 m(non-cube) - 2q)/sqrt q = eps,
         eps = (-1)^{k/2}, for q = 4, ..., 2^14 (data/cubic_gauss_m.txt).
Figure 2 (../figures/subgroup_variance.pdf)
  left : Var_G(Psi)/q^d over all classes, q = 2, 3, 4, 5, 7, 8, 16, with the Keating--Rudnick value d - 1;
  right: Var_H(Psi)/q^d over the Oppermann subgroup (q even), same reference line; the exact values of
         Theorem 7.5 at d = 2, 3 are marked.
"""
import json
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
FIG = os.path.join(HERE, "..", "figures")
COL = {2: "#2a78d6", 4: "#eb6834", 8: "#1baf7a", 16: "#eda100", 3: "#e87ba4", 5: "#008300", 7: "#4a3aa7"}
GREY = "#898781"


def style(ax):
    ax.grid(alpha=0.25, linewidth=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def figure1(rows):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 3.9))
    ax.axhline(0, color=GREY, linewidth=0.8, linestyle=":")
    for q in (2, 4, 8, 16):
        pts = [(r["d"], (r["meanH"] - 1) * q ** (r["d"] / 2)) for r in rows if r["q"] == q]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=COL[q], markersize=4, linewidth=1.3,
                label=f"$q={q}$")
        k = q.bit_length() - 1
        th4 = 0.0 if k % 2 == 0 else (-2 * q * q + 2 * q) / q ** 2
        th5 = (q * q - q) / q ** 2.5
        ax.plot([4, 5], [th4, th5], "o", markerfacecolor="none", markeredgecolor="black", markersize=9,
                markeredgewidth=0.9)
    ax.plot([], [], "o", markerfacecolor="none", markeredgecolor="black", label="Theorem 7.3")
    ax.set_xlabel("$d$")
    ax.set_ylabel(r"$(\mathbb{E}_{H_d}[\Psi]-q^d)\,/\,q^{d/2}$")
    ax.set_title("mean drift on the Oppermann subgroup", fontsize=10)
    ax.set_xlim(1.5, 12.5)
    ax.legend(fontsize=7.5, frameon=False, ncol=2)
    style(ax)

    ks, cube, non = [], [], []
    with open(os.path.join(DATA, "cubic_gauss_m.txt"), encoding="utf-8") as fh:
        for ln in fh:
            p = ln.split()
            if not p or not p[0].isdigit():
                continue
            q = int(p[0])
            mc, mn = int(p[4]), int(p[6])
            ks.append(int(math.log2(q)))
            cube.append((3 * mc - q) / math.sqrt(q))
            non.append((6 * mn - 2 * q) / math.sqrt(q))
    bx.axhline(0, color=GREY, linewidth=0.8, linestyle=":")
    bx.plot(ks, cube, "s-", color=COL[4], markersize=6, linewidth=1.3, label=r"$c$ cube: $(3m(c)-q)/\sqrt{q}$")
    bx.plot(ks, non, "o-", color=COL[2], markersize=6, linewidth=1.3,
            label=r"$c$ non-cube: $(6m(c)-2q)/\sqrt{q}$")
    bx.set_xticks(ks)
    bx.set_xticklabels([f"{k}\n{2 ** k}" for k in ks], fontsize=8)
    bx.set_xlabel(r"$k$  /  $q=2^k$")
    bx.set_ylim(-1.6, 1.6)
    bx.set_title(r"cubic Gauss sum sign $\varepsilon=(-1)^{k/2}$ in $m(c)$", fontsize=10)
    bx.legend(fontsize=7.5, frameon=False, loc="lower center", ncol=1)
    style(bx)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "subgroup_mean_gauss.pdf"))


def figure2(rows):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.0), sharey=True)
    dmax = max(r["d"] for r in rows)
    for a in (ax, bx):
        a.plot([2, dmax], [1, dmax - 1], color=GREY, linewidth=1.0, linestyle="--", label=r"Keating–Rudnick $d-1$")
    for q in (2, 3, 4, 5, 7, 8, 16):
        pts = [(r["d"], r["varG"]) for r in rows if r["q"] == q]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=COL[q], markersize=3.5, linewidth=1.2,
                label=f"$q={q}$")
    for q in (2, 4, 8, 16):
        pts = [(r["d"], r["varH"]) for r in rows if r["q"] == q]
        bx.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=COL[q], markersize=3.5, linewidth=1.2,
                label=f"$q={q}$")
        k = q.bit_length() - 1
        bx.plot([2, 3], [q - 1, (1 - 1 / q) * (1 if k % 2 else 3)], "o", markerfacecolor="none",
                markeredgecolor="black", markersize=9, markeredgewidth=0.9)
    bx.plot([], [], "o", markerfacecolor="none", markeredgecolor="black", label="Theorem 7.5")
    ax.set_title(r"all classes: $\mathrm{Var}_{G_d}(\Psi)/q^d$", fontsize=10)
    bx.set_title(r"Oppermann subgroup: $\mathrm{Var}_{H_d}(\Psi)/q^d$", fontsize=10)
    for a in (ax, bx):
        a.set_xlabel("$d$")
        a.set_xlim(1.5, dmax + 0.5)
        a.set_ylim(0, 21)
        style(a)
    ax.set_ylabel("normalised variance")
    ax.legend(fontsize=7, frameon=False, ncol=2, loc="upper left")
    bx.legend(fontsize=7, frameon=False, ncol=2, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "subgroup_variance.pdf"))


def main():
    with open(os.path.join(DATA, "subgroup_variance.json")) as fh:
        rows = json.load(fh)
    figure1(rows)
    figure2(rows)
    print("wrote figures")


if __name__ == "__main__":
    main()
