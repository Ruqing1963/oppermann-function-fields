"""
make_figure.py -- normalised minimum of the irreducible counts over Oppermann sub-intervals, min / (q^d / 2d),
against d, from ../data/hayes_oppermann_census.txt.  Two panels: odd q, even q (empty cases plotted at 0).
Output: ../figures/oppermann_min_ratio.pdf
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]


def load():
    rows = {}
    with open(os.path.join(HERE, "..", "data", "hayes_oppermann_census.txt"), encoding="utf-8-sig") as fh:
        for ln in fh:
            parts = ln.split()
            if len(parts) < 10 or not parts[0].isdigit():
                continue
            q, d, ratio = int(parts[0]), int(parts[1]), float(parts[5])
            rows.setdefault(q, {})[d] = ratio
    return rows


def main():
    rows = load()
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
    for ax, qs, title in ((axes[0], [3, 5, 7, 9, 11, 13], "odd $q$ (all classes)"),
                          (axes[1], [2, 4, 8, 16], "even $q$ (Oppermann classes)")):
        for i, q in enumerate(qs):
            if q not in rows:
                continue
            ds = sorted(rows[q])
            ax.plot(ds, [rows[q][d] for d in ds], marker="o", markersize=4, linewidth=1.6, color=SERIES[i],
                    label=f"$q={q}$")
        ax.axhline(1.0, color="#898781", linewidth=0.8, linestyle=":")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("$d$")
        ax.grid(alpha=0.25, linewidth=0.5)
        ax.legend(fontsize=8, frameon=False, loc="lower right")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel(r"$\min_{f,c}\,\pi(O_{f,c})\,/\,(q^d/2d)$")
    axes[0].set_ylim(-0.03, 1.05)
    fig.tight_layout()
    out = os.path.join(HERE, "..", "figures", "oppermann_min_ratio.pdf")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
