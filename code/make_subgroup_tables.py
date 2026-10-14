"""
make_subgroup_tables.py -- LaTeX tables of Section 7 (the Oppermann subgroup), generated from ../data:
  ../paper/tab_subgroup_means.tex     subgroup means (verify_subgroup_mean.txt)
  ../paper/tab_subgroup_variance.tex  subgroup and global variances, characteristic 2 (subgroup_variance.json)
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PAPER = os.path.join(HERE, "..", "paper")


def num(x):
    s = f"{x:,}"
    return s.replace(",", r"\,")


def frac(s):
    if "/" in s:
        a, b = s.split("/")
        return rf"{num(int(a))}/{b}"
    return num(int(s))


def means():
    rows = []
    with open(os.path.join(DATA, "verify_subgroup_mean.txt"), encoding="utf-8") as fh:
        for ln in fh:
            p = ln.split()
            if not p or not p[0].isdigit():
                continue
            rows.append(p)
    out = [r"\begin{table}[tbp]", r"\centering",
           r"\caption{Means over the Oppermann subgroup $H_d$ (exact, from the Hayes FFT) against the closed forms of "
           r"Corollary~\ref{cor:d23} and Theorem~\ref{thm:d45}. Every entry agrees.}\label{tab:means}",
           r"\begin{tabular}{rrrrrr}", r"\toprule",
           r"$q$ & $d$ & $|H_d|$ & $\mathbb{E}_{H_d}[\Psi]$ & $\mathbb{E}_{H_d}[\pi]$ & $q^d$\\", r"\midrule"]
    last = None
    for p in rows:
        q, d = int(p[0]), int(p[1])
        if last is not None and q != last:
            out.append(r"\addlinespace")
        last = q
        out.append(f"{q} & {d} & {num(int(p[2]))} & {frac(p[3])} & {frac(p[5])} & {num(q ** d)}\\\\")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(out)


def variance():
    with open(os.path.join(DATA, "subgroup_variance.json")) as fh:
        rows = json.load(fh)
    out = [r"\begin{table}[tbp]", r"\centering",
           r"\caption{Normalised variances in characteristic $2$: over the Oppermann subgroup $H_d$ and over all "
           r"classes $G_d$, for $\Psi$ and for $2d\pi$, with the Keating--Rudnick value $d-1$. The entries at $d=2,3$ "
           r"agree exactly with Theorem~\ref{thm:var23}. For $q=2$ only $d\le12$ is shown.}\label{tab:var}",
           r"\begin{tabular}{rrrrrrrr}", r"\toprule",
           r"$q$ & $d$ & $|H_d|$ & $\frac{\mathrm{Var}_{H}(\Psi)}{q^d}$ & $\frac{\mathrm{Var}_{H}(2d\pi)}{q^d}$ & "
           r"$\frac{\mathrm{Var}_{G}(\Psi)}{q^d}$ & $\frac{\mathrm{Var}_{G}(2d\pi)}{q^d}$ & $d-1$\\", r"\midrule"]
    last = None
    for r in rows:
        if r["q"] % 2 or (r["q"] == 2 and r["d"] > 12):
            continue
        if last is not None and r["q"] != last:
            out.append(r"\addlinespace")
        last = r["q"]
        out.append(f"{r['q']} & {r['d']} & {num(r['H'])} & {r['varH']:.3f} & {r['varHpi']:.3f} & {r['varG']:.3f} & "
                   f"{r['varGpi']:.3f} & {r['d'] - 1}\\\\")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(out)


def main():
    for name, f in (("tab_subgroup_means.tex", means), ("tab_subgroup_variance.tex", variance)):
        with open(os.path.join(PAPER, name), "w", encoding="utf-8") as fh:
            fh.write(f() + "\n")
    print("wrote tables")


if __name__ == "__main__":
    main()
