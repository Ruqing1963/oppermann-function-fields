"""
make_tables.py -- LaTeX tables for the paper from ../data/hayes_oppermann_census.txt -> ../paper/census_tables.tex
(the last few degrees for each q; all degrees are in the data file).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    rows = {}
    with open(os.path.join(HERE, "..", "data", "hayes_oppermann_census.txt"), encoding="utf-8-sig") as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 10 or not p[0].isdigit():
                continue
            q, d = int(p[0]), int(p[1])
            rows.setdefault(q, {})[d] = dict(ncl=int(p[2]), mn=int(p[3]), mean=float(p[4]), ratio=float(p[5]),
                                             empty=int(p[6]))
    return rows


def table(rows, qs, last, label, caption):
    out = [r"\begin{table}[!htbp]", r"\centering", rf"\caption{{{caption}}}\label{{{label}}}",
           r"\begin{tabular}{rrrrrr}", r"\toprule",
           r"$q$ & $d$ & classes & $\min\pi$ & mean $\pi$ & $\min\pi/(q^d/2d)$\\", r"\midrule"]
    for q in qs:
        ds = sorted(rows[q])
        dmax = ds[-1]
        show = [d for d in ds if d > dmax - last]
        for i, d in enumerate(show):
            r = rows[q][d]
            qcol = f"{q}" if i == 0 else ""
            mn = r"\textbf{0}" if r["mn"] == 0 else f"{r['mn']:,}"
            ratio = r"\textbf{0}" if r["mn"] == 0 else f"{r['ratio']:.4f}"
            out.append(f"{qcol} & {d} & {r['ncl']:,} & {mn} & {r['mean']:,.1f} & {ratio}\\\\"
                       .replace(",", "\\,"))
        if q != qs[-1]:
            out.append(r"\addlinespace")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(out)


def main():
    rows = load()
    odd = table(rows, [3, 5, 7, 9, 11, 13], 4, "tab:odd",
                r"Odd $q$: irreducible counts $\pi$ over all $q^d$ Hayes classes of level $d$ in degree $2d$ (all are "
                r"Oppermann sub-intervals). Last four degrees for each $q$; all smaller degrees are in the data file, "
                r"with minimum $\ge1$ throughout.")
    even = table(rows, [2, 4, 8, 16], 4, "tab:even",
                 r"Even $q$: irreducible counts over the $q^{\lceil d/2\rceil}$ Oppermann classes (Lemma~\ref{lem:even}). "
                 r"Last four degrees for each $q$; the only empty classes in the whole range are those of "
                 r"Theorem~\ref{thm:main2}.")
    with open(os.path.join(HERE, "..", "paper", "census_tables.tex"), "w", encoding="utf-8") as fh:
        fh.write(odd + "\n\n" + even + "\n")
    print("wrote census_tables.tex")


if __name__ == "__main__":
    main()
