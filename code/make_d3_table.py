"""
make_d3_table.py -- Table 2: the d = 3 class counts in characteristic 2:
all characteristic-2 Oppermann classes C_{e,c} of degree 6, q = 2^k, k <= 14.

Values from Theorem 5.2 and Corollary 5.4:
  k odd : q(q^2-1)/6 for every class;
  k even: q(q-1)^2/6 (c = 0),  q(q^2+1-2 eps sqrt q)/6 (c a cube),  q(q^2+1+eps sqrt q)/6 (c a non-cube),
          eps = (-1)^{k/2}.
For q <= 256 these agree with the exact Hayes-FFT counts of every class (char2_d3_all.py, data/char2_d3_all.txt);
for q = 1024, 4096, 16384 the two values of m(c) were checked directly (cubic_gauss_m.py).
Output: ../paper/tab_d3.tex
"""
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def num(x):
    return f"{x:,}".replace(",", r"\,")


def main():
    rows = [r"\begin{table}[tbp]", r"\centering",
            r"\caption{Every characteristic-$2$ Oppermann class of degree $6$: the number of irreducibles in $C_{e,c}$ "
            r"(Theorem~\ref{thm:p3zero}, Corollary~\ref{cor:d3}). For $k$ odd all classes have the same count. For "
            r"$q\le256$ every class was also counted exactly with the Hayes FFT; for $q\ge1024$ the two values of $m(c)$ "
            r"were checked directly.}\label{tab:d3}",
            r"\begin{tabular}{rrrrrr}", r"\toprule",
            r"$q$ & $\varepsilon$ & $c=0$ & $c$ cube & $c$ non-cube & $q^3/6$\\", r"\midrule"]
    for k in range(1, 15):
        if k > 8 and k % 2:
            continue
        q = 1 << k
        if k % 2:
            v = q * (q * q - 1) // 6
            rows.append(f"{num(q)} & -- & \\multicolumn{{3}}{{c}}{{{num(v)}}} & {num(round(q ** 3 / 6))}\\\\")
        else:
            s = math.isqrt(q)
            eps = (-1) ** (k // 2)
            c0 = q * (q - 1) ** 2 // 6
            cube = q * (q * q + 1 - 2 * eps * s) // 6
            non = q * (q * q + 1 + eps * s) // 6
            rows.append(f"{num(q)} & ${'+' if eps > 0 else '-'}1$ & {num(c0)} & {num(cube)} & {num(non)} & "
                        f"{num(round(q ** 3 / 6))}\\\\")
    rows += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    with open(os.path.join(HERE, "..", "paper", "tab_d3.tex"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(rows) + "\n")
    print("wrote tab_d3.tex")


if __name__ == "__main__":
    main()
