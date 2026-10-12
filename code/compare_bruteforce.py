"""
compare_bruteforce.py -- class-by-class comparison of the Hayes-FFT counts (hayes_oppermann.py) with the brute-force
Ben-Or counts (verify_pointwise_oppermann.py) for every short interval {prefix + s : deg s < d} of degree 2d.
Usage: python compare_bruteforce.py [max_log2_size=22]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fqlib  # noqa: E402
import hayes_fft_fq as H  # noqa: E402
from hayes_oppermann import run  # noqa: E402
from verify_pointwise_oppermann import FIELDS, irreducible_flags  # noqa: E402


def main():
    max_log2 = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    for q in FIELDS:
        d = 1
        while 2 * d * np.log2(q) <= max_log2:
            _, cnt, ctx = run(q, d)
            A, M, NEG, INV = fqlib.field_tables(q)
            flags = irreducible_flags(2 * d, q, A, M, NEG, INV)
            brute = flags.reshape(q ** d, q ** d).sum(axis=1, dtype=np.int64)      # brute[prefix]
            mism = 0
            for pidx in range(q ** d):
                series = [1] + [0] * d
                r = pidx
                for i in range(d):                     # digit i = coefficient of t^{d+i} = series index d - i
                    series[d - i] = r % q
                    r //= q
                if cnt[H.class_index(series, ctx)] != brute[pidx]:
                    mism += 1
            print(f"q = {q:2d}  d = {d:2d}  classes {q ** d:8d}  mismatches {mism}")
            sys.stdout.flush()
            d += 1


if __name__ == "__main__":
    main()
