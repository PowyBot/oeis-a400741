"""Compare a(n) with the naive guess that the q points of the line
alpha^q + alpha = 1 are primitive as often as random nonzero elements of
GF(q^2), i.e. a(n) ~ q * phi(q^2-1)/(q^2-1) / (2n).

usage: python3 heuristic.py results/orbit_count.txt
"""
import sys
from fractions import Fraction


def prime_factors(x):
    out, d = [], 2
    while d * d <= x:
        if x % d == 0:
            out.append(d)
            while x % d == 0:
                x //= d
        d += 1 if d == 2 else 2
    return out + ([x] if x > 1 else [])


for line in open(sys.argv[1]):
    n, value = map(int, line.split()[:2])
    q = 1 << n
    density = Fraction(1)
    for p in prime_factors(q - 1) + prime_factors(q + 1):
        density *= Fraction(p - 1, p)
    expected = q * density / (2 * n)
    print(f"{n:3d} {value:>14d} {float(expected):16.1f} {value / float(expected):7.4f}")
