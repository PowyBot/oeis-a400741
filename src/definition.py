"""A400741 by its literal definition, for small n.

Counts monic g of degree n over GF(2) such that g(x^2+x) is primitive.
Polynomials are Python ints (bit i = coefficient of x^i). Slow but shares
nothing with the C programs, so it checks them on n <= 12 or so.

usage: python3 definition.py N   -> prints a(1..N)
"""
import sys


def mulmod(a, b, f, deg):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if (a >> deg) & 1:
            a ^= f
    return r


def powmod(a, e, f, deg):
    r = 1
    while e:
        if e & 1:
            r = mulmod(r, a, f, deg)
        a = mulmod(a, a, f, deg)
        e >>= 1
    return r


def prime_factors(x):
    out, d = [], 2
    while d * d <= x:
        if x % d == 0:
            out.append(d)
            while x % d == 0:
                x //= d
        d += 1
    return out + ([x] if x > 1 else [])


def compose_with_x2_plus_x(g):
    """g(x^2 + x) as a polynomial, by Horner's rule."""
    def mul(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            a <<= 1
            b >>= 1
        return r

    r = 0
    for i in range(g.bit_length() - 1, -1, -1):
        r = mul(r, 0b110) ^ ((g >> i) & 1)
    return r


def is_primitive(f, deg, cofactors):
    """x has order 2^deg - 1 modulo f. That already forces f irreducible:
    otherwise the unit group of GF(2)[x]/f has no element of that order."""
    order = (1 << deg) - 1
    if powmod(2, order, f, deg) != 1:
        return False
    return all(powmod(2, c, f, deg) != 1 for c in cofactors)


def a(n):
    deg = 2 * n
    order = (1 << deg) - 1
    cofactors = [order // p for p in prime_factors(order)]
    count = 0
    for low in range(1 << n):
        f = compose_with_x2_plus_x((1 << n) | low)
        count += is_primitive(f, deg, cofactors)
    return count


if __name__ == "__main__":
    print(" ".join(str(a(n)) for n in range(1, int(sys.argv[1]) + 1)))
