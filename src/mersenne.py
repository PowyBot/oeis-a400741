"""a(n) when 2^n - 1 is a Mersenne prime.

Then every element of GF(2^n) other than 0 and 1 is primitive. In the orbit
formulation (see orbit_count.c) V_k(t) = u^k + u^-k is never 0 (u^k != 1) and
never 1 (that would make u^k a cube root of unity, but it has order 2^n + 1),
so every orbit counts:
    a(n) = phi(2^n + 1) / (2n).
This prints the formula's value for the Mersenne exponents below 128, with
2^n + 1 factored by Pollard-Brent. Only the factorization is probabilistic
(Miller-Rabin on the factors); the values for n <= 38 are checked against
the full count in results/.

usage: python3 mersenne.py
"""
import math
import random

MERSENNE_EXPONENTS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127]


def is_probable_prime(n):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def brent(n):
    """R. P. Brent, "An improved Monte Carlo factorization algorithm",
    BIT 20 (1980) 176-184, doi:10.1007/BF01933190."""
    if n % 2 == 0:
        return 2
    while True:
        y, c, m = random.randrange(1, n), random.randrange(1, n), 128
        g = r = q = 1
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
        if g != n:
            return g


def factor(n, out):
    if n == 1:
        return
    if is_probable_prime(n):
        out[n] = out.get(n, 0) + 1
        return
    for p in range(2, 10000):
        if n % p == 0:
            factor(p, out)
            factor(n // p, out)
            return
    d = brent(n)
    factor(d, out)
    factor(n // d, out)


for n in MERSENNE_EXPONENTS:
    primes = {}
    factor((1 << n) + 1, primes)
    phi = 1
    for p, e in primes.items():
        phi *= (p - 1) * p ** (e - 1)
    assert phi % (2 * n) == 0
    print(n, phi // (2 * n), " 2^n+1 = " + " * ".join(
        f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(primes.items())))
