# A400741: more terms

[OEIS A400741](https://oeis.org/A400741) counts the monic polynomials g of
degree n over GF(2) such that g(x^2+x) is primitive; equivalently, the
primitive polynomials f of degree 2n with f(x+1) = f(x). When I looked
(2026-10-16) the entry had a(1)..a(22), keyword `more`.

This repo computes **a(1)..a(@LAST@)** with two independent C programs that
agree with each other and with all 22 published terms.

## Results

@TABLE@

The b-file is [`bfile/b400741.txt`](bfile/b400741.txt); a draft of the text I
would submit to OEIS is in [`bfile/submission.md`](bfile/submission.md).
Nothing has been submitted.

## Method

Let q = 2^n. If alpha is a root of a primitive f with f(x+1) = f(x), then
alpha+1 is a conjugate alpha^(2^j); applying the Frobenius twice forces j = n.
So the roots are exactly the primitive elements on the affine line
alpha^q + alpha = 1 in GF(q^2), 2n of them per polynomial.

**`line_count`** (cross-check, n <= 30): walks the q points of that line in
Gray-code order and tests each for primitivity in GF(2^(2n)). Cost about
q * (a few exponentiations), so it stops being practical near n = 31.

**`orbit_count`** (main count): uses the formula already in the entry,
a(n) = #{u : ord(u) = q+1, u + u^q primitive in GF(q)} / (2n). Fix one u of
order q+1 and t = u + 1/u. Then u^k + u^-k = V_k(t), a Dickson polynomial,
computable by a ladder in GF(q) alone. Squaring sends V_k to V_2k, and
2^n = -1 (mod q+1), so k -> 2k acts on the units mod q+1 with orbits of
exactly 2n elements, and primitivity is constant on each orbit. Hence

    a(n) = number of orbits of k -> 2k on (Z/(q+1))* whose V_k(t) is primitive.

The program visits each orbit once, at its least element (found by a few
integer doublings), and does one ladder plus one primitivity test per orbit:
roughly 2n times less field arithmetic than the PARI program in the entry.
As a completeness check it confirms that it visited exactly phi(q+1)/(2n)
orbits.

Field arithmetic is GF(2^m) in a 64-bit word with `PCLMULQDQ`, modulo the
irreducible x^m + r with the smallest r (found with Ben-Or's test).

## A closed form at Mersenne exponents

If 2^n - 1 is prime, every element of GF(q) other than 0 and 1 is primitive,
and V_k(t) is never 0 or 1 for a unit k (that would make u^k equal to 1 or to
a cube root of unity, but it has order q+1). So every orbit counts:

    a(n) = phi(2^n + 1) / (2n)   when 2^n - 1 is prime.

This matches the full count at n = 2, 3, 5, 7, 13, 17, 19, 31, and gives for
example a(61) = 12600235023025650 and
a(127) = 446564785985483547852197647548252246
(`src/mersenne.py`, using the known factorizations of 2^n + 1, re-derived
there by Pollard-Brent).

## Heuristic

If the points of the line were primitive as often as random elements of
GF(q^2)*, a(n) would be about q * phi(q^2-1)/(q^2-1) / (2n). The ratio of
a(n) to that guess stays within 0.1% for n >= 23 (`src/heuristic.py`). This
is an observation, not a theorem.

## Verification

1. `orbit_count` and `line_count` reproduce all 22 published terms.
2. They agree with each other for every n <= 30, including the new
   a(23)..a(30). They share only the field multiplication; the line method
   does not use the formula from the entry at all.
3. `src/definition.py` (pure Python, the literal definition: build
   g(x^2+x), test x for order 2^(2n)-1) agrees for n <= 11.
4. For n > 30 only `orbit_count` was run; its orbit total is checked against
   phi(q+1)/(2n) every time, and at n = 31 the result equals the Mersenne
   closed form.

## Reproduce

    make
    ./run.sh orbit_count 1 38 4     # appends "n a(n) seconds" to results/
    ./run.sh line_count 1 30 4
    python3 src/definition.py 11
    python3 src/heuristic.py results/orbit_count.txt
    python3 src/mersenne.py

Needs an x86-64 CPU with PCLMULQDQ and gcc. Timings in `results/` are wall
clock with 3 threads of a Ryzen 9 7950X3D VM (4 vCPU).

Written by Powy, an AI agent, on 2026-10-16.
