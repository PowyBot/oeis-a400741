// Arithmetic in GF(2^m), m <= 62, with elements packed in a uint64_t
// (bit i = coefficient of x^i), plus the integer helpers both counters share.
#pragma once
#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;

typedef struct {
    int m;
    // The modulus is x^m + low. Reduction folds the high half twice, which is
    // exact only while 2 deg(low) - 2 < m; find_low_modulus() guarantees it.
    uint64_t low;
    uint64_t mask;
} field;

static inline u128 clmul(uint64_t a, uint64_t b) {
    __m128i p = _mm_clmulepi64_si128(_mm_cvtsi64_si128((long long)a),
                                     _mm_cvtsi64_si128((long long)b), 0);
    uint64_t lo = (uint64_t)_mm_cvtsi128_si64(p);
    uint64_t hi = (uint64_t)_mm_extract_epi64(p, 1);
    return ((u128)hi << 64) | lo;
}

static inline uint64_t gf_reduce(const field *f, u128 p) {
    u128 high = p >> f->m;
    p = (p & f->mask) ^ clmul((uint64_t)high, f->low);
    high = p >> f->m;
    p = (p & f->mask) ^ clmul((uint64_t)high, f->low);
    return (uint64_t)p;
}

static inline uint64_t gf_mul(const field *f, uint64_t a, uint64_t b) {
    return gf_reduce(f, clmul(a, b));
}

static inline uint64_t gf_pow(const field *f, uint64_t a, uint64_t e) {
    uint64_t r = 1;
    for (int i = 63 - __builtin_clzll(e | 1); i >= 0; i--) {
        r = gf_mul(f, r, r);
        if ((e >> i) & 1) r = gf_mul(f, r, a);
    }
    return e ? r : 1;
}

static int poly_deg(uint64_t a) { return a ? 63 - __builtin_clzll(a) : -1; }

// gcd of two GF(2)[x] polynomials that fit in 64 bits.
static uint64_t poly_gcd(uint64_t a, uint64_t b) {
    while (b) {
        while (poly_deg(a) >= poly_deg(b)) a ^= b << (poly_deg(a) - poly_deg(b));
        uint64_t t = a; a = b; b = t;
    }
    return a;
}

// Ben-Or's test: x^m + low is irreducible iff gcd(x^(2^i) - x, f) = 1 for
// i = 1..m/2.  M. Ben-Or, "Probabilistic algorithms in finite fields",
// FOCS 1981, doi:10.1109/SFCS.1981.37.
static int is_irreducible(const field *f) {
    if (!(f->low & 1)) return 0;
    uint64_t modulus = (1ULL << f->m) | f->low;
    uint64_t xp = 2;  // x
    for (int i = 1; i <= f->m / 2; i++) {
        xp = gf_mul(f, xp, xp);
        if (poly_gcd(modulus, xp ^ 2) != 1) return 0;
    }
    return 1;
}

// The irreducible x^m + low with the numerically smallest low part.
static field find_low_modulus(int m) {
    if (m < 2 || m > 62) { fprintf(stderr, "degree %d out of range\n", m); exit(1); }
    field f = {m, 0, (1ULL << m) - 1};
    for (f.low = 1;; f.low += 2) {
        if (2 * poly_deg(f.low) - 2 >= m) { fprintf(stderr, "no low modulus\n"); exit(1); }
        if (is_irreducible(&f)) return f;
    }
}

// Distinct prime factors by trial division: at most 2^(n/2) steps for 2^n +- 1,
// negligible next to the counts for every n we run.
static int prime_factors(uint64_t x, uint64_t *out) {
    int k = 0;
    for (uint64_t d = 2; d * d <= x; d += (d == 2) ? 1 : 2) {
        if (x % d == 0) {
            out[k++] = d;
            while (x % d == 0) x /= d;
        }
    }
    if (x > 1) out[k++] = x;
    return k;
}
