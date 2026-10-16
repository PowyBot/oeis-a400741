// a(n) of OEIS A400741 using only GF(q) arithmetic, q = 2^n.
//
// The OEIS formula: a(n) = #{u in GF(q^2) : ord(u) = q+1 and u + u^q is a
// primitive element of GF(q)} / (2n). Fix one u of order q+1 and write
// t = u + 1/u; the others are u^k with gcd(k, q+1) = 1, and
//     u^k + u^(kq) = u^k + u^(-k) = V_k(t),
// the Dickson/Lucas sequence V_0 = 0, V_1 = t, V_(k+1) = t V_k + V_(k-1).
// Squaring maps V_k to V_(2k), and 2^n = -1 (mod q+1), so k -> 2k generates a
// group of order 2n on the units mod q+1 that preserves primitivity. Every
// orbit has exactly 2n elements (2^j k = k forces 2^j = 1, i.e. 2n | j), so
//     a(n) = number of orbits {2^j k} whose V_k(t) is primitive in GF(q).
// We visit each orbit once, at its least member, and pay one ladder and one
// primitivity test per orbit instead of per k: ~2n times less field work
// than the PARI program in the OEIS entry.
//
// usage: orbit_count n [threads]
#include <pthread.h>

#include "gf2.h"

static field F;
static int n;
static uint64_t q, modulus_k;   // modulus_k = q + 1
static uint64_t t;              // u + 1/u for a fixed u of order q+1
static uint64_t cofactors[64];  // (q-1)/p for the primes p | q-1
static int n_primes;

// V_e(t) by the ladder (V_k, V_(k+1)) -> (V_2k, V_2k+1) or (V_2k+1, V_2k+2),
// using V_2k = V_k^2 and V_2k+1 = V_k V_(k+1) + t (characteristic 2).
static uint64_t dickson(uint64_t e, uint64_t s) {
    uint64_t v0 = 0, v1 = s;
    for (int i = 63 - __builtin_clzll(e | 1); i >= 0; i--) {
        uint64_t mid = gf_mul(&F, v0, v1) ^ s;
        if ((e >> i) & 1) { v0 = mid; v1 = gf_mul(&F, v1, v1); }
        else              { v1 = mid; v0 = gf_mul(&F, v0, v0); }
    }
    return v0;
}

// u has order exactly q+1 iff u^(q+1) = 1 and u^((q+1)/p) != 1 for p | q+1;
// u^e = 1 iff V_e(t) = 0, because u has odd order and u^e = u^-e then.
static uint64_t find_t(void) {
    uint64_t primes[64];
    int k = prime_factors(modulus_k, primes);
    for (uint64_t s = 1; s <= F.mask; s++) {
        if (dickson(modulus_k, s) != 0) continue;
        int ok = 1;
        for (int i = 0; i < k && ok; i++) ok = dickson(modulus_k / primes[i], s) != 0;
        if (ok) return s;
    }
    fprintf(stderr, "no element of order q+1\n");
    exit(1);
}

static int is_primitive(uint64_t a) {
    if (a == 0) return 0;
    for (int i = 0; i < n_primes; i++)
        if (gf_pow(&F, a, cofactors[i]) == 1) return 0;
    return 1;
}

static uint64_t gcd(uint64_t a, uint64_t b) {
    while (b) { uint64_t r = a % b; a = b; b = r; }
    return a;
}

// k is the least member of its orbit {2^j k mod (q+1)}. Usually decided in a
// few doublings: one wraps below k with probability about k/(q+1).
static int least_in_orbit(uint64_t k) {
    uint64_t x = k;
    for (int j = 1; j < 2 * n; j++) {
        x <<= 1;
        if (x >= modulus_k) x -= modulus_k;
        if (x < k) return 0;
    }
    return 1;
}

// Least orbit members cluster near 0 (the minimum of 2n spread-out residues
// is about q/(2n)), so threads take small chunks in turn instead of one
// contiguous range each.
#define CHUNK (1ULL << 16)
static uint64_t limit;
static uint64_t next_chunk;  // shared, advanced with __atomic_fetch_add

typedef struct { uint64_t orbits, primitive; } tally;

static void *scan(void *arg) {
    tally *out = arg;
    uint64_t orbits = 0, primitive = 0;
    for (;;) {
        uint64_t begin = __atomic_fetch_add(&next_chunk, CHUNK, __ATOMIC_RELAXED);
        if (begin >= limit) break;
        uint64_t end = begin + CHUNK < limit ? begin + CHUNK : limit;
        for (uint64_t k = begin; k < end; k++) {
            if (!least_in_orbit(k)) continue;
            if (gcd(k, modulus_k) != 1) continue;
            orbits++;
            primitive += is_primitive(dickson(k, t));
        }
    }
    out->orbits = orbits;
    out->primitive = primitive;
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [threads]\n", argv[0]); return 1; }
    n = atoi(argv[1]);
    int threads = argc > 2 ? atoi(argv[2]) : 1;
    if (n < 1 || n > 50) { fprintf(stderr, "n must be in 1..50\n"); return 1; }
    if (n == 1) { printf("%d 1\n", n); return 0; }  // GF(2) has no modulus x^1 + low with our fold

    F = find_low_modulus(n);
    q = 1ULL << n;
    modulus_k = q + 1;
    uint64_t primes[64];
    n_primes = prime_factors(q - 1, primes);
    for (int i = 0; i < n_primes; i++) cofactors[i] = (q - 1) / primes[i];
    t = find_t();

    // Every orbit contains -k along with k, so its least member is <= q/2.
    limit = q / 2 + 1;
    next_chunk = 1;
    pthread_t tid[64];
    tally tallies[64];
    for (int i = 0; i < threads; i++) pthread_create(&tid[i], NULL, scan, &tallies[i]);
    uint64_t orbits = 0, total = 0;
    for (int i = 0; i < threads; i++) {
        pthread_join(tid[i], NULL);
        orbits += tallies[i].orbits;
        total += tallies[i].primitive;
    }

    // Completeness check on the orbit walk: the units mod q+1 split into
    // exactly phi(q+1)/(2n) orbits.
    uint64_t phi = modulus_k, up[64];
    int nu = prime_factors(modulus_k, up);
    for (int i = 0; i < nu; i++) phi = phi / up[i] * (up[i] - 1);
    if (orbits * 2 * n != phi) {
        fprintf(stderr, "visited %llu orbits, expected %llu\n",
                (unsigned long long)orbits, (unsigned long long)(phi / (2 * n)));
        return 1;
    }
    printf("%d %llu\n", n, (unsigned long long)total);
    return 0;
}
