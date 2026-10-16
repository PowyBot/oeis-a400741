// a(n) of OEIS A400741 straight from the root side of the definition.
//
// f of degree 2n is primitive with f(x+1) = f(x) iff a root alpha of f is a
// primitive element of GF(2^(2n)) and alpha+1 is one of its conjugates
// alpha^(2^j). Applying the Frobenius twice gives alpha^(2^(2j)) = alpha, so
// j = n: the roots are exactly the primitive alpha on the affine line
//     L = { alpha : alpha^q + alpha = 1 },  q = 2^n,
// and each f owns 2n of them. So a(n) = #(primitive elements of L) / (2n).
//
// We walk all q points of L in Gray-code order (one XOR per step) and test
// each for primitivity. Nothing here shares code paths with orbit_count.c
// beyond the field multiplication.
//
// usage: line_count n [threads]
#include <pthread.h>
#include <string.h>

#include "gf2.h"

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static field F;
static uint64_t cofactors[64];  // (2^(2n) - 1) / p for each prime p
static int n_primes;
static uint64_t beta;           // some point of L
static uint64_t subfield_basis[64];  // GF(2)-basis of GF(q) inside GF(q^2)
static int n;

static int is_primitive(uint64_t a) {
    for (int i = 0; i < n_primes; i++)
        if (gf_pow(&F, a, cofactors[i]) == 1) return 0;
    return 1;
}

// trace(a) = a^q + a, the GF(2)-linear map whose kernel is GF(q) and whose
// fibre over 1 is L.
static uint64_t relative_trace(uint64_t a) {
    uint64_t b = a;
    for (int i = 0; i < n; i++) b = gf_mul(&F, b, b);
    return b ^ a;
}

// Gaussian elimination over GF(2) on the images of the monomial basis:
// yields a basis of the kernel (= GF(q)) and one preimage of 1 (= beta).
static void find_line(void) {
    int m = F.m;
    uint64_t image[64], preimage[64];
    for (int i = 0; i < m; i++) { image[i] = relative_trace(1ULL << i); preimage[i] = 1ULL << i; }
    int rank = 0, pivot_bit[64];
    for (int bit = m - 1; bit >= 0; bit--) {
        int r = -1;
        for (int i = rank; i < m; i++) if ((image[i] >> bit) & 1) { r = i; break; }
        if (r < 0) continue;
        uint64_t t = image[r]; image[r] = image[rank]; image[rank] = t;
        t = preimage[r]; preimage[r] = preimage[rank]; preimage[rank] = t;
        for (int i = 0; i < m; i++)
            if (i != rank && ((image[i] >> bit) & 1)) { image[i] ^= image[rank]; preimage[i] ^= preimage[rank]; }
        pivot_bit[rank++] = bit;
    }
    if (rank != n) { fprintf(stderr, "trace rank %d != %d\n", rank, n); exit(1); }
    for (int i = 0; i < n; i++) subfield_basis[i] = preimage[rank + i];
    // Rows are now reduced, so 1 is the sum of the rows whose pivot is bit 0.
    beta = 0;
    for (int i = 0; i < rank; i++) if (pivot_bit[i] == 0) beta = preimage[i];
    if (relative_trace(beta) != 1) { fprintf(stderr, "no point on L\n"); exit(1); }
}

typedef struct { uint64_t begin, end, count; } job;

static void *walk(void *arg) {
    job *jb = arg;
    uint64_t gray = jb->begin ^ (jb->begin >> 1);
    uint64_t a = beta;
    for (int i = 0; i < n; i++) if ((gray >> i) & 1) a ^= subfield_basis[i];
    uint64_t count = 0;
    for (uint64_t s = jb->begin; s < jb->end; s++) {
        count += is_primitive(a);
        // Gray code: step s -> s+1 flips the bit at the trailing zeros of s+1.
        a ^= subfield_basis[__builtin_ctzll(s + 1)];
    }
    jb->count = count;
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n [threads]\n", argv[0]); return 1; }
    n = atoi(argv[1]);
    int threads = argc > 2 ? atoi(argv[2]) : 1;
    if (n < 1 || n > 31) { fprintf(stderr, "n must be in 1..31\n"); return 1; }

    F = find_low_modulus(2 * n);
    uint64_t q = 1ULL << n, order = q * q - 1, primes[64];
    // q^2 - 1 = (q-1)(q+1) with coprime factors, so factor them separately.
    n_primes = prime_factors(q - 1, primes);
    n_primes += prime_factors(q + 1, primes + n_primes);
    qsort(primes, n_primes, sizeof primes[0], cmp_u64);  // small p first rejects most non-primitive points early
    for (int i = 0; i < n_primes; i++) cofactors[i] = order / primes[i];
    find_line();

    pthread_t tid[64];
    job jobs[64];
    for (int t = 0; t < threads; t++) {
        jobs[t].begin = q / threads * t;
        jobs[t].end = (t == threads - 1) ? q : q / threads * (t + 1);
        pthread_create(&tid[t], NULL, walk, &jobs[t]);
    }
    uint64_t total = 0;
    for (int t = 0; t < threads; t++) { pthread_join(tid[t], NULL); total += jobs[t].count; }
    if (total % (2 * n)) { fprintf(stderr, "count %llu not divisible by 2n\n", (unsigned long long)total); return 1; }
    printf("%d %llu\n", n, (unsigned long long)(total / (2 * n)));
    return 0;
}
