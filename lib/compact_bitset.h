#ifndef COMPACT_BITSET_H
#define COMPACT_BITSET_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef _MSC_VER
    #include <intrin.h>
    #define popcount64(x) __popcnt64(x)
    #define ctz64(x) _tzcnt_u64(x)
#else
    #define popcount64(x) __builtin_popcountll(x)
    #define ctz64(x) __builtin_ctzl(x)
#endif

#define WORDSZ 64
#define ALIGNMENT 32

typedef struct {
    uint64_t* bits;
    int nbits;
    int nwords;
} CompactBitset;

static uint64_t* cb_aligned_alloc(size_t bytes) {
#ifdef _WIN32
    return (uint64_t*)_aligned_malloc(bytes, ALIGNMENT);
#else
    void* ptr = NULL;
    if (posix_memalign(&ptr, ALIGNMENT, bytes) != 0)
        return NULL;
    return (uint64_t*)ptr;
#endif
}

static void cb_aligned_free(uint64_t* ptr) {
#ifdef _WIN32
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

static CompactBitset* cb_create(int nbits) {
    CompactBitset* cb = (CompactBitset*)malloc(sizeof(CompactBitset));
    if (!cb) return NULL;

    cb->nbits = nbits;
    cb->nwords = (nbits + WORDSZ - 1) / WORDSZ;

    size_t total_bytes = cb->nwords * sizeof(uint64_t);
    cb->bits = cb_aligned_alloc(total_bytes);
    if (!cb->bits) {
        free(cb);
        return NULL;
    }

    memset(cb->bits, 0, total_bytes);
    return cb;
}

static void cb_destroy(CompactBitset* cb) {
    if (!cb) return;
    cb_aligned_free(cb->bits);
    free(cb);
}

static void cb_set_bit(CompactBitset* cb, int bit) {
    if (bit < 0 || bit >= cb->nbits) {
        fprintf(stderr, "Error: bit index %d out of range [0, %d)\n", bit, cb->nbits);
        exit(EXIT_FAILURE);
    }
    cb->bits[bit / WORDSZ] |= ((uint64_t)1 << (bit % WORDSZ));
}

// static void cb_clear_bit(CompactBitset* cb, int bit) {
//     if (bit < 0 || bit >= cb->nbits) return;
//     cb->bits[bit / WORDSZ] &= ~((uint64_t)1 << (bit % WORDSZ));
// }

static bool cb_get_bit(const CompactBitset* cb, int bit) {
    if (bit < 0 || bit >= cb->nbits) {
        fprintf(stderr, "Error: bit index %d out of range [0, %d)\n", bit, cb->nbits);
        exit(EXIT_FAILURE);
    }
    return (cb->bits[bit / WORDSZ] & ((uint64_t)1 << (bit % WORDSZ))) != 0;
}

// static int cb_count_ones(const CompactBitset* cb) {
//     int sum = 0;
//     for (int i = 0; i < cb->nwords; ++i) {
//         sum += popcount64(cb->bits[i]);
//     }
//     return sum;
// }

// static void cb_print(const CompactBitset* cb) {
//     for (int i = 0; i < cb->nbits; ++i) {
//         printf("%d", cb_get_bit(cb, i));
//     }
//     printf("\n");
// }

// static void cb_print_pos(const CompactBitset* cb) {
//     for (int i = 0; i < cb->nbits; ++i) {
//         if (cb_get_bit(cb, i)) {
//             printf("%d ", i);
//         }
//     }
//     printf("\n");
// }

// static void cb_print_pos_to_file(const CompactBitset* cb, FILE* fout) {
//     for (int i = 0; i < cb->nbits; ++i) {
//         if (cb_get_bit(cb, i)) {
//             fprintf(fout, "%d ", i);
//         }
//     }
//     fprintf(fout, "\n");
// }

#endif // COMPACT_BITSET_H
