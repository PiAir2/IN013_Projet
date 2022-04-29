#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include <limits.h>

#define NB_P 754974721
typedef unsigned int Uint;
// barrett algorithm with avx
//gcc -Wall -mavx2 -o c c.c

void vect_mult(Uint *res, Uint *tab1, Uint *tab2) {
    int s = 0;
    int t = 30;
    int q = 1;
    __m256i p = _mm256_set1_epi32(NB_P);
    __m256i ql = _mm256_set1_epi32(q);
    __m256i x = _mm256_loadu_si256((__m256i *) tab1);
    __m256i y = _mm256_loadu_si256((__m256i *) tab2);
    
    __m256i al = _mm256_mul_epu32(x, y);
    __m256i bl = _mm256_srli_epi64(al, s);
    __m256i cl = _mm256_srli_epi64(_mm256_mul_epi32(bl, ql), t);
    __m256i xh = _mm256_srli_si256(x, 2);
    __m256i yh = _mm256_srli_si256(y, 2);
    __m256i ah = _mm256_mul_epu32(xh, yh);
    __m256i bh = _mm256_srli_epi64(ah, s);
    __m256i ch = _mm256_srli_epi64(_mm256_mul_epi32(bh, ql), t);
    __m256i a = _mm256_mullo_epi32(x, y);
    //__m256i a = _mm256_blend_epi32(al, _mm256_slli_si256(ah, 4), 4+8+64+128);
    __m256i c = _mm256_or_si256(cl, _mm256_slli_si256(ch, 2));
    __m256i d = _mm256_sub_epi32(a, _mm256_mullo_epi32(c, p));
    __m256i result = _mm256_min_epu32(d, _mm256_sub_epi32(d, p));
    _mm256_storeu_si256((__m256i *) res, result);
}

// void vect_mult(Uint *res, Uint *tab1, Uint *tab2) {
//     int s = 0;
//     int t = 30;
//     __m256i p = _mm256_set1_epi32(NB_P);
//     __m256i zero = _mm256_set1_epi32(0);
//     __m256i ql = _mm256_set1_epi32(1);
//     __m256i x = _mm256_loadu_si256((__m256i *) tab1);
//     __m256i y = _mm256_loadu_si256((__m256i *) tab2);
    
//     __m256i xl = _mm256_unpacklo_epi32(x, zero);
//     __m256i xh = _mm256_unpackhi_epi32(x, zero);
//     __m256i yl = _mm256_unpacklo_epi32(y, zero);
//     __m256i yh = _mm256_unpackhi_epi32(y, zero);
//     __m256i al = _mm256_mul_epu32(xl, yl);
//     __m256i ah = _mm256_mul_epu32(xh, yh);
//     __m256i bl = _mm256_srli_epi64(al, s);
//     __m256i bh = _mm256_srli_epi64(ah, s);
//     __m256i cl = _mm256_srli_epi64(_mm256_mullo_epi32(bl, ql), t);
//     __m256i ch = _mm256_srli_epi64(_mm256_mullo_epi32(bh, ql), t);
//     __m256i c = _mm256_packus_epi32(cl, ch);
//     __m256i d = _mm256_sub_epi32(_mm256_mullo_epi32(x, y), _mm256_mullo_epi32(c, p));
//     __m256i result = _mm256_min_epu32(d, _mm256_sub_epi32(d, p));
//     _mm256_storeu_si256((__m256i *) res, result);
// }

// void vect_mult(Uint *res, Uint *tab1, Uint *tab2) {
//     __m256i zero = _mm256_set1_epi32(0);
//     __m256i a = _mm256_loadu_si256((__m256i *) tab1);
//     __m256i b = _mm256_loadu_si256((__m256i *) tab2);
//     __m256i a_low = _mm256_unpacklo_epi32(a, zero);
//     __m256i a_high = _mm256_unpackhi_epi32(a, zero);
//     __m256i b_low = _mm256_unpacklo_epi32(b, zero);
//     __m256i b_high = _mm256_unpackhi_epi32(b, zero);
//     __m256i x_low = _mm256_mul_epu32(a_low, b_low);
//     __m256i x_high = _mm256_mul_epu32(a_high, b_high);
//     __m256i result = _mm256_packs_epi32(x_low, x_high);

//     _mm256_storeu_si256((__m256i *) res, result);
// }

int main() { // 2 / 566qqchose
    Uint tab1[8] = {754974718, 2, 3, 4, 5, 6, 7, 8};
    Uint tab2[8] = {754974719, 10, 11, 12, 13, 14, 15, 16};
    Uint res[8];
    vect_mult(res, tab1, tab2);
    printf("%ld\n",  (754974718l*754974719l)%754974721);
    for (int i = 0; i < 8; i++) {
        printf("%d ", res[i]);
    }
    printf("\n");
    return 0;
}
