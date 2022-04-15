#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include <limits.h>
#include "polynome.h"

//gcc -Wall -mavx2 -o c c.c
int main() {
    int i = 0;
    int p = 23;

    __m256 invP = _mm256_set1_ps(1.0f/p);
    __m256i intP = _mm256_set1_epi32(p);

    __m256i a, b, x, tmp;
    __m256 fx, fx_invP;
    __m256i int_fx_invP, int_fx_invP_p;

    Uint tab1[8] = {2147483647, 0, 10, 0, 10, 0, 10, 0};
    Uint tab2[8] = {200, 0, 10, 0, 10, 0, 10, 0};
    Uint res[8];
    a = _mm256_loadu_si256((__m256i *) &tab1[i]);
    b = _mm256_loadu_si256((__m256i *) &tab2[i]);
    x = _mm256_mul_epu32(a, b);
    fx = _mm256_cvtepi32_ps(x);
    fx_invP = _mm256_mul_ps(fx, invP);
    int_fx_invP = _mm256_cvtps_epi32(fx_invP);
    int_fx_invP_p = _mm256_mullo_epi32(int_fx_invP, intP);
    tmp = _mm256_sub_epi32(x, int_fx_invP_p);
    _mm256_storeu_si256((__m256i *) &res[i], tmp);
    printf("%u %u %u %u %u %u %u %u\n", res[0], res[1], res[2], res[3], res[4], res[5], res[6], res[7]);
    return 0;
}
