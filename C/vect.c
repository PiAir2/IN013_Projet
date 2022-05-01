#include <immintrin.h>
#include "polynome.h"
#include "vect.h"

void vect_mod_add(Uint *res1, Uint *tab1, Uint *tab2) {
    __m256i p = _mm256_set1_epi32(NB_P);
    __m256i a = _mm256_loadu_si256((__m256i *) tab1);
    __m256i b = _mm256_loadu_si256((__m256i *) tab2);
    __m256i x = _mm256_add_epi32(a, b);
    __m256i result = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
    // __m256i y = _mm256_sub_epi32(x, p);
    // __m256i mask = _mm256_cmpgt_epi32(p, x);
    // __m256i result = _mm256_blendv_epi8(y, x, mask);
    _mm256_storeu_si256((__m256i *) res1, result);
}

void vect_mod_sub(Uint *res1, Uint *tab1, Uint *tab2) {
    __m256i p = _mm256_set1_epi32(NB_P);
    __m256i a = _mm256_loadu_si256((__m256i *) tab1);
    __m256i b = _mm256_loadu_si256((__m256i *) tab2);
    __m256i x = _mm256_sub_epi32(a, b);
    __m256i result = _mm256_min_epu32(x, _mm256_add_epi32(x, p));
    // __m256i y = _mm256_add_epi32(x, p);
    // __m256i mask = _mm256_cmpgt_epi32(p, y);
    // __m256i result = _mm256_blendv_epi8(x, y, mask);
    _mm256_storeu_si256((__m256i *) res1, result);
}

void vect_mod_mult(Uint *res, Uint *tab1, Uint *tab2) {
    __m256i x = _mm256_loadu_si256((__m256i *) tab1);
    __m256i y = _mm256_loadu_si256((__m256i *) tab2);
    __m256i p = _mm256_set1_epi32(NB_P);
    __m256i q = _mm256_set1_epi32(NB_Q);

    __m256i a_low = _mm256_mul_epu32(x, y);
    __m256i b_low = _mm256_srli_epi64(a_low, NB_S);
    __m256i c_low = _mm256_srli_epi64(_mm256_mul_epu32(b_low, q), NB_T);

    __m256i a_high = _mm256_mul_epu32(_mm256_srli_si256(x, 4), _mm256_srli_si256(y, 4));
    __m256i b_high = _mm256_srli_epi64(a_high, NB_S);
    __m256i c_high = _mm256_srli_epi64(_mm256_mul_epu32(b_high, q), NB_T);

    __m256i d_low = _mm256_sub_epi64(a_low, _mm256_mul_epu32(c_low, p));
    __m256i d_high = _mm256_sub_epi64(a_high, _mm256_mul_epu32(c_high, p));
    __m256i d = _mm256_or_si256(d_low, _mm256_slli_si256(d_high, 4));
    
    __m256i result = _mm256_min_epu32(d, _mm256_sub_epi32(d, p));
    _mm256_storeu_si256((__m256i *) res, result);
}

void vect_mod_add_sub_eval(Uint *res_add, Uint *res_sub, Uint *tab1, Uint *tab2) {
    __m256i x, result;
    __m256i p = _mm256_set1_epi32(NB_P);
    __m256i a = _mm256_loadu_si256((__m256i *) tab1);
    __m256i b = _mm256_loadu_si256((__m256i *) tab2);

    // add
    x = _mm256_add_epi32(a, b);
    result = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
    _mm256_storeu_si256((__m256i *) res_add, result);

    // sub
    x = _mm256_sub_epi32(a, b);
    result = _mm256_min_epu32(x, _mm256_add_epi32(x, p));
    _mm256_storeu_si256((__m256i *) res_sub, result);
}

void vect_mod_mult_eval(Uint *res, Uint *tab1, Uint *tab2, Uint i, Uint pas, __m256i *u) {
    __m256i x = _mm256_loadu_si256((__m256i *) tab1);
    __m256i tmp = _mm256_set1_epi32(i);
    tmp = _mm256_add_epi32(tmp, *u);
    __m256i tmp2 = _mm256_set1_epi32(pas);
    tmp = _mm256_mullo_epi32(tmp , tmp2);
    __m256i y = _mm256_i32gather_epi32((const int *) tab2, tmp, 4);

    __m256i p = _mm256_set1_epi32(NB_P);
    __m256i q = _mm256_set1_epi32(NB_Q);

    __m256i a_low = _mm256_mul_epu32(x, y);
    __m256i b_low = _mm256_srli_epi64(a_low, NB_S);
    __m256i c_low = _mm256_srli_epi64(_mm256_mul_epu32(b_low, q), NB_T);

    __m256i a_high = _mm256_mul_epu32(_mm256_srli_si256(x, 4), _mm256_srli_si256(y, 4));
    __m256i b_high = _mm256_srli_epi64(a_high, NB_S);
    __m256i c_high = _mm256_srli_epi64(_mm256_mul_epu32(b_high, q), NB_T);

    __m256i d_low = _mm256_sub_epi64(a_low, _mm256_mul_epu32(c_low, p));
    __m256i d_high = _mm256_sub_epi64(a_high, _mm256_mul_epu32(c_high, p));
    __m256i d = _mm256_or_si256(d_low, _mm256_slli_si256(d_high, 4));
    
    __m256i result = _mm256_min_epu32(d, _mm256_sub_epi32(d, p));
    _mm256_storeu_si256((__m256i *) res, result);
}

// __m256i mod_x(Uint *res, __m256i x, Uint i, Uint p) {
//     __m256 float_p = _mm256_set1_ps(p);
//     __m256i int_p = _mm256_set1_epi32(p);

//     __m256i tmp;
//     __m256 x_div_p, float_x;
//     __m256i int_x_div_p;
//     float_x = _mm256_cvtepi32_ps(x); // float(x)
//     x_div_p = _mm256_div_ps(float_x, float_p); // float(x)/float(p)
//     x_div_p = _mm256_floor_ps(x_div_p); // lower((float(x)/float(p)))
//     int_x_div_p = _mm256_cvtps_epi32(x_div_p); // int(float(x)/float(p))
//     tmp = _mm256_mullo_epi32(int_x_div_p, int_p); // int(float(x)/float(p))*p
//     tmp = _mm256_sub_epi32(x, tmp); // x - int(float(x)*(1/float(p)))*p = x%p
//     _mm256_storeu_si256((__m256i *) &res[i], tmp);

//     return tmp;
// }