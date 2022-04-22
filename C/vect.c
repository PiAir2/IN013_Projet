#include <immintrin.h>
#include "vect.h"


__m256i mod_x(Uint *res, __m256i x, Uint i, Uint p) {
    __m256 float_p = _mm256_set1_ps(p);
    __m256i int_p = _mm256_set1_epi32(p);

    __m256i tmp;
    __m256 x_div_p, float_x;
    __m256i int_x_div_p;
    float_x = _mm256_cvtepi32_ps(x); // float(x)
    x_div_p = _mm256_div_ps(float_x, float_p); // float(x)/float(p)
    x_div_p = _mm256_floor_ps(x_div_p); // lower((float(x)/float(p)))
    int_x_div_p = _mm256_cvtps_epi32(x_div_p); // int(float(x)/float(p))
    tmp = _mm256_mullo_epi32(int_x_div_p, int_p); // int(float(x)/float(p))*p
    tmp = _mm256_sub_epi32(x, tmp); // x - int(float(x)*(1/float(p)))*p = x%p
    _mm256_storeu_si256((__m256i *) &res[i], tmp);

    return tmp;
}


void vect_add(Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    __m256i a, b, x;
    for (Uint i = 0; i < taille; i += 8) {
        a = _mm256_loadu_si256((__m256i *) &tab1[i]);
        b = _mm256_loadu_si256((__m256i *) &tab2[i]);
        x = _mm256_add_epi32(a, b); // x = a + b
        x = mod_x(res, x, i, p);
    }
}


void vect_sub(Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    __m256i a, b, x;
    for (Uint i = 0; i < taille; i += 8) {
        a = _mm256_loadu_si256((__m256i *) &tab1[i]);
        b = _mm256_loadu_si256((__m256i *) &tab2[i]);
        x = _mm256_sub_epi32(a, b); // x = a - b
        mod_x(res, x, i, p);
    }
}

void vect_add_sub(Uint *res1, Uint *res2, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    __m256i a, b, x;
    for (Uint i = 0; i < taille; i += 8) {
        a = _mm256_loadu_si256((__m256i *) &tab1[i]);
        b = _mm256_loadu_si256((__m256i *) &tab2[i]);
        x = _mm256_add_epi32(a, b); // x = a + b
        mod_x(res1, x, i, p);
        x = _mm256_sub_epi32(a, b); // x = a - b
        mod_x(res2, x, i, p);
    }
}


//void vect_add((Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
//    __m256i a, b, x;
//    for (Uint i = 0; i < taille; i += 8) {
//        a = _mm256_loadu_si256((__m256i *) &tab1[i]);
//        b = _mm256_loadu_si256((__m256i *) &tab2[i]);
//        x = _mm256_add_epi32(a, b); // x = a + b
//        
//    }
//}
