#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include <limits.h>
#include "polynome.h"

clock_t temps_initial;
clock_t temps_final;
double temps_cpu;

//x % p = x - int(float(x)/float(p))*p

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

//Problème : fonction ne marche pas pour des nombres grands, même problème que mod_mult() :
//comment faire pour que les multiplications avec _mm256_mullo_epi retournent des entiers 64 bits ?
void vect_mult(Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    __m256 invP = _mm256_set1_ps(1.0f/p);
    __m256i intP = _mm256_set1_epi32(p);

    __m256i a, b, x, tmp;
    __m256 fx, fx_invP;
    __m256i int_fx_invP, int_fx_invP_p;
    for (Uint i = 0; i < taille; i += 8) {
        a = _mm256_loadu_si256((__m256i *) &tab1[i]);
        b = _mm256_loadu_si256((__m256i *) &tab2[i]);
        x = _mm256_mullo_epi32(a, b);
        fx = _mm256_cvtepi32_ps(x);
        fx_invP = _mm256_mul_ps(fx, invP);
        int_fx_invP = _mm256_cvtps_epi32(fx_invP);
        int_fx_invP_p = _mm256_mullo_epi32(int_fx_invP, intP);
        tmp = _mm256_sub_epi32(x, int_fx_invP_p);
        _mm256_storeu_si256((__m256i *) &res[i], tmp);
    }      
}


void add_norm(Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    for (int i = 0; i < taille; i++) {
        res[i] = mod_add(tab1[i], tab2[i], p);
    }
}

void sub_norm(Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    for (int i = 0; i < taille; i++) {
        res[i] = mod_sub(tab1[i], tab2[i], p);
    }
}

void mult_norm(Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    for (int i = 0; i < taille; i++) {
        res[i] = mod_mult(tab1[i], tab2[i], p);
    }
}

void test_add(Uint *res, Uint *resn, Uint *t1, Uint *t2, Uint taille, Uint p) {
    temps_initial = clock();
    for (int i = 0; i < 10000000; i++) {
        vect_add(res, t1, t2, taille, p);
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    for (int i = 0; i < taille; i++) printf("%d ", res[i]);
    printf("\nadd tmps vect : %f\n", temps_cpu);
    temps_initial = clock();
    for (int i = 0; i < 10000000; i++) {
        add_norm(resn, t1, t2, taille, p);
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    for (int i = 0; i < taille; i++) printf("%d ", resn[i]);
    printf("\nadd tmps norm : %f\n\n", temps_cpu);
}

void test_sub(Uint *res, Uint *resn, Uint *t1, Uint *t2, Uint taille, Uint p) {
    temps_initial = clock();
    for (int i = 0; i < 10000000; i++) {
        vect_sub(res, t1, t2, taille, p);
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    for (int i = 0; i < taille; i++) printf("%d ", res[i]);
    printf("\nsub tmps vect : %f\n", temps_cpu);
    temps_initial = clock();
    for (int i = 0; i < 10000000; i++) {
        sub_norm(resn, t1, t2, taille, p);
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    for (int i = 0; i < taille; i++) printf("%d ", resn[i]);
    printf("\nsub tmps norm : %f\n\n", temps_cpu);
}

void test_mult(Uint *res, Uint *resn, Uint *t1, Uint *t2, Uint taille) {
    temps_initial = clock();
    for (int i = 0; i < 10000000; i++) {
        vect_mult(res, t1, t2, taille, NB_P);
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    for (int i = 0; i < taille; i++) printf("%d ", res[i]);
    printf("\nmult tmps vect : %f\n", temps_cpu);
    temps_initial = clock();
    for (int i = 0; i < 10000000; i++) {
        mult_norm(resn, t1, t2, taille, NB_P);
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    for (int i = 0; i < taille; i++) printf("%d ", resn[i]);
    printf("\nmult tmps norm : %f\n\n", temps_cpu);
}

int main() {
    srand(time(NULL));

    int p = NB_P;
    //int coeff = 1073741824;
    int taille = 32;
    Uint *t1 = (Uint *) malloc(sizeof(Uint)*taille);
    Uint *t2 = (Uint *) malloc(sizeof(Uint)*taille);

    for (int i = 0; i < taille; i++) {
        t1[i] = rand()%p;
        t2[i] = rand()%p;
        //t1[i] = coeff;
        //t2[i] = coeff;
    }
    /*printf("tab1 : ");
    for (int i = 0; i < taille; i++) {
        printf("%d ", t1[i]);
    }
    printf("\ntab2 : ");
    for (int i = 0; i < taille; i++) {
        printf("%d ", t2[i]);
    }
    printf("\n\n");
    */

    Uint *res = (Uint *) malloc(sizeof(Uint)*taille);
    Uint *resn = (Uint *) malloc(sizeof(Uint)*taille);

    test_add(res, resn, t1, t2, taille, p);
    test_sub(res, resn, t1, t2, taille, p);
    //test_mult(res, resn, t1, t2, taille);

    printf("\n");
    return 0;
}
