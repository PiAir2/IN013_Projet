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
//      = x - int(float(x)*(1/float(p)))*p

//Problème : retourne parfois le résultat négatif du modulo pour certains éléments.
void vect_add(Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
                       //__mi256i *tab1, __mi256i *tab2 pour avoir à cast juste à l'appel, plus rapide ?

    __m256 invP = _mm256_set1_ps(1.0f/p);
    __m256i intP = _mm256_set1_epi32(p);

    __m256i a, b, x, tmp;
    __m256 fx, fx_invP;
    __m256i int_fx_invP, int_fx_invP_p;
    for (Uint i = 0; i < taille; i += 8) {
        a = _mm256_loadu_si256((__m256i *) &tab1[i]);
        b = _mm256_loadu_si256((__m256i *) &tab2[i]);
        x = _mm256_add_epi32(a, b); // x = a + b
        fx = _mm256_cvtepi32_ps(x); // float(x)
        fx_invP = _mm256_mul_ps(fx, invP); // float(x)*(1/float(p))
        int_fx_invP = _mm256_cvtps_epi32(fx_invP); // int(float(x)*(1/float(p)))
        int_fx_invP_p = _mm256_mullo_epi32(int_fx_invP, intP); // int(float(x)*(1/float(p)))*p
        tmp = _mm256_sub_epi32(x, int_fx_invP_p); // x - int(float(x)*(1/float(p)))*p = (a+b)%p
        _mm256_storeu_si256((__m256i *) &res[i], tmp);

        /* //première version : environ 35% plus rapide avec p = NB_P et taille = 80
        a = _mm256_loadu_si256((__m256i *) &tab1[i]);
        b = _mm256_loadu_si256((__m256i *) &tab2[i]);
        tmp = _mm256_add_epi32(a, b);
        _mm256_storeu_si256((__m256i *) &res[i], tmp);
        for (Uint j = i; j < i+8; j++) { //Comment gérer les modulos efficacement ?
            if (res[j] >= p) {
                res[j] = res[j] - p;
            }
        }*/
    }
}

//Problème : marche seulement si tous les éléments de tab1 >= tous les éléments de tab2,
//sinon retourne le résulat négatif du modulo quand tab1[i] < tab2[i]
void vect_sub(/*int*/ Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    __m256 invP = _mm256_set1_ps(1.0f/p);
    __m256i intP = _mm256_set1_epi32(p);

    __m256i a, b, x, tmp;
    __m256 fx, fx_invP;
    __m256i int_fx_invP, int_fx_invP_p;
    for (Uint i = 0; i < taille; i += 8) {
        a = _mm256_loadu_si256((__m256i *) &tab1[i]);
        b = _mm256_loadu_si256((__m256i *) &tab2[i]);
        x = _mm256_sub_epi32(a, b); // x = a - b
        fx = _mm256_cvtepi32_ps(x); // float(x)
        fx_invP = _mm256_mul_ps(fx, invP); // float(x)*(1/float(p))
        int_fx_invP = _mm256_cvtps_epi32(fx_invP); // int(float(x)*(1/float(p)))
        int_fx_invP_p = _mm256_mullo_epi32(int_fx_invP, intP); // int(float(x)*(1/float(p)))*p
        tmp = _mm256_sub_epi32(x, int_fx_invP_p); // x - int(float(x)*(1/float(p)))*p = (a+b)%p
        _mm256_storeu_si256((__m256i *) &res[i], tmp);

        /* //première version : environ 40% plus rapide avec p = NB_P et taille = 80
        a = _mm256_loadu_si256((__m256i *) &tab1[i]);
        b = _mm256_loadu_si256((__m256i *) &tab2[i]);
        tmp = _mm256_sub_epi32(a, b);
        _mm256_storeu_si256((__m256i *) &res[i], tmp);
        for (Uint j = i; j < i+8; j++) {
            if (res[j] < 0) {
                res[j] = res[j] + p;
            }
        }*/
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
        vect_sub(/*(int *)*/ res, t1, t2, taille, p);
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

    int p = 23;
    int taille = 16;
    Uint *t1 = (Uint *) malloc(sizeof(Uint)*taille);
    Uint *t2 = (Uint *) malloc(sizeof(Uint)*taille);

    for (int i = 0; i < taille; i++) {
        t1[i] = rand()%p;
        t2[i] = rand()%p;
    }
    printf("tab1 : ");
    for (int i = 0; i < taille; i++) {
        printf("%d ", t1[i]);
    }
    printf("\ntab2 : ");
    for (int i = 0; i < taille; i++) {
        printf("%d ", t2[i]);
    }
    printf("\n\n");

    Uint *res = (Uint *) malloc(sizeof(Uint)*taille);
    Uint *resn = (Uint *) malloc(sizeof(Uint)*taille);

    test_add(res, resn, t1, t2, taille, p);
    test_sub(res, resn, t1, t2, taille, p);
    //test_mult(res, resn, t1, t2, taille);

    printf("\n");
    return 0;
}