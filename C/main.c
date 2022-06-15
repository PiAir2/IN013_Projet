#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include <limits.h>
#include "polynome.h"
#include "vect.h"

clock_t temps_initial;
clock_t temps_final;
double temps_cpu;

void add_norm(Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    for (int i = 0; i < taille; i++) {
        res[i] = mod_add(tab1[i], tab2[i]);
    }
}
void sub_norm(Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    for (int i = 0; i < taille; i++) {
        res[i] = mod_sub(tab1[i], tab2[i]);
    }
}
void mult_norm(Uint *res, Uint *tab1, Uint *tab2, Uint taille, Uint p) {
    for (int i = 0; i < taille; i++) {
        res[i] = mod_mult(tab1[i], tab2[i]);
    }
}


void test_add(Uint *res, Uint *resn, Uint *t1, Uint *t2, Uint taille, Uint p, int repet) {
    temps_initial = clock();
    for (int i = 0; i < repet; i++) {
        for (int j = 0; j < taille; j+=8) {
            vect_mod_add(&res[j], &t1[j], &t2[j]);
        }
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    //for (int i = 0; i < taille; i++) printf("%d ", res[i]);
    printf("\nadd tmps vect : %f\n", temps_cpu);
    temps_initial = clock();
    for (int i = 0; i < repet; i++) {
        add_norm(resn, t1, t2, taille, p);
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    //for (int i = 0; i < taille; i++) printf("%d ", resn[i]);

    // for (int i = 0; i < taille; i++) assert(res[i] == resn[i]);
    printf("\nadd tmps norm : %f\n\n", temps_cpu);
}

void test_sub(Uint *res, Uint *resn, Uint *t1, Uint *t2, Uint taille, Uint p, int repet) {
    temps_initial = clock();
    for (int i = 0; i < repet; i++) {
        for (int j = 0; j < taille; j+=8) {
            vect_mod_sub(&res[j], &t1[j], &t2[j]);
        }
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    //for (int i = 0; i < taille; i++) printf("%d ", res[i]);
    printf("\nsub tmps vect : %f\n", temps_cpu);
    temps_initial = clock();
    for (int i = 0; i < repet; i++) {
        sub_norm(resn, t1, t2, taille, p);
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    //for (int i = 0; i < taille; i++) printf("%d ", resn[i]);

    // for (int i = 0; i < taille; i++) assert(res[i] == resn[i]);
    printf("\nsub tmps norm : %f\n\n", temps_cpu);
}

void test_mult(Uint *res, Uint *resn, Uint *t1, Uint *t2, Uint taille, int repet) {
    temps_initial = clock();
    for (int i = 0; i < repet; i++) {
        for (int j = 0; j < taille; j+=8) {
            vect_mod_mult(&res[j], &t1[j], &t2[j]);
        }
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    //for (int i = 0; i < taille; i++) printf("%d ", res[i]);
    printf("\nmult tmps vect : %f\n", temps_cpu);
    temps_initial = clock();
    for (int i = 0; i < repet; i++) {
        mult_norm(resn, t1, t2, taille, NB_P);
    }
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    //for (int i = 0; i < taille; i++) printf("%d ", resn[i]);

    // for (int i = 0; i < taille; i++) assert(res[i] == resn[i]);
    printf("\nmult tmps norm : %f\n\n", temps_cpu);
}

int main() {
    srand(time(NULL));

    int p = NB_P;
    //int coeff = 1073741824;
    int taille = 80000;
    Uint *t1 = (Uint *) malloc(sizeof(Uint)*taille);
    Uint *t2 = (Uint *) malloc(sizeof(Uint)*taille);

    for (int i = 0; i < taille; i++) {
        t1[i] = rand()%p;
        t2[i] = rand()%p;
        // t1[i] = coeff;
        // t2[i] = coeff;
    }

    // printf("tab1 : ");
    // for (int i = 0; i < taille; i++) printf("%d ", t1[i]);
    // printf("\ntab2 : ");
    // for (int i = 0; i < taille; i++) printf("%d ", t2[i]);
    // printf("\n\n");
    
    Uint *res = (Uint *) malloc(sizeof(Uint)*taille);
    Uint *resn = (Uint *) malloc(sizeof(Uint)*taille);

    int repet = 1;
    test_add(res, resn, t1, t2, taille, p, repet); // echauffement du CPU.

    test_add(res, resn, t1, t2, taille, p, repet);
    test_sub(res, resn, t1, t2, taille, p, repet);
    test_mult(res, resn, t1, t2, taille, repet);

    free(t1);
    free(t2);
    free(res);
    free(resn);
    printf("\n");
    return 0;
}
