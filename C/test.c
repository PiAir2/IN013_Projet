#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include "polynome.h"

clock_t temps_initial;
clock_t temps_final;
double temps_cpu;
double temps_tot_eval_malloc = 0;
double temps_tot_eval = 0;
double temps_tot_vect_eval = 0;

Poly test_naif(Poly P, Poly Q, int deg, double *temps) {
    temps_initial = clock();
    Poly R = prod_poly_naif(P, Q);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    *temps = temps_cpu;
    printf("Degré = %d, Naif : %f\n", deg, temps_cpu);
    return R;
}

Poly test_karatsuba(Poly P, Poly Q, int deg, double *temps) {
    temps_initial = clock();
    Poly R = prod_poly_karatsuba(P, Q);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    *temps = temps_cpu;
    printf("Degré = %d Karatsuba : %f\n", deg, temps_cpu);
    return R;
}

void compare_naif_karatsuba() {
    FILE *f = fopen("./gnuplot/temps/temps_naif_kara.txt", "w");
    Poly P, Q, R1, R2;
    double temps_naif, temps_karatsuba;
    for (int i = 32; i < 100000; i *= 2) {
        P = gen_poly(i);
        Q = gen_poly(i);

        R1 = test_naif(P, Q, i, &temps_naif);
        R2 = test_karatsuba(P, Q, i, &temps_karatsuba);

        fprintf(f, "%d %.10f %.10f\n", i, temps_naif, temps_karatsuba);
        //assert(compare_poly(R1, R2));
    }
    liberer_poly(P);
    liberer_poly(Q);
    liberer_poly(R1);
    liberer_poly(R2);
    fclose(f);
}

void verif(Uint *res, Poly P, Uint *racines) {
    for (int i = 0; i < P.deg+1; i++) {
        int horner_x = horner(P, racines[i]);
        // if (res[i] != horner_x) printf("i = %d, vect_eval[i] = %d, horner[i] = %d\n", i, res[i], horner_x);
        assert(res[i] == horner_x);
    }
    // afficher_poly(P);
}

void test_eval_malloc(int deg, Uint racine, int aff, int v) {
    Poly P = gen_poly(deg);
    Poly P_cpy = copy_poly(P, 0, P.deg);
    Uint *racines = get_racines(racine, P.deg+1);
    temps_initial = clock();
    Uint *res = eval_malloc(P, racines);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    temps_tot_eval_malloc += temps_cpu;

    if (aff == 1) printf("Degré = %d | Temps eval() = %f\n", P.deg, temps_cpu);
    if (v == 1) verif(res, P_cpy, racines);
    
    liberer_poly(P);
}

void test_eval(int deg, Uint racine, int aff, int v) {
    Poly P = gen_poly(deg);
    Poly P_cpy = copy_poly(P, 0, P.deg);
    Uint *racines = get_racines(racine, P.deg+1);
    temps_initial = clock();
    Uint *tmp_coeffs = (Uint *) malloc(sizeof(Uint)*(deg+1));
    Uint *res = eval(P.coeffs, P.deg+1, tmp_coeffs, racines, 1);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    temps_tot_eval += temps_cpu;

    if (aff == 1) printf("Degré = %d | Temps eval() = %f\n", P.deg, temps_cpu);
    if (v == 1) verif(res, P_cpy, racines);

    liberer_poly(P);
}

void test_vect_eval(Uint deg, Uint racine, int aff, int v) {
    Poly P = gen_poly(deg);
    Poly P_cpy = copy_poly(P, 0, P.deg);
    Uint *racines = get_racines(racine, P.deg+1);
    //Uint *res = eval(P, racines); ancien eval()
    temps_initial = clock();
    Uint *tmp_coeffs = (Uint *) malloc(sizeof(Uint)*(deg+1));
    Uint *tmp_sub = (Uint *) malloc(sizeof(Uint)*(deg+1));
    __m256i u = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    Uint *res = vect_eval(P.coeffs, P.deg+1, tmp_coeffs, racines, 1, tmp_sub, &u);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    temps_tot_vect_eval += temps_cpu;

    if (aff == 1) printf("Degré = %d | Temps vect_eval() = %f\n", P.deg, temps_cpu);
    if (v == 1) verif(res, P_cpy, racines);

    liberer_poly(P);
}

int main() {
    srand(time(NULL));

    // Uint racine = 2;
    // Uint ordre_racine = (NB_P-1)/2;
    Uint racine = 11;
    Uint ordre_racine = NB_P-1;

    //compare_naif_karatsuba();

    // 65535, 16777215, 33554431, 67108863;
    Uint deg = 16777215; //2^k - 1
    Uint rac = mod_pow(racine, ordre_racine/(deg+1));
    
    // for (int i = 0; i < 10; i++) {
    //     test_eval(deg, rac);
    //     test_vect_eval(deg, rac);
    // }

    Uint nb_tours = 2;

    for (Uint i = 32768; i <= pow(2, 24); i = i*2) {
        Uint deg = i-1;
        Uint rac = mod_pow(racine, ordre_racine/(deg+1));
        for (Uint j = 0; j < nb_tours; j++) {
            test_eval_malloc(deg, rac, 0, 0);
            test_eval(deg, rac, 0, 0);
            test_vect_eval(deg, rac, 0, 0);
        }
        printf(">>>> Degré = %d <<<<\n", deg);
        printf("=====> Temps moyen eval_malloc : %f\n", temps_tot_eval_malloc/nb_tours);
        printf("=====> Temps moyen eval : %f\n", temps_tot_eval/nb_tours);
        printf("=====> Temps moyen vect_eval : %f\n", temps_tot_vect_eval/nb_tours);
    }
    
    printf("\n");
    return 0;
}
