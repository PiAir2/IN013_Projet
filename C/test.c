#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "polynome.h"

clock_t temps_initial;
clock_t temps_final;
double temps_cpu;

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

void test_eval(int deg, Uint racine) {
    Poly P = gen_poly(deg);
    Poly P_cpy = copy_poly(P, 0, P.deg);
    Uint *racines = get_racines(racine, P.deg+1);
    //Uint *res = eval(P, racines); ancien eval()
    temps_initial = clock();
    Uint *tmp_coeffs = (Uint *) malloc(sizeof(Uint)*(deg+1));
    Uint *res = eval(P.coeffs, P.deg+1, tmp_coeffs, racines, 1);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    printf("Degré = %d | Temps eval() = %f\n", P.deg, temps_cpu);

    // for (int i = 0; i < P.deg+1; i++) {
    //     // printf("%d %d %d\n", res[i], horner(P_cpy, racines[i]), racines[i]);
    //     assert(res[i] == horner(P_cpy, racines[i]));
    // }
    // afficher_poly(P);
    
    liberer_poly(P);
}

void test_vect_eval(int deg, Uint racine) {
    Poly P = gen_poly(deg);
    Poly P_cpy = copy_poly(P, 0, P.deg);
    Uint *racines = get_racines(racine, P.deg+1);
    //Uint *res = eval(P, racines); ancien eval()
    temps_initial = clock();
    Uint *tmp_coeffs = (Uint *) malloc(sizeof(Uint)*(deg+1));
    Uint *tmp_sub = (Uint *) malloc(sizeof(Uint)*(deg+1));
    Uint *res = vect_eval(P.coeffs, P.deg+1, tmp_coeffs, racines, 1, tmp_sub);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    printf("Degré = %d | Temps vect_eval() = %f\n", P.deg, temps_cpu);

    // for (int i = 0; i < P.deg+1; i++) {
    //     // printf("%d %d %d\n", res[i], horner(P_cpy, racines[i]), racines[i]);
    //     assert(res[i] == horner(P_cpy, racines[i]));
    // }
    // afficher_poly(P);
    
    liberer_poly(P);
}

int main() {
    //srand(time(NULL));

    //Uint racine = 2;
    //Uint ordre_racine = (NB_P-1)/2;
    Uint racine = 11;
    Uint ordre_racine = NB_P-1;

    //compare_naif_karatsuba();

    // 65535, 16777215;
    int deg = 16777215; //2^k - 1
    Uint rac = mod_pow(racine, ordre_racine/(deg+1));
    for (int i = 0; i < 10; i++) {
        test_eval(deg, rac);
        test_vect_eval(deg, rac);
    }
    
    printf("\n");
    return 0;
}
