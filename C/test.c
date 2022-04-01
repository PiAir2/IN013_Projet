#include "polynome.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

clock_t temps_initial;
clock_t temps_final;
double temps_cpu;

Poly test_naif(Poly P, Poly Q, int deg) {
    temps_initial = clock();
    Poly R = prod_poly_naif(P, Q);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    printf("Degré = %d, Naif : %f\n", deg, temps_cpu);
    return R;
}

Poly test_karatsuba(Poly P, Poly Q, int deg) {
    temps_initial = clock();
    Poly R = prod_poly_karatsuba(P, Q);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    printf("Degré = %d Karatsuba : %f\n", deg, temps_cpu);
    return R;
}

void compare_naif_karatsuba() {
    Poly P, Q, R1, R2;
    for (int i = 32; i < 100000; i *= 2) {
        P = gen_poly(i);
        Q = gen_poly(i);

        R1 = test_naif(P, Q, i);
        R2 = test_karatsuba(P, Q, i);
        //assert(compare_poly(R1, R2)); //test validité
    }
    liberer_poly(P);
    liberer_poly(Q);
    liberer_poly(R1);
    liberer_poly(R2);
}

void test_eval(int deg, Uint racine) {
    Poly P = gen_poly(deg);
    Poly P_cpy = copy_poly(P, 0, P.deg);
    temps_initial = clock();
    Uint *racines = get_racines(racine, P.deg+1);
    //Uint *res = eval(P, racines); ancien eval()
    Uint *tmp_coeffs = (Uint *) malloc(sizeof(Uint)*(deg+1));
    Uint *res = eval(P.coeffs, P.deg, tmp_coeffs, racines, 1);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    printf("Degré = %d : %f\n", P.deg, temps_cpu);

    
    for (int i = 0; i < P.deg+1; i++) {
        //printf("%d %d %d\n", res[i], horner(P_cpy, racines[i]), racines[i]);
        assert(res[i] == horner(P_cpy, racines[i]));
    }
    //afficher_poly(P);
    
    liberer_poly(P);
}

int main() {
    //srand(time(NULL));

    Uint racine = 2;
    Uint ordre_racine = (NB_P-1)/2;

    //compare_naif_karatsuba();

    int deg = 16777215;
    test_eval(deg, mod_pow(racine, ordre_racine/(deg+1)));
    
    printf("\n");
    return 0;
}
