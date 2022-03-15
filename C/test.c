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
    printf("Degré = %d Karatsuba : %f\n", deg,  temps_cpu);
    return R;
}

void compare_naif_karatsuba() {
    Poly P, Q, R1, R2;
    for (int i = 32; i < 60000; i += 1000) {
        P = gen_poly(i);
        Q = gen_poly(i);

        R1 = test_naif(P, Q, i);
        R2 = test_karatsuba(P, Q, i);
        assert(compare_poly(R1, R2));
    }
    liberer_poly(P);
    liberer_poly(Q);
    liberer_poly(R1);
    liberer_poly(R2);
}

int main() {
    srand(time(NULL));

    //compare_naif_karatsuba();

    //printf("%ld\n", modpow(2, (long) 2013265920/4));

    /*
    Poly P = gen_poly(3);
    afficher_poly(P);
    printf("%ld\n", horner(P, 3));
    */

    Poly P = gen_poly(256);
    afficher_poly(P);
    long racine = 2;
    long *res = eval_P(P, get_racines(racine, P.deg));
    for (int i = 0; i < 1; i++) {
        printf("%ld | ", res[i]);
    }
    printf("\n\n");
    
    printf("%ld\n", horner(P, 1));

    printf("\n");

    return 0;
}