#include "polynome.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

int main() {
    srand(time(NULL));
    clock_t temps_initial;
    clock_t temps_final;
    double temps_cpu;

    Poly P;
    Poly Q;
    Poly R1;
    Poly R2;
    for (int i = 32; i < 100000; i*=2) {
        P = gen_poly(i);
        Q = gen_poly(i);
        /*
        temps_initial = clock();
        R1 = prod_poly_naif(P, Q);
        temps_final = clock();
        temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
        printf("%d Naif : %f\n", i, temps_cpu);
        */

        temps_initial = clock();
        R2 = prod_poly_karatsuba(P, Q);
        temps_final = clock();
        temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
        printf("%d Karatsuba : %f\n", i,  temps_cpu);

        /*if (compare_poly(R1, R2) == 0) {
            printf("RANG : %d\n", i);
            printf("COEFF 0 : %ld %ld\n", (P.coeffs)[0], (Q.coeffs)[0]);
            break;
        } */
    }

    liberer_poly(P);
    liberer_poly(Q);
    liberer_poly(R1);
    liberer_poly(R2);
    return 0;
}