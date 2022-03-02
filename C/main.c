#include "polynome.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

int main() {
    srand(time(NULL));
    Poly P = gen_poly(20);
    Poly Q = gen_poly(20);
    Poly R1 = prod_poly_naif(P, Q);
    Poly R2 = prod_poly_karatsuba(P, Q);

    assert(compare_poly(R1, R2));

    liberer_poly(P);
    liberer_poly(Q);
    liberer_poly(R1);
    liberer_poly(R2);
    return 0;
}