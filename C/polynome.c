#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "polynome.h"

void afficher_poly(Poly P) {
    printf("%ld", (P.coeffs)[0]);
    for (int i = 1; i <= P.deg; i++) {
        printf(" + %ldx^%d", (P.coeffs)[i], i);
    }
    printf("\n");
}

Poly creer_poly(int deg) {
    Poly P;
    P.deg = deg;
    P.coeffs = (long *) malloc(sizeof(long)*(deg+1));
    return P;
}

Poly gen_poly(int deg) {
    Poly P;
    P.coeffs = (long *) malloc(sizeof(long)*(deg+1));
    P.deg = deg;
    for (int i = 0; i <= deg; i++) {
        (P.coeffs)[i] = rand() % NB_P;
    }
    return P;
}

Poly liberer_poly(Poly P) {
    free(P.coeffs);
}

int compare_poly(Poly P, Poly Q) {
    if (!(P.deg == Q.deg)) {
        return 0;
    }
    for (int i = 0; i <= P.deg; i++) {
        if ((P.coeffs)[i] != (Q.coeffs)[i]) {
            printf("Problème : C1 = %ld, C2 = %ld, pow = x^%d\n", (P.coeffs)[i], (Q.coeffs)[i], i);
            return 0;
        }
    }
    return 1;
}

Poly prod_poly_naif(Poly P, Poly Q) {
    Poly R;
    R.coeffs = (long *) malloc(sizeof(long)*(P.deg + Q.deg + 1));
    R.deg = P.deg + Q.deg;
    for (int i = 0; i <= P.deg; i++) {
        for (int j = 0; j <= Q.deg; j++) {
            (R.coeffs)[i + j] = ((R.coeffs)[i+j] + (P.coeffs)[i]*(Q.coeffs)[j]) % NB_P;
        }
    }
    return R;
}

Poly copy_poly(Poly P, int deb, int fin) {
    Poly R;
    R.coeffs = (long *) malloc(sizeof(long)*(fin - deb + 1));
    R.deg = fin - deb;
    for (int i = deb; i <= fin; i++) {
        (R.coeffs)[i-deb] = (P.coeffs)[i];
    }
    return R;
}

Poly somme_poly(Poly P, Poly Q, int deb) {
    Poly R = copy_poly(P, 0, P.deg);
    for (int i = 0; i <= Q.deg; i++) {
        (R.coeffs)[deb + i] = ((R.coeffs)[deb + i] + (Q.coeffs)[i]) % NB_P;
    }
    return R;
}

Poly oppose_poly(Poly P) {
    Poly R = copy_poly(P, 0, P.deg);
    for (int i = 0; i <= R.deg; i++) {
        (R.coeffs)[i] = -(R.coeffs)[i] + NB_P;
    }
    return R;
}

Poly prod_poly_karatsuba(Poly P, Poly Q) {
    if (P.deg <= 130) {
        return prod_poly_naif(P, Q);
    }
    int k = (int) ceil((P.deg+1) / 2.0);
    Poly P1 = copy_poly(P, 0, k-1);
    Poly P2 = copy_poly(P, k, P.deg);
    Poly Q1 = copy_poly(Q, 0, k-1);
    Poly Q2 = copy_poly(Q, k, Q.deg);
    Poly E1 = prod_poly_karatsuba(P1, Q1);
    Poly E2 = prod_poly_karatsuba(P2, Q2);
    Poly E3 = prod_poly_karatsuba(somme_poly(P1, P2, 0), somme_poly(Q1, Q2, 0));

    Poly R;
    R.deg = P.deg + Q.deg;
    R.coeffs = (long *) calloc(sizeof(long), R.deg + 1);

    R = somme_poly(R, E1, 0);
    R = somme_poly(R, E2, 2*((int) ((P.deg+2)/2.0)));
    R = somme_poly(R, somme_poly(E3, oppose_poly(somme_poly(E1, E2, 0)), 0), k);
    
    return R;
}

long modpow(long x, long n) {
    if (n == 0) {
        return 1;
    }
    if (n == 1) {
        return x;
    }
    long res;
    if (n % 2 == 0) {
        res = modpow(x, n/2);
        return (res*res) % NB_P;
    }
    res = modpow(x, n-1);
    return (res*x) % NB_P;
}

long horner(Poly P, long x) {
    long res = (P.coeffs)[P.deg];
    for (int i = P.deg-1; i >= 0; i--) {
        res = (res*x + (P.coeffs)[i]) % NB_P;
    }
    return res;
}

long *get_racines(long racine, int n) {
    long *racines = (long *) malloc(sizeof(long) * n);
    racines[0] = 1;
    for (int i = 1; i < n; i++) {
        racines[i] = racines[i-1]*racine;
    }
    return racines;
}

long *eval_P(Poly P, long *racines) {
    if (P.deg == 0) {
        long *tmp = (long *) malloc(sizeof(long));
        tmp[0] = (P.coeffs)[0];
        return tmp;
    }
    int k = (P.deg + 1)/2;
    Poly R0 = creer_poly(k-1);
    Poly R1 = creer_poly(k-1);
    
    long *racines_bis = (long *) malloc(sizeof(long) * k);
    for (int i = 0; i < k; i++) {
        (R0.coeffs)[i] = ((P.coeffs)[i] + (P.coeffs)[i + k]) % NB_P;
        (R1.coeffs)[i] = (((P.coeffs)[i] + (P.coeffs)[i + k]) * racines[i]) % NB_P;
        racines_bis[i] = racines[2*i];
    }
    long *r0 = eval_P(R0, racines_bis);
    long *r1 = eval_P(R1, racines_bis);
    long *res = (long *) malloc(sizeof(long) * 2*k);
    for (int i = 0; i < k; i++) {
        res[2*i] = r0[i];
        res[2*i+1] = r1[i];
    }
    return res;
}