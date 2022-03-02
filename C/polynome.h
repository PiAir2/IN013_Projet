#ifndef POLYNOME_H
#define POLYNOME_H

#define NB_P 1283

typedef struct _Poly {
    long *coeffs;
    int deg;
} Poly;

void afficher_poly(Poly P);
Poly gen_poly(int deg);
Poly liberer_poly(Poly P);
int compare_poly(Poly P, Poly Q);
Poly prod_poly_naif(Poly P, Poly Q);
Poly copy_poly(Poly P, int deb, int fin);
Poly poly_somme(Poly P, Poly Q, int deb);
Poly oppose_poly(Poly P);
Poly prod_poly_karatsuba(Poly P, Poly Q);

#endif