#ifndef POLYNOME_H
#define POLYNOME_H

#define NB_P 2013265921

typedef struct _Poly {
    unsigned int *coeffs;
    int deg;
} Poly;

void afficher_poly(Poly P);
Poly gen_poly(unsigned int deg);
Poly liberer_poly(Poly P);
int compare_poly(Poly P, Poly Q);
Poly prod_poly_naif(Poly P, Poly Q);
Poly copy_poly(Poly P, int deb, int fin);
Poly poly_somme(Poly P, Poly Q, int deb);
Poly oppose_poly(Poly P);
Poly prod_poly_karatsuba(Poly P, Poly Q);
unsigned int mod_pow(unsigned int x, unsigned int n);
unsigned int horner(Poly P, unsigned int x);
unsigned int *get_racines(unsigned int racine, int n);
unsigned int *eval(Poly P, unsigned int *racines);
unsigned int mod_add(unsigned int a, unsigned int b, unsigned int p);
unsigned int mod_sub(unsigned int a, unsigned int b, unsigned int p);
unsigned int mod_mult(unsigned int a, unsigned int b, unsigned int p);
unsigned int inv(unsigned int a, unsigned int p);

#endif