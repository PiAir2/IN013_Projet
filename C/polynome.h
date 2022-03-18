#ifndef POLYNOME_H
#define POLYNOME_H

#define NB_P 2013265921

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
long mod_pow(long x, long n);
long horner(Poly P, long x);
long get_racine(long racine, long ordre_racine, int deg);
long *get_racines(long racine, int n);
long *eval(Poly P, long *racines);
long mod_add(long a, long b, long p);
long mod_sub(long a, long b, long p);
long mod_mult(long a, long b, long p);
long inv(long a, long p);

#endif