#ifndef POLYNOME_H
#define POLYNOME_H

//#define NB_P // tests
//#define NB_P 2013265921 // racine = 2, ordre = (NB_P-1)/2
#define NB_P 754974721 // racine = 11, ordre = NB_P-1

typedef unsigned int Uint;

typedef struct _Poly {
    Uint *coeffs;
    Uint deg;
} Poly;

void afficher_poly(Poly P);
Poly gen_poly(Uint deg);
void liberer_poly(Poly P);
int compare_poly(Poly P, Poly Q);
Poly prod_poly_naif(Poly P, Poly Q);
Poly copy_poly(Poly P, Uint deb, Uint fin);
Poly poly_somme(Poly P, Poly Q, Uint deb);
Poly oppose_poly(Poly P);
Poly prod_poly_karatsuba(Poly P, Poly Q);
Uint mod_pow(Uint x, Uint n);
Uint horner(Poly P, Uint x);
Uint *get_racines(Uint racine, Uint n);
Uint mod_add(Uint a, Uint b, Uint p);
Uint mod_sub(Uint a, Uint b, Uint p);
Uint mod_mult(Uint a, Uint b, Uint p);
Uint inv(Uint a, Uint p);
//Uint *eval(Poly P, Uint *racines);
Uint *eval(Uint *coeffs, Uint taille, Uint *tmp_coeffs, Uint *racines, Uint pas_rac);
Uint *vect_eval(Uint *coeffs, Uint taille, Uint *tmp_coeffs, Uint *racines, Uint pas_rac, Uint *tmp_sub);

#endif