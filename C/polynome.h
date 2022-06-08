#ifndef POLYNOME_H
#define POLYNOME_H

#include <math.h>
#include <immintrin.h>

// #define NB_P 2013265921 // racine = 2, ordre = (P-1)/2, 31 bits = R
#define NB_P 754974721 // racine = 11, ordre = P-1, 30 bits = R
                       // 754974721 = 1 + 2^24 * 3^2 * 5 => degré max du polynôme = 2^24

#define NB_R log(NB_P)/log(2) // 2^(R-1) < P <=  2^R
// choix de T et S tels que T >= R et S+T < n+R-1, ici n=32
#define NB_T 30
#define NB_S 27
#define NB_Q (int) (pow(2,NB_T+NB_S)/NB_P) // Q = floor(2^(S+T)/P)

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
Uint *eval_malloc(Poly P, Uint *racines);
Uint *eval(Uint *coeffs, Uint taille, Uint *tmp_coeffs, Uint *racines, Uint pas_rac);
Uint *vect_eval(Uint *coeffs, Uint taille, Uint *tmp_coeffs, Uint *racines, Uint pas_rac, Uint *tmp_sub);
float *get_racines_inverse(Uint racine, Uint n);
Uint *eval_inv(Uint *coeffs, Uint deg, Uint *tmp_coeffs, float *racines, Uint pas_rac);
Poly FFT(Poly P, Poly Q);

#endif