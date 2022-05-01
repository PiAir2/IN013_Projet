#ifndef VECT_H
#define VECT_H

#include <immintrin.h>

typedef unsigned int Uint;

void vect_mod_add(Uint *res1, Uint *tab1, Uint *tab2);
void vect_mod_sub(Uint *res1, Uint *tab1, Uint *tab2);
void vect_mod_mult(Uint *res, Uint *tab1, Uint *tab2);
void vect_mod_add_sub_eval(Uint *res_add, Uint *res_sub, Uint *tab1, Uint *tab2);
void vect_mod_mult_eval(Uint *res, Uint *tab1, Uint *tab2, Uint i, Uint pas_rac, __m256i *u);
// __m256i mod_x(Uint *res, __m256i x, Uint i, Uint p);

#endif
