#ifndef __QUAD_H__
#define __QUAD_H__

// Definitions
#include "basic_definitions.h"

// Functions
void get_quad_rule(const size_t N, real_t *x, real_t *w);
complex_t quad_rule(complex_t func(complex_t, void*), const real_t a, const real_t b, void *args,
	real_t *x_quad, real_t *w_quad, const size_t N);
complex_t quadL(complex_t func(complex_t, void*), const real_t a, const real_t b, void *args,
	real_t *x_quad, real_t *w_quad, const size_t N, const real_t tol, const size_t k_max, 
	size_t *flag, const complex_t I_p, size_t k);

#endif