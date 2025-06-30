#ifndef __BESSEL_H__
#define __BESSEL_H__

// Definitions
#include "basic_definitions.h"

// Functions
complex_t besselj(real_t n, complex_t z);
complex_t bessely(real_t n, complex_t z);
complex_t besseli(real_t n, complex_t z);
complex_t besselk(real_t n, complex_t z);
complex_t besselh(int_t k, real_t n, complex_t z);

#endif