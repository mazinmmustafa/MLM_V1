#ifndef __SPECIAL_FUNCTIONS_HPP__
#define __SPECIAL_FUNCTIONS_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "amos.hpp"

// Definitions

// Functions
complex_t besselj(const real_t n, const complex_t z);
complex_t bessely(const real_t n, const complex_t z);
complex_t besseli(const real_t n, const complex_t z);
complex_t besselk(const real_t n, const complex_t z);
complex_t besselh1(const real_t n, const complex_t z);
complex_t besselh2(const real_t n, const complex_t z);

// 

complex_t sphj(const real_t n, const complex_t z);
complex_t sphy(const real_t n, const complex_t z);
complex_t sphh1(const real_t n, const complex_t z);
complex_t sphh2(const real_t n, const complex_t z);

#endif