#ifndef __MATH_UTILITIES_HPP__
#define __MATH_UTILITIES_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"

// Definitions

// Functions
real_t deg_to_rad(const real_t x);
real_t rad_to_deg(const real_t x);
//
complex_t sinc(const complex_t z);
real_t sinc(const real_t x);
complex_t sinc_dx(const complex_t z);
real_t sinc_dx(const real_t x);

#endif