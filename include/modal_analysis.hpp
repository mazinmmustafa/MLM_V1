#ifndef __MODAL_ANALYSIS_HPP__
#define __MODAL_ANALYSIS_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "special_functions.hpp"
#include "math_utilities.hpp"
#include "quadl.hpp"
#include "vector.hpp"

// Definitions
struct contour_t{
    real_t x_1=0.0, y_1=0.0, x_2=0.0, y_2=0.0;
};

struct function_args_t{
    complex_t (*func)(const complex_t, void*);
    contour_t contour;
    real_t k=0.0;
    void *args=null;
    complex_t dummy=0.0;
};

// Functions
complex_t function_derivative(complex_t func(const complex_t, void*), const complex_t z, void *args);
complex_t compute_mu_k(complex_t func(const complex_t, void*), void *args_, 
    contour_t contour, const real_t k, quadl_t quadl);
void find_polynomial_roots(const size_t N, const real_t C[], complex_t *roots);

#endif