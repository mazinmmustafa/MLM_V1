#ifndef __ENGINE_1d_HPP__
#define __ENGINE_1d_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "special_functions.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "math_utilities.hpp"
#include "quadl.hpp"
#include "shape.hpp"
#include "projection.hpp"

// Definitions
class engine_1d_t{
    private:
        real_t k, lambda;
        matrix_t<complex_t> Z_mn, V_m, I_n;
    public:
        shape_t shape;
        engine_1d_t(){}
        ~engine_1d_t(){}
};

// Functions
real_t I_1_1d(const basis_1d_t bn, const vector_t<real_t> r, const real_t a);
vector_t<real_t> I_2_1d(const basis_1d_t bn, const vector_t<real_t> r, const real_t a);
vector_t<real_t> I_3_1d(const basis_1d_t bn, const vector_t<real_t> r, const real_t a);

complex_t Z_mn_1d(basis_1d_t bm, basis_1d_t bn, const real_t k, 
    const real_t lambda, quadl_domain_t quadl, const real_t a, const real_t eta);

#endif