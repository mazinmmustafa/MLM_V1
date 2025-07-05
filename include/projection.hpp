#ifndef __PROJECTION_HPP__
#define __PROJECTION_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "math_utilities.hpp"
#include "vector.hpp"
#include "mesh.hpp"
#include "shape.hpp"

// Definitions
struct projection_para_1d_t{
    real_t P_m, P_p, l_m, l_p, P0;
    vector_t<real_t> l_unit, P0_unit;
};

// Functions
projection_para_1d_t projection_1d(const vector_t<real_t> v1, const vector_t<real_t> v2, 
    const vector_t<real_t> p, const real_t a);

#endif