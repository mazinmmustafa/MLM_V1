#ifndef __TESTBENCH_HPP__
#define __TESTBENCH_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "special_functions.hpp"
#include "matrix.hpp"
#include "math_utilities.hpp"
#include "quadl.hpp"
#include "engine.hpp"

// Definitions

// Functions
void test_configuration();
void test_TLGFs();
void test_integrands();
void test_DGFs_Paulus();
void test_DGFs_Chew();
void test_DGFs_Gold_Kretschmann();
void test_DGFs_Gold_Kretschmann_cut();
void test_DGFs_Paulus_near_field();
void test_planar_wave_Chew_cut();
void test_Gold_Kretschmann_reflection();
void test_planar_wave_Gold_Kretschmann_near_field();
void test_plasmonic_WG_far_field();
void test_Chew_far_field();

#endif