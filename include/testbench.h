#ifndef __TESTBENCH_H__
#define __TESTBENCH_H__

// libraries
#include "basic_definitions.h"
#include "utilities.h"
#include "layers.h"
#include "bessel.h"
#include "quad.h"
#include "transmission_line.h"
//
#include "testbench_gold_Kretschmann.h"

// definitions

// Functions
void test_basic_definitions();
void test_utilities();
void test_layers();
void test_bessel();
void test_quad();
void test_transmission_line();

#endif