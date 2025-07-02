#ifndef __TRANSMISSION_LINE_H__
#define __TRANSMISSION_LINE_H__

// Definitions
#include "basic_definitions.h"
typedef enum sheet_t{sheet_I, sheet_II, sheet_III, sheet_IV}sheet_t;

// Functions
complex_t sqrt_Riemann(const complex_t z, const sheet_t sheet);

#endif