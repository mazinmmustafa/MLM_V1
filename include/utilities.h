#ifndef __UTILITIES_H__
#define __UTILITIES_H__

// libraries
#include "basic_definitions.h"

// definitions
#define assert_error(condition, msg) __assert_error((condition), (msg), __FILE__, __LINE__);

// Functions
void __assert_error(const size_t condition, const char *msg, const char *file, const size_t line);
void progress_bar(const size_t n, const size_t Ns, const char *msg);
void printc(const complex_t z);

real_t linspace(const real_t x_min, const real_t x_max, const size_t Ns, const size_t i);
real_t logspace(const real_t x_min, const real_t x_max, const size_t Ns, const size_t i);

#endif