#ifndef __BASIC_DEFINITIONS_HPP__
#define __BASIC_DEFINITIONS_HPP__

// librariries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <limits.h>

// definitions
#define true 1
#define false 0
#define null NULL

#ifdef _WIN64
#define __windows__
#endif

#define int_t int32_t
#define real_t double
#define complex_t std::complex<real_t>

#undef I

// constants

// math
const real_t pi=3.14159265358979;

// physics
const real_t c_0=299792458.0; // m/s
const real_t mu_0=1.25663706127E-6; // H/m
const real_t eps_0=8.8541878188E-12; // F/m
const real_t eta_0=sqrt(mu_0/eps_0); // ohm
const real_t h_bar=6.582119569E-16; // eV.s
const real_t q_e=1.602176634E-19; // C
const real_t m_e=9.1093837139E-31; // kg

// engineering
namespace units{
    const real_t m=1.0E+0;
    const real_t cm=1.0E-2;
    const real_t mm=1.0E-3;
    const real_t um=1.0E-6;
    const real_t nm=1.0E-9;
    const real_t THz=1.0E+12;
    const real_t GHz=1.0E+9;
    const real_t MHz=1.0E+6;
    const real_t kHz=1.0E+3;
    const real_t Hz=1.0E+0;
}

#endif