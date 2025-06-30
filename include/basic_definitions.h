#ifndef __BASIC_DEFINITIONS_H__
#define __BASIC_DEFINITIONS_H__

// librariries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <stdint.h>
#include <limits.h>
#include <unistd.h>

// definitions
#define true 1
#define false 0
#define null NULL

#ifdef _WIN64
#define __windows__
#endif

#define int_t int
#define real_t double
#define complex_t double complex

// Constatns

// math
#define _pi 3.14159265358979 
#define _1j _Complex_I

// physics
#define _c_0 299792458.0 // m/s
#define _mu_0 1.25663706127E-6 // H/m
#define _eps_0 8.8541878188E-12 // F/m
#define _eta_0 376.730313412 // ohm
#define _h_bar 6.582119569E-16 // eV.s
#define _q_e 1.602176634E-19 // C
#define _m_e 9.1093837139E-31 // kg

// Engineering
#define _m 1.0E+0
#define _cm 1.0E-2
#define _mm 1.0E-3
#define _um 1.0E-6
#define _nm 1.0E-9
#define _THz 1.0E+12
#define _GHz 1.0E+9
#define _MHz 1.0E+6
#define _kHz 1.0E+3
#define _Hz 1.0E+0

#endif