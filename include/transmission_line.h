#ifndef __TRANSMISSION_LINE_H__
#define __TRANSMISSION_LINE_H__

// Definitions
#include "basic_definitions.h"
#include "layers.h"

typedef enum sheet_t{sheet_I, sheet_II, sheet_III, sheet_IV}sheet_t;
typedef struct spectral_parameters_t{
    complex_t k_z, Theta, Z_e, Z_h;
}spectral_parameters_t;
typedef struct TL_parameters_t{
    complex_t Gamma_ij_e, Gamma_ij_h;
    complex_t tau_ij_e, tau_ij_h;
    complex_t Omega_ij_e, Omega_ij_h;
}TL_parameters_t;

// Functions
complex_t sqrt_Riemann(const complex_t z, const sheet_t sheet);
spectral_parameters_t spectral_components(const complex_t k_rho, config_t *config, 
	const size_t n, const sheet_t sheet);

#endif