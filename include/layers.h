#ifndef __LAYERS_H__
#define __LAYERS_H__

// libraries
#include "basic_definitions.h"
#include "utilities.h"

// definitions
typedef struct layer_t{
    complex_t mu, eps, k, eta, lambda;
    real_t z_u, z_d, d;
}layer_t;

typedef struct config_t{
    size_t N;
    complex_t gamma_u, gamma_d;
    real_t freq, omega;
    layer_t *layers;
    size_t is_allocated;
    real_t k_0, lambda_0;
}config_t;
#define default_config_t {0, 0.0, 0.0, 0.0, 0.0, null, false, 0.0, 0.0}

// Functions
void init_layer_t(layer_t *layer, const complex_t mu, const complex_t eps,
    const real_t z_u, const real_t z_d);
void init_config_t(config_t *config, const size_t N, const complex_t gamma_u, const complex_t gamma_d, 
    const real_t freq);
void clear_config_t(config_t *config);
void set_layer_config_t(config_t *config, const size_t n, const complex_t mu, const complex_t eps,
    const real_t z_d, const real_t z_u);
void check_config_t(config_t *config);
void print_config_t(config_t *config);

#endif