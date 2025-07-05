#ifndef __ENGINE_HPP__
#define __ENGINE_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "special_functions.hpp"
#include "math_utilities.hpp"
#include "quadl.hpp"

// Definitions
struct layer_t{
    size_t n=0;
    real_t z_min=0.0, z_max=0.0;
    real_t d=0.0;
    complex_t k=0.0, eta=0.0;
    complex_t mu=1.0, eps=1.0;
};

enum sheet_t{sheet_I, sheet_II, sheet_III, sheet_IV};

struct k_rho_paramters_t{
    complex_t k_z=0.0, Theta=0.0, Z_e=0.0, Z_h=0.0; 
};

struct Gamma_t{
    complex_t Gamma_e=0.0, Gamma_h=0.0;
};

class configuration_t{
    private:
        
    public:
        size_t N=1;
        real_t freq=0.0, omega=0.0, k_0=0.0, lambda_0=0.0; 
        complex_t Gamma_u=0.0, Gamma_d=0.0;
        size_t is_allocated=false;
        layer_t *layers=null;

        configuration_t();
        ~configuration_t();
        void set(const size_t N, const real_t freq);
        void set_boundary_conditions(const complex_t Gamma_u, const complex_t Gamma_d);
        void unset();
        void add_layer(const size_t n, const real_t z_min, const real_t z_max, const complex_t mu, const complex_t eps);
        void log();
        //
        k_rho_paramters_t get_k_rho_parameters(const complex_t k_rho, const size_t n, const sheet_t sheet);
        Gamma_t refl_u(const complex_t k_rho, const size_t n, const sheet_t sheet);
        Gamma_t refl_d(const complex_t k_rho, const size_t n, const sheet_t sheet);
};

// Functions
complex_t sqrt_m(const complex_t z, const sheet_t sheet);

#endif