#ifndef __ENGINE_HPP__
#define __ENGINE_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "special_functions.hpp"
#include "math_utilities.hpp"
#include "quadl.hpp"
#include "vector.hpp"

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

struct tau_t{
    complex_t tau_e=0.0, tau_h=0.0;
};

struct TLGF_t{
    complex_t V_e=0.0, I_e=0.0;
    complex_t V_h=0.0, I_h=0.0;
};

enum source_t{v_source, i_source};

struct Greens_functions_t{
    complex_t xx=0.0, xy=0.0, xz=0.0;
    complex_t yx=0.0, yy=0.0, yz=0.0;
    complex_t zx=0.0, zy=0.0, zz=0.0;
};

struct position_t{
    real_t x=0.0, y=0.0, z=0.0, rho=0.0, phi=0.0, r=0.0, theta=0.0;
};

struct cartesian_t:position_t{
    cartesian_t(const real_t x, const real_t y, const real_t z){
        this->x = x;
        this->y = y;
        this->z = z;
        //
        this->r = sqrt(x*x+y*y+z*z);
        this->rho = sqrt(x*x+y*y);
        this->phi = atan2(y, x);
        this->theta = acos(z/r);
    }
};   

struct cylindrical_t:position_t{
    cylindrical_t(const real_t rho, const real_t phi, const real_t z){
        this->rho = rho;
        this->phi = phi;
        this->z = z;
        //
        this->x = rho*cos(phi);
        this->y = rho*sin(phi);
        this->r = sqrt(rho*rho+z*z);
        this->theta = acos(z/r);
    }
};   

struct spherical_t:position_t{
    spherical_t(const real_t r, const real_t theta, const real_t phi){
        this->r = r;
        this->theta = theta;
        this->phi = phi;
        //
        this->x = r*sin(theta)*cos(phi);
        this->y = r*sin(theta)*sin(phi);
        this->r = r*cos(theta);
        this->rho = r*sin(theta);
    }
};  

struct integrand_args{
    position_t r, r_;
};

struct dipole_t{
    position_t r;
    real_t theta=0.0, phi=0.0;
    complex_t Il=0.0;
    complex_t x=0.0, y=0.0, z=0.0;
    dipole_t(const position_t r, const real_t theta, const real_t phi, const complex_t Il){
        this->r = r;
        this->theta = theta;
        this->phi = phi;
        this->Il = Il;
        this->x = Il*sin(theta)*cos(phi);
        this->y = Il*sin(theta)*sin(phi);
        this->z = Il*cos(theta);
    }
};

struct near_field_t{
    complex_t x=0.0, y=0.0, z=0.0;
};

struct incident_plane_wave_field_E_0_t{
    complex_t TM=0.0, TE=0.0;
};

struct near_field_plane_wave_t{
    near_field_t E, H;
};

struct far_field_components_t{
    complex_t theta=0.0, phi=0.0;
};

struct far_field_t{
    far_field_components_t E, H;
};

enum termination_t{PEC, PMC, INF};

enum polarization_t{e, h};


class configuration_t{
    private:
        
    public:
        size_t N=1;
        real_t freq=0.0, omega=0.0, k_0=0.0, lambda_0=0.0; 
        complex_t Gamma_u=0.0, Gamma_d=0.0;
        size_t is_allocated=false;
        layer_t *layers=null;
        real_t k_min=+1.0;
        //
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
        tau_t tau_u(const complex_t k_rho, const complex_t V_m_e, const complex_t V_n_e,
            const size_t m, const size_t n, const sheet_t sheet);
        tau_t tau_d(const complex_t k_rho, const complex_t V_m_e, const complex_t V_n_e,
            const size_t m, const size_t n, const sheet_t sheet);
        size_t find_layer(const real_t z);
        TLGF_t TLGF_r(const complex_t k_rho, const real_t z, const real_t z_, const source_t source, const sheet_t sheet);
        TLGF_t TLGF(const complex_t k_rho, const real_t z, const real_t z_, const source_t source, const sheet_t sheet);
        complex_t detour(complex_t &k_rho, const real_t rho, const real_t distance);
        complex_t hankel_transform(complex_t func(complex_t, void*), void *args, quadl_t quadl);
        // EJ
        Greens_functions_t G_EJ_0(const position_t r, const position_t r_);
        complex_t G_EJ_xx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EJ_xy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EJ_xz(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EJ_yx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EJ_yy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EJ_yz(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EJ_zx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EJ_zy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EJ_zz(const position_t r, const position_t r_, quadl_t quadl);
        // EM
        Greens_functions_t G_EM_0(const position_t r, const position_t r_);
        complex_t G_EM_xx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EM_xy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EM_xz(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EM_yx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EM_yy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EM_yz(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EM_zx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EM_zy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_EM_zz(const position_t r, const position_t r_, quadl_t quadl);
        // E-field
        near_field_t compute_E_J_near_field(const position_t r, const dipole_t dipole_J, quadl_t quadl);
        near_field_t compute_E_M_near_field(const position_t r, const dipole_t dipole_M, quadl_t quadl);
        // HJ
        Greens_functions_t G_HJ_0(const position_t r, const position_t r_);
        complex_t G_HJ_xx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HJ_xy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HJ_xz(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HJ_yx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HJ_yy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HJ_yz(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HJ_zx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HJ_zy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HJ_zz(const position_t r, const position_t r_, quadl_t quadl);
        // HM
        Greens_functions_t G_HM_0(const position_t r, const position_t r_);
        complex_t G_HM_xx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HM_xy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HM_xz(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HM_yx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HM_yy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HM_yz(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HM_zx(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HM_zy(const position_t r, const position_t r_, quadl_t quadl);
        complex_t G_HM_zz(const position_t r, const position_t r_, quadl_t quadl);
        // H-field
        near_field_t compute_H_J_near_field(const position_t r, const dipole_t dipole_J, quadl_t quadl);
        near_field_t compute_H_M_near_field(const position_t r, const dipole_t dipole_M, quadl_t quadl);
        // plane wave
        near_field_plane_wave_t compute_plane_wave(const position_t r, const real_t theta_i, const real_t phi_i, 
            const incident_plane_wave_field_E_0_t E_0);
        Gamma_t compute_Gamma(const real_t theta_i, const sheet_t sheet);
        // far field
        far_field_t compute_far_field_J(const real_t r, const real_t theta_s, const real_t phi_s, 
            const dipole_t dipole_J);
        far_field_t compute_far_field_M(const real_t r, const real_t theta_s, const real_t phi_s, 
            const dipole_t dipole_M);
        // modal analysis
};

struct modal_args_t{
    termination_t term_top, term_bottom;
    polarization_t pol;
    sheet_t sheet;
    configuration_t config;
};

// Functions
complex_t sqrt_m(const complex_t z, const sheet_t sheet);


complex_t dispersion_function_on_sheet(const complex_t k_rho, void *args_);
complex_t dispersion_function_all_sheets(const complex_t k_rho, void *args_);

#endif