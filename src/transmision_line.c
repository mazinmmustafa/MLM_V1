//
#include "transmission_line.h"

complex_t sqrt_Riemann(const complex_t z, const sheet_t sheet){
	complex_t w=csqrt(z);
	if (sheet==sheet_I){
		return cimag(w)<=0.0 ? +w : -w;
	}
	if (sheet==sheet_II){
		return cimag(w)<=0.0 ? -w : +w;
	}
	return 0.0;
}

spectral_parameters_t spectral_components(const complex_t k_rho, config_t *config, 
	const size_t n, const sheet_t sheet){
	assert_error(n<config->N, "invalid layer number");
	spectral_parameters_t spectral;
	const complex_t k=config->layers[n].k;
	const real_t d=config->layers[n].d;
	const size_t N=config->N;
	const real_t omega=config->omega;
	const complex_t mu=config->layers[n].mu; 
	const complex_t eps=config->layers[n].eps;
	complex_t k_z=0.0;
	if (sheet==sheet_I){
		if (n==0){
			k_z = sqrt_Riemann(k*k-k_rho*k_rho, sheet_I);
		}
		if (n==(N-1)){
			k_z = sqrt_Riemann(k*k-k_rho*k_rho, sheet_I);
		}
	}else
	if (sheet==sheet_II){
		if (n==0){
			k_z = sqrt_Riemann(k*k-k_rho*k_rho, sheet_II);
		}
		if (n==(N-1)){
			k_z = sqrt_Riemann(k*k-k_rho*k_rho, sheet_I);
		}
	}else
	if (sheet==sheet_III){
		if (n==0){
			k_z = sqrt_Riemann(k*k-k_rho*k_rho, sheet_I);
		}
		if (n==(N-1)){
			k_z = sqrt_Riemann(k*k-k_rho*k_rho, sheet_II);
		}
	}else
	if (sheet==sheet_IV){
		if (n==0){
			k_z = sqrt_Riemann(k*k-k_rho*k_rho, sheet_II);
		}
		if (n==(N-1)){
			k_z = sqrt_Riemann(k*k-k_rho*k_rho, sheet_II);
		}
	}else{
		k_z = sqrt_Riemann(k*k-k_rho*k_rho, sheet_I);
	} 
	spectral.k_z = k_z;
	spectral.Theta = k_z*d;
	spectral.Z_e = k_z/(omega*_eps_0*eps);
	spectral.Z_h = (omega*_mu_0*mu)/k_z;
	return spectral;
}

TL_parameters_t basic_reflection_TL(const complex_t k_rho, config_t *config, const size_t i, const size_t j, 
	const sheet_t sheet){
	TL_parameters_t parameters; 
	spectral_parameters_t para_i=spectral_components(k_rho, config, i, sheet);
	spectral_parameters_t para_j=spectral_components(k_rho, config, j, sheet);
	parameters.Gamma_ij_e = (para_i.Z_e-para_j.Z_e)/(para_i.Z_e+para_j.Z_e);
	parameters.Gamma_ij_h = (para_i.Z_h-para_j.Z_h)/(para_i.Z_h+para_j.Z_h);
	parameters.tau_ij_e = 1.0+parameters.Gamma_ij_e;
	parameters.tau_ij_h = 1.0+parameters.Gamma_ij_h;
	return parameters;
}
