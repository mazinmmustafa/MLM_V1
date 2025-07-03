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
	complex_t k_z=sqrt_Riemann(k*k-k_rho*k_rho, sheet_I);
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
	}
	spectral.k_z = k_z;
	spectral.Theta = k_z*d;
	spectral.Z_e = k_z/(omega*_eps_0*eps);
	spectral.Z_h = (omega*_mu_0*mu)/k_z;
	return spectral;
}

void Refl(const complex_t k_rho, config_t *config, const size_t n, const sheet_t sheet, 
	complex_t *Gamma_e_u, complex_t *Gamma_h_u, complex_t *Gamma_e_d, complex_t *Gamma_h_d){
	spectral_parameters_t para_m, para_n;
	complex_t Gamma_mn;
	complex_t A, B;
	// Up
	(*Gamma_e_u) = config->Gamma_u;
	(*Gamma_h_u) = config->Gamma_u;
	for (size_t i=1; i<=n; i++){
		para_m = spectral_components(k_rho, config, i-1, sheet);
		para_n = spectral_components(k_rho, config, i, sheet);
		// e-mode
		Gamma_mn = (para_m.Z_e-para_n.Z_e)/(para_m.Z_e+para_n.Z_e);
		A = Gamma_mn+(*Gamma_e_u)*cexp(-_1j*2.0*para_m.Theta);
		B = 1.0+Gamma_mn*(*Gamma_e_u)*cexp(-_1j*2.0*para_m.Theta);
		(*Gamma_e_u) = A/B;
		// h-mode
		Gamma_mn = (para_m.Z_h-para_n.Z_h)/(para_m.Z_h+para_n.Z_h);
		A = Gamma_mn+(*Gamma_h_u)*cexp(-_1j*2.0*para_m.Theta);
		B = 1.0+Gamma_mn*(*Gamma_h_u)*cexp(-_1j*2.0*para_m.Theta);
		(*Gamma_h_u) = A/B;
	}
	// Down
	(*Gamma_e_d) = config->Gamma_d;
	(*Gamma_h_d) = config->Gamma_d;
	for (int_t i=config->N-2; i>=(int_t)n; i--){
		para_m = spectral_components(k_rho, config, (size_t)i+1, sheet);
		para_n = spectral_components(k_rho, config, (size_t)i, sheet);
		// e-mode
		Gamma_mn = (para_m.Z_e-para_n.Z_e)/(para_m.Z_e+para_n.Z_e);
		A = Gamma_mn+(*Gamma_e_d)*cexp(-_1j*2.0*para_m.Theta);
		B = 1.0+Gamma_mn*(*Gamma_e_d)*cexp(-_1j*2.0*para_m.Theta);
		(*Gamma_e_d) = A/B;
		// h-mode
		Gamma_mn = (para_m.Z_h-para_n.Z_h)/(para_m.Z_h+para_n.Z_h);
		A = Gamma_mn+(*Gamma_h_d)*cexp(-_1j*2.0*para_m.Theta);
		B = 1.0+Gamma_mn*(*Gamma_h_d)*cexp(-_1j*2.0*para_m.Theta);
		(*Gamma_h_d) = A/B;
	}
}

void tau(const complex_t k_rho, config_t *config, const size_t n, const sheet_t sheet, 
	const complex_t V_m_e, const complex_t V_m_h, 
	complex_t *V_e_u, complex_t *V_h_u, complex_t *V_e_d, complex_t *V_h_d){
	spectral_parameters_t para_m, para_n;
	complex_t Gamma_m_e, Gamma_m_h, Gamma_n_e, Gamma_n_h;
	complex_t dummy_1, dummy_2;
	complex_t A, B;
	// Up
	(*V_e_u) = config->Gamma_u;
	(*V_h_u) = config->Gamma_u;
	for (int_t i=config->N-2; i>=(int_t)n; i--){
		para_m = spectral_components(k_rho, config, (size_t)i+1, sheet);
		para_n = spectral_components(k_rho, config, (size_t)i, sheet);
		// e-mode
		Refl(k_rho, config, (size_t)i+1, sheet, &Gamma_m_e, &Gamma_m_h, &dummy_1, &dummy_2);
		A = (1.0+Gamma_m_e)*cexp(-_1j*para_m.Theta)/(1.0+cexp(-_1j*2.0*para_m.Theta));
		B = 1.0+Gamma_mn*(*Gamma_e_u)*cexp(-_1j*2.0*para_m.Theta);
		(*Gamma_e_u) = A/B;
		// h-mode
		Gamma_mn = (para_m.Z_h-para_n.Z_h)/(para_m.Z_h+para_n.Z_h);
		A = Gamma_mn+(*Gamma_h_u)*cexp(-_1j*2.0*para_m.Theta);
		B = 1.0+Gamma_mn*(*Gamma_h_u)*cexp(-_1j*2.0*para_m.Theta);
		(*Gamma_h_u) = A/B;
	}
	// Down
	(*V_e_d) = config->Gamma_d;
	(*V_h_d) = config->Gamma_d;
	for (int_t i=config->N-2; i>=(int_t)n; i--){
		para_m = spectral_components(k_rho, config, (size_t)i+1, sheet);
		para_n = spectral_components(k_rho, config, (size_t)i, sheet);
		// e-mode
		Gamma_mn = (para_m.Z_e-para_n.Z_e)/(para_m.Z_e+para_n.Z_e);
		A = Gamma_mn+(*Gamma_e_d)*cexp(-_1j*2.0*para_m.Theta);
		B = 1.0+Gamma_mn*(*Gamma_e_d)*cexp(-_1j*2.0*para_m.Theta);
		(*Gamma_e_d) = A/B;
		// h-mode
		Gamma_mn = (para_m.Z_h-para_n.Z_h)/(para_m.Z_h+para_n.Z_h);
		A = Gamma_mn+(*Gamma_h_d)*cexp(-_1j*2.0*para_m.Theta);
		B = 1.0+Gamma_mn*(*Gamma_h_d)*cexp(-_1j*2.0*para_m.Theta);
		(*Gamma_h_d) = A/B;
	}
}

