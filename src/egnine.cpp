//
#include "engine.hpp"

configuration_t::configuration_t(){}
    
configuration_t::~configuration_t(){}

void configuration_t::set(const size_t N, const real_t freq){
    assert_error(this->is_allocated==false, "configuration is already set");
    assert_error(N>0, "invalid number of layers");
    assert_error(freq>0.0, "invalid operating frequency");
    this->N = N;
    this->freq = freq;
    this->omega = 2.0*pi*freq;
    this->lambda_0 = c_0/freq;
    this->k_0 = 2.0*pi/this->lambda_0;
    this->layers = (layer_t*)calloc(N, sizeof(layer_t));
    this->is_allocated = true;
}

void configuration_t::set_boundary_conditions(const complex_t Gamma_u, const complex_t Gamma_d){
    this->Gamma_u = Gamma_u;
    this->Gamma_d = Gamma_d;
}

void configuration_t::unset(){
    assert_error(this->is_allocated==true, "configuration is already set");
    this->N = 0;
    this->freq = 0.0;
    this->omega = 0.0;
    this->k_0 = 0.0;
    this->lambda_0 = 0.0;
    this->Gamma_d = 0.0;
    this->Gamma_u = 0.0;
    free(this->layers);
    this->is_allocated = false;
}

void configuration_t::add_layer(const size_t n, const real_t z_min, const real_t z_max, const complex_t mu, const complex_t eps){
    assert_error(this->is_allocated==true, "configuration is not set yet");
    assert_error(n<this->N, "invalid layer number");
    assert_error(z_max>=z_min, "invalid layer description");
    this->layers[n].z_min = z_min;
    this->layers[n].z_max = z_max;
    this->layers[n].d = z_max-z_min;
    this->layers[n].mu = mu;
    this->layers[n].eps = eps;
    this->layers[n].k = this->k_0*sqrt(mu)*sqrt(eps);
    this->layers[n].eta = eta_0*sqrt(mu)/sqrt(eps);
}

void configuration_t::log(){
    assert_error(this->is_allocated==true, "configuration is not set yet");
    file_t file;
    file.open("data/configuration.txt", 'w');
    file.write("================================================================================\n");
    file.write("number of layers\t\t\t: %2zu\n", this->N);
    file.write("operating frequency\t\t\t: %21.14E\n", this->freq);
    file.write("operating radial frequency\t: %21.14E\n", this->omega);
    file.write("free-space wavenumber\t\t: %21.14E\n", this->k_0);
    file.write("free-space wavelength\t\t: %21.14E\n", this->lambda_0);
    file.write("boundary Gamma up\t\t\t: %21.14E, %21.14E\n", real(this->Gamma_u), imag(this->Gamma_u));
    file.write("boundary Gamma down\t\t\t: %21.14E, %21.14E\n", real(this->Gamma_d), imag(this->Gamma_d));
    for (size_t n=0; n<this->N; n++){
        file.write("================================================================================\n");
        file.write("layer %zu:\n", n);
        file.write("z_max\t\t\t\t: %21.14E\n", this->layers[n].z_max);
        file.write("z_min\t\t\t\t: %21.14E\n", this->layers[n].z_min);
        file.write("layer thickness\t\t: %21.14E\n", this->layers[n].d);
        file.write("mu\t\t\t\t\t: %21.14E, %21.14E\n", real(this->layers[n].mu), imag(this->layers[n].mu));
        file.write("eps\t\t\t\t\t: %21.14E, %21.14E\n", real(this->layers[n].eps), imag(this->layers[n].eps));
        file.write("k\t\t\t\t\t: %21.14E, %21.14E\n", real(this->layers[n].k), imag(this->layers[n].k));
        file.write("eta\t\t\t\t\t: %21.14E, %21.14E\n", real(this->layers[n].eta), imag(this->layers[n].eta));
    }
    file.write("================================================================================\n");
    file.close();
}

//

complex_t sqrt_m(const complex_t z, const sheet_t sheet){
    const complex_t w=sqrt(z);
    if (sheet==sheet_I){
        return imag(w)<=0.0 ? +w : -w;
    }
    if (sheet==sheet_II){
        return imag(w)>=0.0 ? +w : -w;
    }
    return 0.0;
}

k_rho_paramters_t configuration_t::get_k_rho_parameters(const complex_t k_rho, const size_t n, const sheet_t sheet){
    k_rho_paramters_t para;
    const complex_t k=this->layers[n].k;
    para.k_z = sqrt_m(k*k-k_rho*k_rho, sheet_I);
    if (sheet==sheet_I){
        if (n==0){
            para.k_z = sqrt_m(k*k-k_rho*k_rho, sheet_I);
        }
        if (n==this->N-1){
            para.k_z = sqrt_m(k*k-k_rho*k_rho, sheet_I);
        }
    }
    if (sheet==sheet_II){
        if (n==0){
            para.k_z = sqrt_m(k*k-k_rho*k_rho, sheet_II);
        }
        if (n==this->N-1){
            para.k_z = sqrt_m(k*k-k_rho*k_rho, sheet_I);
        }
    }
    if (sheet==sheet_III){
        if (n==0){
            para.k_z = sqrt_m(k*k-k_rho*k_rho, sheet_I);
        }
        if (n==this->N-1){
            para.k_z = sqrt_m(k*k-k_rho*k_rho, sheet_II);
        }
    }
    if (sheet==sheet_IV){
        if (n==0){
            para.k_z = sqrt_m(k*k-k_rho*k_rho, sheet_II);
        }
        if (n==this->N-1){
            para.k_z = sqrt_m(k*k-k_rho*k_rho, sheet_II);
        }
    }
    para.Theta = para.k_z*this->layers[n].d;
    const complex_t mu=this->layers[n].mu;
    const complex_t eps=this->layers[n].eps;
    const real_t omega=this->omega;
    para.Z_e = para.k_z/(omega*eps_0*eps);
    para.Z_h = (omega*mu_0*mu)/para.k_z;
    return para;
}

Gamma_t configuration_t::refl_u(const complex_t k_rho, const size_t n, const sheet_t sheet){
    complex_t Gamma_e=this->Gamma_u;
    complex_t Gamma_h=this->Gamma_u;
    Gamma_t Gamma;
    const complex_t j=complex_t(0.0, 1.0);
    for (size_t i=1; i<=n; i++){
        k_rho_paramters_t para_m=configuration_t::get_k_rho_parameters(k_rho, i-1, sheet);
        k_rho_paramters_t para_n=configuration_t::get_k_rho_parameters(k_rho, i, sheet);
        const complex_t Z_e_m=para_m.Z_e;
        const complex_t Z_h_m=para_m.Z_h;
        const complex_t Z_e_n=para_n.Z_e;
        const complex_t Z_h_n=para_n.Z_h;
        const complex_t Gamma_mn_e=(Z_e_m-Z_e_n)/(Z_e_m+Z_e_n);
        const complex_t Gamma_mn_h=(Z_h_m-Z_h_n)/(Z_h_m+Z_h_n);
        const complex_t Theta_m=para_m.Theta;
        complex_t A, B;
        A = Gamma_mn_e+Gamma_e*exp(-j*2.0*Theta_m);
        B = 1.0+Gamma_mn_e*Gamma_e*exp(-j*2.0*Theta_m);
        Gamma_e = A/B;
        A = Gamma_mn_h+Gamma_h*exp(-j*2.0*Theta_m);
        B = 1.0+Gamma_mn_h*Gamma_h*exp(-j*2.0*Theta_m);
        Gamma_h = A/B;
    }
    Gamma.Gamma_e = Gamma_e;
    Gamma.Gamma_h = Gamma_h;
    return Gamma;
}

Gamma_t configuration_t::refl_d(const complex_t k_rho, const size_t n, const sheet_t sheet){
    complex_t Gamma_e=this->Gamma_d;
    complex_t Gamma_h=this->Gamma_d;
    Gamma_t Gamma;
    const complex_t j=complex_t(0.0, 1.0);
    for (int_t i=(int_t)this->N-2; i>=(int_t)n; i--){
        k_rho_paramters_t para_m=configuration_t::get_k_rho_parameters(k_rho, (size_t)i+1, sheet);
        k_rho_paramters_t para_n=configuration_t::get_k_rho_parameters(k_rho, (size_t)i, sheet);
        const complex_t Z_e_m=para_m.Z_e;
        const complex_t Z_h_m=para_m.Z_h;
        const complex_t Z_e_n=para_n.Z_e;
        const complex_t Z_h_n=para_n.Z_h;
        const complex_t Gamma_mn_e=(Z_e_m-Z_e_n)/(Z_e_m+Z_e_n);
        const complex_t Gamma_mn_h=(Z_h_m-Z_h_n)/(Z_h_m+Z_h_n);
        const complex_t Theta_m=para_m.Theta;
        complex_t A, B;
        A = Gamma_mn_e+Gamma_e*exp(-j*2.0*Theta_m);
        B = 1.0+Gamma_mn_e*Gamma_e*exp(-j*2.0*Theta_m);
        Gamma_e = A/B;
        A = Gamma_mn_h+Gamma_h*exp(-j*2.0*Theta_m);
        B = 1.0+Gamma_mn_h*Gamma_h*exp(-j*2.0*Theta_m);
        Gamma_h = A/B;
    }
    Gamma.Gamma_e = Gamma_e;
    Gamma.Gamma_h = Gamma_h;
    return Gamma;
}