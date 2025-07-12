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
    for (size_t i=0; i<this->N; i++){
        if (abs(real(this->layers[i].k))>this->k_min*this->k_0){
            this->k_min = abs(real(this->layers[i].k))/this->k_0;
        }
    }
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


tau_t configuration_t::tau_u(const complex_t k_rho, const complex_t V_m_e, const complex_t V_n_e,
    const size_t m, const size_t n, const sheet_t sheet){
    tau_t tau;
    tau.tau_e = V_m_e;
    tau.tau_h = V_n_e;
    Gamma_t Gamma;
    const complex_t j=complex_t(0.0, 1.0);
    for (int_t i=(int_t)m-2; i>=(int_t)n; i--){
        k_rho_paramters_t para_m=configuration_t::get_k_rho_parameters(k_rho, (size_t)i+1, sheet);
        k_rho_paramters_t para_n=configuration_t::get_k_rho_parameters(k_rho, (size_t)i, sheet);
        const complex_t Theta_m=para_m.Theta;
        const complex_t Theta_n=para_n.Theta;
        Gamma_t Gamma;
        Gamma = configuration_t::refl_u(k_rho, (size_t)i+1, sheet);
        const complex_t Gamma_m_e=Gamma.Gamma_e;
        const complex_t Gamma_m_h=Gamma.Gamma_h;
        Gamma = configuration_t::refl_u(k_rho, (size_t)i, sheet);
        const complex_t Gamma_n_e=Gamma.Gamma_e;
        const complex_t Gamma_n_h=Gamma.Gamma_h;
        complex_t A, B;
        A = (1.0+Gamma_m_e)*exp(-j*Theta_m);
        B = 1.0+Gamma_n_e*exp(-j*2.0*Theta_n);
        tau.tau_e = tau.tau_e*A/B;
        A = (1.0+Gamma_m_h)*exp(-j*Theta_m);
        B = 1.0+Gamma_n_h*exp(-j*2.0*Theta_n);
        tau.tau_h = tau.tau_h*A/B;
    }
    return tau;
}


tau_t configuration_t::tau_d(const complex_t k_rho, const complex_t V_m_e, const complex_t V_n_e,
    const size_t m, const size_t n, const sheet_t sheet){
    tau_t tau;
    tau.tau_e = V_m_e;
    tau.tau_h = V_n_e;
    Gamma_t Gamma;
    const complex_t j=complex_t(0.0, 1.0);
    for (size_t i=m+2; i<=n; i++){
        k_rho_paramters_t para_m=configuration_t::get_k_rho_parameters(k_rho, i-1, sheet);
        k_rho_paramters_t para_n=configuration_t::get_k_rho_parameters(k_rho, i, sheet);
        const complex_t Theta_m=para_m.Theta;
        const complex_t Theta_n=para_n.Theta;
        Gamma_t Gamma;
        Gamma = configuration_t::refl_d(k_rho, i-1, sheet);
        const complex_t Gamma_m_e=Gamma.Gamma_e;
        const complex_t Gamma_m_h=Gamma.Gamma_h;
        Gamma = configuration_t::refl_d(k_rho, i, sheet);
        const complex_t Gamma_n_e=Gamma.Gamma_e;
        const complex_t Gamma_n_h=Gamma.Gamma_h;
        complex_t A, B;
        A = (1.0+Gamma_m_e)*exp(-j*Theta_m);
        B = 1.0+Gamma_n_e*exp(-j*2.0*Theta_n);
        tau.tau_e = tau.tau_e*A/B;
        A = (1.0+Gamma_m_h)*exp(-j*Theta_m);
        B = 1.0+Gamma_n_h*exp(-j*2.0*Theta_n);
        tau.tau_h = tau.tau_h*A/B;
    }
    return tau;
}

size_t configuration_t::find_layer(const real_t z){
    for (size_t i=0; i<this->N; i++){
        if ((z>=this->layers[i].z_min)&&(z<=this->layers[i].z_max)){
            return i;
        }
    }
    assert_error(false, "invalid layer");
    return 0;
}

TLGF_t configuration_t::TLGF_r(const complex_t k_rho, const real_t z, const real_t z_, const source_t source, const sheet_t sheet){
    const size_t m=configuration_t::find_layer(z_);
    const size_t n=configuration_t::find_layer(z);
    //
    complex_t v=0.0, i=0.0;
    if (source==v_source){
        v = 1.0;
    }
    if (source==i_source){
        i = 1.0;
    }
    TLGF_t TLGF;
    //
    complex_t Gamma_e_d_, Gamma_h_d_, Gamma_e_u_, Gamma_h_u_;
    Gamma_e_d_ = Gamma_h_d_ = Gamma_e_u_ = Gamma_h_u_ = 0.0;
    const complex_t j=complex_t(0.0, 1.0);
    real_t z_d=this->layers[m].z_min;
    real_t z_u=this->layers[m].z_max;
    //
    k_rho_paramters_t para_m=configuration_t::get_k_rho_parameters(k_rho, m, sheet);
    complex_t k_z=para_m.k_z;
    complex_t Theta_m=para_m.Theta;
    Gamma_t Gamma_u, Gamma_d;
    Gamma_u = configuration_t::refl_u(k_rho, m, sheet);
    Gamma_d = configuration_t::refl_d(k_rho, m, sheet);
    complex_t D_e=1.0-Gamma_d.Gamma_e*Gamma_u.Gamma_e*exp(-j*2.0*Theta_m); 
    complex_t D_h=1.0-Gamma_d.Gamma_h*Gamma_u.Gamma_h*exp(-j*2.0*Theta_m);
    complex_t A, B;
    A = (1.0-Gamma_d.Gamma_e*exp(-j*2.0*k_z*(z_-z_d)))*v;
    B = (1.0+Gamma_d.Gamma_e*exp(-j*2.0*k_z*(z_-z_d)))*i*para_m.Z_e;
    complex_t V_e_p=(+A+B)/(2.0*D_e);
    A = (1.0-Gamma_d.Gamma_h*exp(-j*2.0*k_z*(z_-z_d)))*v;
    B = (1.0+Gamma_d.Gamma_h*exp(-j*2.0*k_z*(z_-z_d)))*i*para_m.Z_h;
    complex_t V_h_p=(+A+B)/(2.0*D_h);
    A = (1.0-Gamma_u.Gamma_e*exp(+j*2.0*k_z*(z_-z_u)))*v;
    B = (1.0+Gamma_u.Gamma_e*exp(+j*2.0*k_z*(z_-z_u)))*i*para_m.Z_e;
    complex_t V_e_n=(-A+B)/(2.0*D_e);
    A = (1.0-Gamma_u.Gamma_h*exp(+j*2.0*k_z*(z_-z_u)))*v;
    B = (1.0+Gamma_u.Gamma_h*exp(+j*2.0*k_z*(z_-z_u)))*i*para_m.Z_h;
    complex_t V_h_n=(-A+B)/(2.0*D_h);
    Gamma_e_u_ = Gamma_u.Gamma_e*exp(+j*2.0*k_z*(z_-z_u));
    Gamma_h_u_ = Gamma_u.Gamma_h*exp(+j*2.0*k_z*(z_-z_u));
    Gamma_e_d_ = Gamma_d.Gamma_e*exp(-j*2.0*k_z*(z_-z_d));
    Gamma_h_d_ = Gamma_d.Gamma_h*exp(-j*2.0*k_z*(z_-z_d));
    //
    if (n==m){
        k_rho_paramters_t para_m=configuration_t::get_k_rho_parameters(k_rho, m, sheet);
        complex_t k_z=para_m.k_z;
        complex_t Theta_m=para_m.Theta;
        complex_t Gamma_e_d, Gamma_h_d, Gamma_e_u, Gamma_h_u;
        Gamma_e_d = Gamma_h_d = Gamma_e_u = Gamma_h_u = 0.0;
        complex_t D_e=1.0-Gamma_e_d*Gamma_e_u*exp(-j*2.0*Theta_m); 
        complex_t D_h=1.0-Gamma_h_d*Gamma_h_u*exp(-j*2.0*Theta_m);
        complex_t A, B;
        A = (1.0-Gamma_e_d*exp(-j*2.0*k_z*(z_-z_d)))*v;
        B = (1.0+Gamma_e_d*exp(-j*2.0*k_z*(z_-z_d)))*i*para_m.Z_e;
        complex_t V_e_p_0=(+A+B)/(2.0*D_e);
        A = (1.0-Gamma_h_d*exp(-j*2.0*k_z*(z_-z_d)))*v;
        B = (1.0+Gamma_h_d*exp(-j*2.0*k_z*(z_-z_d)))*i*para_m.Z_h;
        complex_t V_h_p_0=(+A+B)/(2.0*D_h);
        A = (1.0-Gamma_e_u*exp(+j*2.0*k_z*(z_-z_u)))*v;
        B = (1.0+Gamma_e_u*exp(+j*2.0*k_z*(z_-z_u)))*i*para_m.Z_e;
        complex_t V_e_n_0=(-A+B)/(2.0*D_e);
        A = (1.0-Gamma_h_u*exp(+j*2.0*k_z*(z_-z_u)))*v;
        B = (1.0+Gamma_h_u*exp(+j*2.0*k_z*(z_-z_u)))*i*para_m.Z_h;
        complex_t V_h_n_0=(-A+B)/(2.0*D_h);
        if (z>=z_){
            TLGF.V_e = +((V_e_p-V_e_p_0)*exp(-j*k_z*(z-z_))+V_e_p*Gamma_e_u_*exp(+j*k_z*(z-z_)));
            TLGF.I_e = +((V_e_p-V_e_p_0)*exp(-j*k_z*(z-z_))-V_e_p*Gamma_e_u_*exp(+j*k_z*(z-z_)))/para_m.Z_e;
            TLGF.V_h = +((V_h_p-V_h_p_0)*exp(-j*k_z*(z-z_))+V_h_p*Gamma_h_u_*exp(+j*k_z*(z-z_)));
            TLGF.I_h = +((V_h_p-V_h_p_0)*exp(-j*k_z*(z-z_))-V_h_p*Gamma_h_u_*exp(+j*k_z*(z-z_)))/para_m.Z_h;
        }else{
            TLGF.V_e = +((V_e_n-V_e_n_0)*exp(+j*k_z*(z-z_))+V_e_n*Gamma_e_d_*exp(-j*k_z*(z-z_)));
            TLGF.I_e = -((V_e_n-V_e_n_0)*exp(+j*k_z*(z-z_))-V_e_n*Gamma_e_d_*exp(-j*k_z*(z-z_)))/para_m.Z_e;
            TLGF.V_h = +((V_h_n-V_h_n_0)*exp(+j*k_z*(z-z_))+V_h_n*Gamma_h_d_*exp(-j*k_z*(z-z_)));
            TLGF.I_h = -((V_h_n-V_h_n_0)*exp(+j*k_z*(z-z_))-V_h_n*Gamma_h_d_*exp(-j*k_z*(z-z_)))/para_m.Z_h;
        }
    }
    //
    if (n<m){
        k_rho_paramters_t para_n=configuration_t::get_k_rho_parameters(k_rho, m-1, sheet);
        complex_t Theta_n=para_n.Theta;
        Gamma_t Gamma_u;
        Gamma_u = configuration_t::refl_u(k_rho, m-1, sheet);
        real_t z_u=this->layers[m].z_max;
        V_e_p = V_e_p*(exp(-j*k_z*(z_u-z_))+Gamma_e_u_*exp(+j*k_z*(z_u-z_)))/(1.0+Gamma_u.Gamma_e*exp(-j*2.0*Theta_n));
        V_h_p = V_h_p*(exp(-j*k_z*(z_u-z_))+Gamma_h_u_*exp(+j*k_z*(z_u-z_)))/(1.0+Gamma_u.Gamma_h*exp(-j*2.0*Theta_n));
        para_n = configuration_t::get_k_rho_parameters(k_rho, n, sheet);
        k_z = para_n.k_z;
        Theta_n = para_n.Theta;
        z_u = this->layers[n].z_max;
        Gamma_u = configuration_t::refl_u(k_rho, n, sheet);
        tau_t tau;
        tau = tau_u(k_rho, V_e_p, V_h_p, m, n, sheet);
        TLGF.V_e = tau.tau_e*exp(-j*Theta_n)*(exp(-j*k_z*(z-z_u))+Gamma_u.Gamma_e*exp(+j*k_z*(z-z_u)));
        TLGF.I_e = tau.tau_e*exp(-j*Theta_n)*(exp(-j*k_z*(z-z_u))-Gamma_u.Gamma_e*exp(+j*k_z*(z-z_u)))/para_n.Z_e;
        TLGF.V_h = tau.tau_h*exp(-j*Theta_n)*(exp(-j*k_z*(z-z_u))+Gamma_u.Gamma_h*exp(+j*k_z*(z-z_u)));
        TLGF.I_h = tau.tau_h*exp(-j*Theta_n)*(exp(-j*k_z*(z-z_u))-Gamma_u.Gamma_h*exp(+j*k_z*(z-z_u)))/para_n.Z_h;
    }
    //
    if (n>m){
        k_rho_paramters_t para_n=configuration_t::get_k_rho_parameters(k_rho, m+1, sheet);
        complex_t Theta_n=para_n.Theta;
        Gamma_t Gamma_d;
        Gamma_d = configuration_t::refl_d(k_rho, m+1, sheet);
        real_t z_d=this->layers[m].z_min;
        V_e_n = V_e_n*(exp(+j*k_z*(z_d-z_))+Gamma_e_d_*exp(-j*k_z*(z_d-z_)))/(1.0+Gamma_d.Gamma_e*exp(-j*2.0*Theta_n));
        V_h_n = V_h_n*(exp(+j*k_z*(z_d-z_))+Gamma_h_d_*exp(-j*k_z*(z_d-z_)))/(1.0+Gamma_d.Gamma_h*exp(-j*2.0*Theta_n));
        para_n = configuration_t::get_k_rho_parameters(k_rho, n, sheet);
        k_z = para_n.k_z;
        Theta_n = para_n.Theta;
        z_d = this->layers[n].z_min;
        Gamma_d = configuration_t::refl_d(k_rho, n, sheet);
        tau_t tau;
        tau = tau_d(k_rho, V_e_n, V_h_n, m, n, sheet);
        TLGF.V_e = +tau.tau_e*exp(-j*Theta_n)*(exp(+j*k_z*(z-z_d))+Gamma_d.Gamma_e*exp(-j*k_z*(z-z_d)));
        TLGF.I_e = -tau.tau_e*exp(-j*Theta_n)*(exp(+j*k_z*(z-z_d))-Gamma_d.Gamma_e*exp(-j*k_z*(z-z_d)))/para_n.Z_e;
        TLGF.V_h = +tau.tau_h*exp(-j*Theta_n)*(exp(+j*k_z*(z-z_d))+Gamma_d.Gamma_h*exp(-j*k_z*(z-z_d)));
        TLGF.I_h = -tau.tau_h*exp(-j*Theta_n)*(exp(+j*k_z*(z-z_d))-Gamma_d.Gamma_h*exp(-j*k_z*(z-z_d)))/para_n.Z_h;
    }
    return TLGF;
}

TLGF_t configuration_t::TLGF(const complex_t k_rho, const real_t z, const real_t z_, const source_t source, const sheet_t sheet){
    const size_t m=configuration_t::find_layer(z_);
    const size_t n=configuration_t::find_layer(z);
    //
    complex_t v=0.0, i=0.0;
    if (source==v_source){
        v = 1.0;
    }
    if (source==i_source){
        i = 1.0;
    }
    TLGF_t TLGF;
    //
    complex_t Gamma_e_d_, Gamma_h_d_, Gamma_e_u_, Gamma_h_u_;
    Gamma_e_d_ = Gamma_h_d_ = Gamma_e_u_ = Gamma_h_u_ = 0.0;
    const complex_t j=complex_t(0.0, 1.0);
    real_t z_d=this->layers[m].z_min;
    real_t z_u=this->layers[m].z_max;
    //
    k_rho_paramters_t para_m=configuration_t::get_k_rho_parameters(k_rho, m, sheet);
    complex_t k_z=para_m.k_z;
    complex_t Theta_m=para_m.Theta;
    Gamma_t Gamma_u, Gamma_d;
    Gamma_u = configuration_t::refl_u(k_rho, m, sheet);
    Gamma_d = configuration_t::refl_d(k_rho, m, sheet);
    complex_t D_e=1.0-Gamma_d.Gamma_e*Gamma_u.Gamma_e*exp(-j*2.0*Theta_m); 
    complex_t D_h=1.0-Gamma_d.Gamma_h*Gamma_u.Gamma_h*exp(-j*2.0*Theta_m);
    complex_t A, B;
    A = (1.0-Gamma_d.Gamma_e*exp(-j*2.0*k_z*(z_-z_d)))*v;
    B = (1.0+Gamma_d.Gamma_e*exp(-j*2.0*k_z*(z_-z_d)))*i*para_m.Z_e;
    complex_t V_e_p=(+A+B)/(2.0*D_e);
    A = (1.0-Gamma_d.Gamma_h*exp(-j*2.0*k_z*(z_-z_d)))*v;
    B = (1.0+Gamma_d.Gamma_h*exp(-j*2.0*k_z*(z_-z_d)))*i*para_m.Z_h;
    complex_t V_h_p=(+A+B)/(2.0*D_h);
    A = (1.0-Gamma_u.Gamma_e*exp(+j*2.0*k_z*(z_-z_u)))*v;
    B = (1.0+Gamma_u.Gamma_e*exp(+j*2.0*k_z*(z_-z_u)))*i*para_m.Z_e;
    complex_t V_e_n=(-A+B)/(2.0*D_e);
    A = (1.0-Gamma_u.Gamma_h*exp(+j*2.0*k_z*(z_-z_u)))*v;
    B = (1.0+Gamma_u.Gamma_h*exp(+j*2.0*k_z*(z_-z_u)))*i*para_m.Z_h;
    complex_t V_h_n=(-A+B)/(2.0*D_h);
    Gamma_e_u_ = Gamma_u.Gamma_e*exp(+j*2.0*k_z*(z_-z_u));
    Gamma_h_u_ = Gamma_u.Gamma_h*exp(+j*2.0*k_z*(z_-z_u));
    Gamma_e_d_ = Gamma_d.Gamma_e*exp(-j*2.0*k_z*(z_-z_d));
    Gamma_h_d_ = Gamma_d.Gamma_h*exp(-j*2.0*k_z*(z_-z_d));
    //
    if (n==m){
        k_rho_paramters_t para_m=configuration_t::get_k_rho_parameters(k_rho, m, sheet);
        complex_t k_z=para_m.k_z;
        complex_t Gamma_e_d, Gamma_h_d, Gamma_e_u, Gamma_h_u;
        Gamma_e_d = Gamma_h_d = Gamma_e_u = Gamma_h_u = 0.0;
        if (z>=z_){
            TLGF.V_e = +(V_e_p*exp(-j*k_z*(z-z_))+V_e_p*Gamma_e_u_*exp(+j*k_z*(z-z_)));
            TLGF.I_e = +(V_e_p*exp(-j*k_z*(z-z_))-V_e_p*Gamma_e_u_*exp(+j*k_z*(z-z_)))/para_m.Z_e;
            TLGF.V_h = +(V_h_p*exp(-j*k_z*(z-z_))+V_h_p*Gamma_h_u_*exp(+j*k_z*(z-z_)));
            TLGF.I_h = +(V_h_p*exp(-j*k_z*(z-z_))-V_h_p*Gamma_h_u_*exp(+j*k_z*(z-z_)))/para_m.Z_h;
        }else{
            TLGF.V_e = +(V_e_n*exp(+j*k_z*(z-z_))+V_e_n*Gamma_e_d_*exp(-j*k_z*(z-z_)));
            TLGF.I_e = -(V_e_n*exp(+j*k_z*(z-z_))-V_e_n*Gamma_e_d_*exp(-j*k_z*(z-z_)))/para_m.Z_e;
            TLGF.V_h = +(V_h_n*exp(+j*k_z*(z-z_))+V_h_n*Gamma_h_d_*exp(-j*k_z*(z-z_)));
            TLGF.I_h = -(V_h_n*exp(+j*k_z*(z-z_))-V_h_n*Gamma_h_d_*exp(-j*k_z*(z-z_)))/para_m.Z_h;
        }
    }
    //
    if (n<m){
        k_rho_paramters_t para_n=configuration_t::get_k_rho_parameters(k_rho, m-1, sheet);
        complex_t Theta_n=para_n.Theta;
        Gamma_t Gamma_u;
        Gamma_u = configuration_t::refl_u(k_rho, m-1, sheet);
        real_t z_u=this->layers[m].z_max;
        V_e_p = V_e_p*(exp(-j*k_z*(z_u-z_))+Gamma_e_u_*exp(+j*k_z*(z_u-z_)))/(1.0+Gamma_u.Gamma_e*exp(-j*2.0*Theta_n));
        V_h_p = V_h_p*(exp(-j*k_z*(z_u-z_))+Gamma_h_u_*exp(+j*k_z*(z_u-z_)))/(1.0+Gamma_u.Gamma_h*exp(-j*2.0*Theta_n));
        para_n = configuration_t::get_k_rho_parameters(k_rho, n, sheet);
        k_z = para_n.k_z;
        Theta_n = para_n.Theta;
        z_u = this->layers[n].z_max;
        Gamma_u = configuration_t::refl_u(k_rho, n, sheet);
        tau_t tau;
        tau = tau_u(k_rho, V_e_p, V_h_p, m, n, sheet);
        TLGF.V_e = tau.tau_e*exp(-j*Theta_n)*(exp(-j*k_z*(z-z_u))+Gamma_u.Gamma_e*exp(+j*k_z*(z-z_u)));
        TLGF.I_e = tau.tau_e*exp(-j*Theta_n)*(exp(-j*k_z*(z-z_u))-Gamma_u.Gamma_e*exp(+j*k_z*(z-z_u)))/para_n.Z_e;
        TLGF.V_h = tau.tau_h*exp(-j*Theta_n)*(exp(-j*k_z*(z-z_u))+Gamma_u.Gamma_h*exp(+j*k_z*(z-z_u)));
        TLGF.I_h = tau.tau_h*exp(-j*Theta_n)*(exp(-j*k_z*(z-z_u))-Gamma_u.Gamma_h*exp(+j*k_z*(z-z_u)))/para_n.Z_h;
    }
    //
    if (n>m){
        k_rho_paramters_t para_n=configuration_t::get_k_rho_parameters(k_rho, m+1, sheet);
        complex_t Theta_n=para_n.Theta;
        Gamma_t Gamma_d;
        Gamma_d = configuration_t::refl_d(k_rho, m+1, sheet);
        real_t z_d=this->layers[m].z_min;
        V_e_n = V_e_n*(exp(+j*k_z*(z_d-z_))+Gamma_e_d_*exp(-j*k_z*(z_d-z_)))/(1.0+Gamma_d.Gamma_e*exp(-j*2.0*Theta_n));
        V_h_n = V_h_n*(exp(+j*k_z*(z_d-z_))+Gamma_h_d_*exp(-j*k_z*(z_d-z_)))/(1.0+Gamma_d.Gamma_h*exp(-j*2.0*Theta_n));
        para_n = configuration_t::get_k_rho_parameters(k_rho, n, sheet);
        k_z = para_n.k_z;
        Theta_n = para_n.Theta;
        z_d = this->layers[n].z_min;
        Gamma_d = configuration_t::refl_d(k_rho, n, sheet);
        tau_t tau;
        tau = tau_d(k_rho, V_e_n, V_h_n, m, n, sheet);
        TLGF.V_e = +tau.tau_e*exp(-j*Theta_n)*(exp(+j*k_z*(z-z_d))+Gamma_d.Gamma_e*exp(-j*k_z*(z-z_d)));
        TLGF.I_e = -tau.tau_e*exp(-j*Theta_n)*(exp(+j*k_z*(z-z_d))-Gamma_d.Gamma_e*exp(-j*k_z*(z-z_d)))/para_n.Z_e;
        TLGF.V_h = +tau.tau_h*exp(-j*Theta_n)*(exp(+j*k_z*(z-z_d))+Gamma_d.Gamma_h*exp(-j*k_z*(z-z_d)));
        TLGF.I_h = -tau.tau_h*exp(-j*Theta_n)*(exp(+j*k_z*(z-z_d))-Gamma_d.Gamma_h*exp(-j*k_z*(z-z_d)))/para_n.Z_h;
    }
    return TLGF;
}

complex_t configuration_t::detour(complex_t &k_rho, const real_t rho, const real_t distance){
    const real_t k_0=this->k_0;
    k_rho = k_0*k_rho;
    const complex_t j=complex_t(0.0, 1.0);
    const real_t a=this->k_0*(1.0+this->k_min);
    // For testing
    // const real_t b=rho>distance ? (this->k_0<1.0/rho ? this->k_0 : 1.0/rho) : this->k_0;
    // const real_t b=this->k_0*1.0E-3 + (0.0*rho*distance);
    // const real_ ratio=1.0+(0.0*distance);
    //
    const real_t b=(0.0*distance) + this->k_0<1.0/rho ? this->k_0 : 1.0/rho;
    const real_t t=real(k_rho);
    const real_t x=t;
    const real_t y=t<a ? b*sin(pi*t/a) : 0.0;
    k_rho = x+j*y;
    return 1.0+j*(pi*b/a)*cos(pi*t/a);
}

Greens_functions_t configuration_t::G_EJ_0(const position_t r, const position_t r_){
    Greens_functions_t DGF;
    const size_t n=configuration_t::find_layer(r_.z);
    const complex_t k=this->layers[n].k;
    const complex_t mu=this->layers[n].mu;
    const real_t omega=this->omega;
    const real_t R=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y)+(r.z-r_.z)*(r.z-r_.z));
    const complex_t j=complex_t(0.0, 1.0);
    const complex_t g=exp(-j*k*R)/(4.0*pi*R);
    const complex_t g1=-(1.0+j*k*R)*g/R;
    const complex_t g2=-(1.0+j*k*R)*g1/R+g/(R*R);
    const complex_t factor=-j*omega*mu_0*mu;
    //
    if (R==0.0){
        DGF.xx = 0.0*factor*(-1.0/(3.0*k*k));
        DGF.yy = 0.0*factor*(-1.0/(3.0*k*k));
        DGF.zz = 0.0*factor*(-1.0/(3.0*k*k));
        DGF.xy = DGF.xz = DGF.yx = DGF.yz = DGF.zx = DGF.zy = 0.0;
    }else{
        DGF.xx = factor*(g+(1.0/(k*k))*(g1/R+((r.x-r_.x)/R)*((r.x-r_.x)/R)*(g2-g1/R)));
        DGF.xy = factor*((1.0/(k*k))*(((r.x-r_.x)/R)*((r.y-r_.y)/R)*(g2-g1/R)));
        DGF.xz = factor*((1.0/(k*k))*(((r.x-r_.x)/R)*((r.z-r_.z)/R)*(g2-g1/R)));
        DGF.yx = factor*((1.0/(k*k))*(((r.y-r_.y)/R)*((r.x-r_.x)/R)*(g2-g1/R)));
        DGF.yy = factor*(g+(1.0/(k*k))*(g1/R+((r.y-r_.y)/R)*((r.y-r_.y)/R)*(g2-g1/R)));
        DGF.yz = factor*((1.0/(k*k))*(((r.y-r_.y)/R)*((r.z-r_.z)/R)*(g2-g1/R)));
        DGF.zx = factor*((1.0/(k*k))*(((r.z-r_.z)/R)*((r.x-r_.x)/R)*(g2-g1/R)));
        DGF.zy = factor*((1.0/(k*k))*(((r.z-r_.z)/R)*(((r.y-r_.y)-r_.y)/R)*(g2-g1/R)));
        DGF.zz = factor*(g+(1.0/(k*k))*(g1/R+((r.z-r_.z)/R)*((r.z-r_.z)/R)*(g2-g1/R)));
    }
    return DGF;
}

Greens_functions_t configuration_t::G_HM_0(const position_t r, const position_t r_){
    Greens_functions_t DGF;
    const size_t n=configuration_t::find_layer(r_.z);
    const complex_t k=this->layers[n].k;
    const complex_t eps=this->layers[n].eps;
    const real_t omega=this->omega;
    const real_t R=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y)+(r.z-r_.z)*(r.z-r_.z));
    const complex_t j=complex_t(0.0, 1.0);
    const complex_t g=exp(-j*k*R)/(4.0*pi*R);
    const complex_t g1=-(1.0+j*k*R)*g/R;
    const complex_t g2=-(1.0+j*k*R)*g1/R+g/(R*R);
    const complex_t factor=-j*omega*eps_0*eps;
    //
    if (R==0.0){
        DGF.xx = 0.0*factor*(-1.0/(3.0*k*k));
        DGF.yy = 0.0*factor*(-1.0/(3.0*k*k));
        DGF.zz = 0.0*factor*(-1.0/(3.0*k*k));
        DGF.xy = DGF.xz = DGF.yx = DGF.yz = DGF.zx = DGF.zy = 0.0;
    }else{
        DGF.xx = factor*(g+(1.0/(k*k))*(g1/R+((r.x-r_.x)/R)*((r.x-r_.x)/R)*(g2-g1/R)));
        DGF.xy = factor*((1.0/(k*k))*(((r.x-r_.x)/R)*((r.y-r_.y)/R)*(g2-g1/R)));
        DGF.xz = factor*((1.0/(k*k))*(((r.x-r_.x)/R)*((r.z-r_.z)/R)*(g2-g1/R)));
        DGF.yx = factor*((1.0/(k*k))*(((r.y-r_.y)/R)*((r.x-r_.x)/R)*(g2-g1/R)));
        DGF.yy = factor*(g+(1.0/(k*k))*(g1/R+((r.y-r_.y)/R)*((r.y-r_.y)/R)*(g2-g1/R)));
        DGF.yz = factor*((1.0/(k*k))*(((r.y-r_.y)/R)*((r.z-r_.z)/R)*(g2-g1/R)));
        DGF.zx = factor*((1.0/(k*k))*(((r.z-r_.z)/R)*((r.x-r_.x)/R)*(g2-g1/R)));
        DGF.zy = factor*((1.0/(k*k))*(((r.z-r_.z)/R)*(((r.y-r_.y)-r_.y)/R)*(g2-g1/R)));
        DGF.zz = factor*(g+(1.0/(k*k))*(g1/R+((r.z-r_.z)/R)*((r.z-r_.z)/R)*(g2-g1/R)));
    }
    return DGF;
}

Greens_functions_t configuration_t::G_HJ_0(const position_t r, const position_t r_){
    Greens_functions_t DGF;
    const size_t n=configuration_t::find_layer(r_.z);
    const complex_t k=this->layers[n].k;
    const real_t R=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y)+(r.z-r_.z)*(r.z-r_.z));
    const complex_t j=complex_t(0.0, 1.0);
    const complex_t g=exp(-j*k*R)/(4.0*pi*R);
    const complex_t g1=-(1.0+j*k*R)*g/R;
    //
    if (R==0.0){
        DGF.xx = 0.0;
        DGF.yy = 0.0;
        DGF.zz = 0.0;
        DGF.xy = DGF.xz = DGF.yx = DGF.yz = DGF.zx = DGF.zy = 0.0;
    }else{
        DGF.xx = 0.0;
        DGF.xy = -((r.z-r_.z)/R)*g1;
        DGF.xz = +((r.y-r_.y)/R)*g1;
        DGF.yx = +((r.z-r_.z)/R)*g1;
        DGF.yy = 0.0;
        DGF.yz = -((r.x-r_.x)/R)*g1;
        DGF.zx = -((r.y-r_.y)/R)*g1;
        DGF.zy = +((r.x-r_.x)/R)*g1;
        DGF.zz = 0.0;
    }
    return DGF;
}

Greens_functions_t configuration_t::G_EM_0(const position_t r, const position_t r_){
    Greens_functions_t DGF;
    const size_t n=configuration_t::find_layer(r_.z);
    const complex_t k=this->layers[n].k;
    const real_t R=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y)+(r.z-r_.z)*(r.z-r_.z));
    const complex_t j=complex_t(0.0, 1.0);
    const complex_t g=exp(-j*k*R)/(4.0*pi*R);
    const complex_t g1=-(1.0+j*k*R)*g/R;
    //
    if (R==0.0){
        DGF.xx = 0.0;
        DGF.yy = 0.0;
        DGF.zz = 0.0;
        DGF.xy = DGF.xz = DGF.yx = DGF.yz = DGF.zx = DGF.zy = 0.0;
    }else{
        DGF.xx = 0.0;
        DGF.xy = +((r.z-r_.z)/R)*g1;
        DGF.xz = -((r.y-r_.y)/R)*g1;
        DGF.yx = -((r.z-r_.z)/R)*g1;
        DGF.yy = 0.0;
        DGF.yz = +((r.x-r_.x)/R)*g1;
        DGF.zx = +((r.y-r_.y)/R)*g1;
        DGF.zy = -((r.x-r_.x)/R)*g1;
        DGF.zz = 0.0;
    }
    return DGF;
}

//

struct integrand_args_t{
    real_t rho=0.0, phi=0.0, z=0.0, z_=0.0;
    configuration_t config;
};

complex_t configuration_t::hankel_transform(complex_t func(complex_t, void*), void *args, quadl_t quadl){
    size_t flag;
    complex_t ans=quadl.integral_1d(func, args, 0.0, k_min, flag);
    if (flag==true){print("no convergence!\n"); exit(1);}
    return ans/lambda_0;
}

// EJ

complex_t G_EJ_1(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, i_source, sheet_I);
    return factor*(TLGF.V_e+TLGF.V_h)*besselj(0.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_EJ_2(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, i_source, sheet_I);
    return factor*(TLGF.V_e-TLGF.V_h)*besselj(2.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_EJ_3(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, v_source, sheet_I);
    return factor*(k_rho_*TLGF.V_e)*besselj(1.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_EJ_4(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, i_source, sheet_I);
    return factor*(k_rho_*TLGF.I_e)*besselj(1.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_EJ_5(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, v_source, sheet_I);
    return factor*(k_rho_*k_rho_*TLGF.I_e)*besselj(0.0, k_rho_*args->rho)*k_rho_;
}

complex_t configuration_t::G_EJ_xx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EJ_1_=configuration_t::hankel_transform(G_EJ_1, &args, quadl);
    complex_t G_EJ_2_=configuration_t::hankel_transform(G_EJ_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-0.5*G_EJ_1_;
    complex_t I_2=+0.5*cos(2.0*phi)*G_EJ_2_;
    complex_t I_=I_1+I_2;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EJ_0(r, r_);
        I_ = I_+DGF.xx;
    }
    return I_;
}

complex_t configuration_t::G_EJ_xy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EJ_2_=configuration_t::hankel_transform(G_EJ_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=+0.5*sin(2.0*phi)*G_EJ_2_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EJ_0(r, r_);
        I_ = I_+DGF.xy;
    }
    return I_;
}

complex_t configuration_t::G_EJ_yx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EJ_2_=configuration_t::hankel_transform(G_EJ_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=+0.5*sin(2.0*phi)*G_EJ_2_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EJ_0(r, r_);
        I_ = I_+DGF.yx;
    }
    return I_;
}

complex_t configuration_t::G_EJ_yy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EJ_1_=configuration_t::hankel_transform(G_EJ_1, &args, quadl);
    complex_t G_EJ_2_=configuration_t::hankel_transform(G_EJ_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-0.5*G_EJ_1_;
    complex_t I_2=-0.5*cos(2.0*phi)*G_EJ_2_;
    complex_t I_=I_1+I_2;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EJ_0(r, r_);
        I_ = I_+DGF.yy;
    }
    return I_;
}

complex_t configuration_t::G_EJ_xz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EJ_3_=configuration_t::hankel_transform(G_EJ_3, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((eta_0)/(j*this->k_0*this->layers[m].eps))*cos(phi)*G_EJ_3_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EJ_0(r, r_);
        I_ = I_+DGF.xz;
    }
    return I_;
}

complex_t configuration_t::G_EJ_yz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EJ_3_=configuration_t::hankel_transform(G_EJ_3, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((eta_0)/(j*this->k_0*this->layers[m].eps))*sin(phi)*G_EJ_3_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EJ_0(r, r_);
        I_ = I_+DGF.yz;
    }
    return I_;
}

complex_t configuration_t::G_EJ_zx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EJ_4_=configuration_t::hankel_transform(G_EJ_4, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((eta_0)/(j*this->k_0*this->layers[n].eps))*cos(phi)*G_EJ_4_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EJ_0(r, r_);
        I_ = I_+DGF.zx;
    }
    return I_;
}

complex_t configuration_t::G_EJ_zy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EJ_4_=configuration_t::hankel_transform(G_EJ_4, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((eta_0)/(j*this->k_0*this->layers[n].eps))*sin(phi)*G_EJ_4_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EJ_0(r, r_);
        I_ = I_+DGF.zy;
    }
    return I_;
}

complex_t configuration_t::G_EJ_zz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EJ_5_=configuration_t::hankel_transform(G_EJ_5, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-((eta_0*eta_0)/(this->k_0*this->k_0*this->layers[m].eps*this->layers[n].eps))*G_EJ_5_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EJ_0(r, r_);
        I_ = I_+DGF.zz;
    }
    return I_;
}

near_field_t configuration_t::compute_E_J_near_field(const position_t r, const dipole_t dipole_J, 
    quadl_t quadl){
    near_field_t E;
    E.x = this->G_EJ_xx(r, dipole_J.r, quadl)*dipole_J.x+this->G_EJ_xy(r, dipole_J.r, quadl)*dipole_J.y+this->G_EJ_xz(r, dipole_J.r, quadl)*dipole_J.z;
    E.y = this->G_EJ_yx(r, dipole_J.r, quadl)*dipole_J.x+this->G_EJ_yy(r, dipole_J.r, quadl)*dipole_J.y+this->G_EJ_yz(r, dipole_J.r, quadl)*dipole_J.z;
    E.z = this->G_EJ_zx(r, dipole_J.r, quadl)*dipole_J.x+this->G_EJ_zy(r, dipole_J.r, quadl)*dipole_J.y+this->G_EJ_zz(r, dipole_J.r, quadl)*dipole_J.z;
    return E;
}

// EM

complex_t G_EM_1(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, v_source, sheet_I);
    return factor*(TLGF.V_e-TLGF.V_h)*besselj(2.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_EM_2(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, v_source, sheet_I);
    return factor*(TLGF.V_e+TLGF.V_h)*besselj(0.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_EM_3(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, i_source, sheet_I);
    return factor*(k_rho_*TLGF.V_h)*besselj(1.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_EM_4(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, v_source, sheet_I);
    return factor*(k_rho_*TLGF.I_e)*besselj(1.0, k_rho_*args->rho)*k_rho_;
}

complex_t configuration_t::G_EM_xx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EM_1_=configuration_t::hankel_transform(G_EM_1, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-0.5*sin(2.0*phi)*G_EM_1_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EM_0(r, r_);
        I_ = I_+DGF.xx;
    }
    return I_;
}

complex_t configuration_t::G_EM_xy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EM_1_=configuration_t::hankel_transform(G_EM_1, &args, quadl);
    complex_t G_EM_2_=configuration_t::hankel_transform(G_EM_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-0.5*G_EM_2_;
    complex_t I_2=+0.5*cos(2.0*phi)*G_EM_1_;
    complex_t I_=I_1+I_2;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EM_0(r, r_);
        I_ = I_+DGF.xy;
    }
    return I_;
}

complex_t configuration_t::G_EM_yx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EM_1_=configuration_t::hankel_transform(G_EM_1, &args, quadl);
    complex_t G_EM_2_=configuration_t::hankel_transform(G_EM_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=+0.5*G_EM_2_;
    complex_t I_2=+0.5*cos(2.0*phi)*G_EM_1_;
    complex_t I_=I_1+I_2;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EM_0(r, r_);
        I_ = I_+DGF.yx;
    }
    return I_;
}

complex_t configuration_t::G_EM_yy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EM_1_=configuration_t::hankel_transform(G_EM_1, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=+0.5*sin(2.0*phi)*G_EM_1_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EM_0(r, r_);
        I_ = I_+DGF.yy;
    }
    return I_;
}

complex_t configuration_t::G_EM_xz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EM_3_=configuration_t::hankel_transform(G_EM_3, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((1.0)/(j*eta_0*this->k_0*this->layers[m].mu))*sin(phi)*G_EM_3_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EM_0(r, r_);
        I_ = I_+DGF.xz;
    }
    return I_;
}

complex_t configuration_t::G_EM_yz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EM_3_=configuration_t::hankel_transform(G_EM_3, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-((1.0)/(j*eta_0*this->k_0*this->layers[m].mu))*cos(phi)*G_EM_3_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EM_0(r, r_);
        I_ = I_+DGF.yz;
    }
    return I_;
}

complex_t configuration_t::G_EM_zx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EM_4_=configuration_t::hankel_transform(G_EM_4, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((j*eta_0)/(this->k_0*this->layers[n].eps))*sin(phi)*G_EM_4_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EM_0(r, r_);
        I_ = I_+DGF.zx;
    }
    return I_;
}

complex_t configuration_t::G_EM_zy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_EM_4_=configuration_t::hankel_transform(G_EM_4, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-((j*eta_0)/(this->k_0*this->layers[n].eps))*cos(phi)*G_EM_4_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_EM_0(r, r_);
        I_ = I_+DGF.zy;
    }
    return I_;
}

complex_t configuration_t::G_EM_zz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    return 0.0*rho*phi*quadl.N;
}

near_field_t configuration_t::compute_E_M_near_field(const position_t r, const dipole_t dipole_M, 
    quadl_t quadl){
    near_field_t E;
    E.x = this->G_EM_xx(r, dipole_M.r, quadl)*dipole_M.x+this->G_EM_xy(r, dipole_M.r, quadl)*dipole_M.y+this->G_EM_xz(r, dipole_M.r, quadl)*dipole_M.z;
    E.y = this->G_EM_yx(r, dipole_M.r, quadl)*dipole_M.x+this->G_EM_yy(r, dipole_M.r, quadl)*dipole_M.y+this->G_EM_yz(r, dipole_M.r, quadl)*dipole_M.z;
    E.z = this->G_EM_zx(r, dipole_M.r, quadl)*dipole_M.x+this->G_EM_zy(r, dipole_M.r, quadl)*dipole_M.y+this->G_EM_zz(r, dipole_M.r, quadl)*dipole_M.z;
    return E;
}

// HM

complex_t G_HM_1(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, v_source, sheet_I);
    return factor*(TLGF.I_h+TLGF.I_e)*besselj(0.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_HM_2(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, v_source, sheet_I);
    return factor*(TLGF.I_h-TLGF.I_e)*besselj(2.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_HM_3(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, i_source, sheet_I);
    return factor*(k_rho_*TLGF.I_h)*besselj(1.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_HM_4(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, v_source, sheet_I);
    return factor*(k_rho_*TLGF.V_h)*besselj(1.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_HM_5(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, i_source, sheet_I);
    return factor*(k_rho_*k_rho_*TLGF.V_h)*besselj(0.0, k_rho_*args->rho)*k_rho_;
}

complex_t configuration_t::G_HM_xx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HM_1_=configuration_t::hankel_transform(G_HM_1, &args, quadl);
    complex_t G_HM_2_=configuration_t::hankel_transform(G_HM_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-0.5*G_HM_1_;
    complex_t I_2=+0.5*cos(2.0*phi)*G_HM_2_;
    complex_t I_=I_1+I_2;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HM_0(r, r_);
        I_ = I_+DGF.xx;
    }
    return I_;
}

complex_t configuration_t::G_HM_xy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HM_2_=configuration_t::hankel_transform(G_HM_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=+0.5*sin(2.0*phi)*G_HM_2_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HM_0(r, r_);
        I_ = I_+DGF.xy;
    }
    return I_;
}

complex_t configuration_t::G_HM_yx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HM_2_=configuration_t::hankel_transform(G_HM_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=+0.5*sin(2.0*phi)*G_HM_2_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HM_0(r, r_);
        I_ = I_+DGF.yx;
    }
    return I_;
}

complex_t configuration_t::G_HM_yy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HM_1_=configuration_t::hankel_transform(G_HM_1, &args, quadl);
    complex_t G_HM_2_=configuration_t::hankel_transform(G_HM_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-0.5*G_HM_1_;
    complex_t I_2=-0.5*cos(2.0*phi)*G_HM_2_;
    complex_t I_=I_1+I_2;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HM_0(r, r_);
        I_ = I_+DGF.yy;
    }
    return I_;
}

complex_t configuration_t::G_HM_xz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HM_3_=configuration_t::hankel_transform(G_HM_3, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((1.0)/(j*eta_0*this->k_0*this->layers[m].mu))*cos(phi)*G_HM_3_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HM_0(r, r_);
        I_ = I_+DGF.xz;
    }
    return I_;
}

complex_t configuration_t::G_HM_yz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HM_3_=configuration_t::hankel_transform(G_HM_3, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((1.0)/(j*eta_0*this->k_0*this->layers[m].mu))*sin(phi)*G_HM_3_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HM_0(r, r_);
        I_ = I_+DGF.yz;
    }
    return I_;
}

complex_t configuration_t::G_HM_zx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HM_4_=configuration_t::hankel_transform(G_HM_4, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((1.0)/(j*eta_0*this->k_0*this->layers[n].mu))*cos(phi)*G_HM_4_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HM_0(r, r_);
        I_ = I_+DGF.zx;
    }
    return I_;
}

complex_t configuration_t::G_HM_zy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HM_4_=configuration_t::hankel_transform(G_HM_4, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((1.0)/(j*eta_0*this->k_0*this->layers[n].mu))*sin(phi)*G_HM_4_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HM_0(r, r_);
        I_ = I_+DGF.zy;
    }
    return I_;
}

complex_t configuration_t::G_HM_zz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HM_5_=configuration_t::hankel_transform(G_HM_5, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-((1.0)/(eta_0*eta_0*this->k_0*this->k_0*this->layers[m].mu*this->layers[n].mu))*G_HM_5_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HM_0(r, r_);
        I_ = I_+DGF.zz;
    }
    return I_;
}

near_field_t configuration_t::compute_H_M_near_field(const position_t r, const dipole_t dipole_M, 
    quadl_t quadl){
    near_field_t H;
    H.x = this->G_HM_xx(r, dipole_M.r, quadl)*dipole_M.x+this->G_HM_xy(r, dipole_M.r, quadl)*dipole_M.y+this->G_HM_xz(r, dipole_M.r, quadl)*dipole_M.z;
    H.y = this->G_HM_yx(r, dipole_M.r, quadl)*dipole_M.x+this->G_HM_yy(r, dipole_M.r, quadl)*dipole_M.y+this->G_HM_yz(r, dipole_M.r, quadl)*dipole_M.z;
    H.z = this->G_HM_zx(r, dipole_M.r, quadl)*dipole_M.x+this->G_HM_zy(r, dipole_M.r, quadl)*dipole_M.y+this->G_HM_zz(r, dipole_M.r, quadl)*dipole_M.z;
    return H;
}

// HJ

complex_t G_HJ_1(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, i_source, sheet_I);
    return factor*(TLGF.I_h-TLGF.I_e)*besselj(2.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_HJ_2(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, i_source, sheet_I);
    return factor*(TLGF.I_h+TLGF.I_e)*besselj(0.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_HJ_3(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, v_source, sheet_I);
    return factor*(k_rho_*TLGF.I_e)*besselj(1.0, k_rho_*args->rho)*k_rho_;
}

complex_t G_HJ_4(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    const real_t rho=args->rho;
    const real_t z=args->z;
    const real_t z_=args->z_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, rho, abs(z-z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, z, z_, i_source, sheet_I);
    return factor*(k_rho_*TLGF.V_h)*besselj(1.0, k_rho_*args->rho)*k_rho_;
}

complex_t configuration_t::G_HJ_xx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HJ_1_=configuration_t::hankel_transform(G_HJ_1, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=+0.5*sin(2.0*phi)*G_HJ_1_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HJ_0(r, r_);
        I_ = I_+DGF.xx;
    }
    return I_;
}

complex_t configuration_t::G_HJ_xy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HJ_1_=configuration_t::hankel_transform(G_HJ_1, &args, quadl);
    complex_t G_HJ_2_=configuration_t::hankel_transform(G_HJ_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=+0.5*G_HJ_2_;
    complex_t I_2=-0.5*cos(2.0*phi)*G_HJ_1_;
    complex_t I_=I_1+I_2;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HJ_0(r, r_);
        I_ = I_+DGF.xy;
    }
    return I_;
}

complex_t configuration_t::G_HJ_yx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HJ_1_=configuration_t::hankel_transform(G_HJ_1, &args, quadl);
    complex_t G_HJ_2_=configuration_t::hankel_transform(G_HJ_2, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-0.5*G_HJ_2_;
    complex_t I_2=-0.5*cos(2.0*phi)*G_HJ_1_;
    complex_t I_=I_1+I_2;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HJ_0(r, r_);
        I_ = I_+DGF.yx;
    }
    return I_;
}

complex_t configuration_t::G_HJ_yy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HJ_1_=configuration_t::hankel_transform(G_HJ_1, &args, quadl);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-0.5*sin(2.0*phi)*G_HJ_1_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HJ_0(r, r_);
        I_ = I_+DGF.yy;
    }
    return I_;
}

complex_t configuration_t::G_HJ_xz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HJ_3_=configuration_t::hankel_transform(G_HJ_3, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-((eta_0)/(j*this->k_0*this->layers[m].eps))*sin(phi)*G_HJ_3_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HJ_0(r, r_);
        I_ = I_+DGF.xz;
    }
    return I_;
}

complex_t configuration_t::G_HJ_yz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HJ_3_=configuration_t::hankel_transform(G_HJ_3, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((eta_0)/(j*this->k_0*this->layers[m].eps))*cos(phi)*G_HJ_3_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HJ_0(r, r_);
        I_ = I_+DGF.yz;
    }
    return I_;
}

complex_t configuration_t::G_HJ_zx(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HJ_4_=configuration_t::hankel_transform(G_HJ_4, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=((1.0)/(j*eta_0*this->k_0*this->layers[n].mu))*sin(phi)*G_HJ_4_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HJ_0(r, r_);
        I_ = I_+DGF.zx;
    }
    return I_;
}

complex_t configuration_t::G_HJ_zy(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    integrand_args_t args={rho, phi, r.z, r_.z, *this};
    complex_t G_HJ_4_=configuration_t::hankel_transform(G_HJ_4, &args, quadl);
    const complex_t j=complex_t(0.0, 1.0);
    const size_t n=this->find_layer(r.z);
    const size_t m=this->find_layer(r_.z);
    complex_t I_1=-((1.0)/(j*eta_0*this->k_0*this->layers[n].mu))*cos(phi)*G_HJ_4_;
    complex_t I_=I_1;
    if (m==n){
        Greens_functions_t DGF=configuration_t::G_HJ_0(r, r_);
        I_ = I_+DGF.zy;
    }
    return I_;
}

complex_t configuration_t::G_HJ_zz(const position_t r, const position_t r_, quadl_t quadl){
    const real_t rho=sqrt((r.x-r_.x)*(r.x-r_.x)+(r.y-r_.y)*(r.y-r_.y));
    const real_t phi=atan2(r.y-r_.y, r.x-r_.x);
    return 0.0*rho*phi*quadl.N;
}

near_field_t configuration_t::compute_H_J_near_field(const position_t r, const dipole_t dipole_J, 
    quadl_t quadl){
    near_field_t H;
    H.x = this->G_HJ_xx(r, dipole_J.r, quadl)*dipole_J.x+this->G_HJ_xy(r, dipole_J.r, quadl)*dipole_J.y+this->G_HJ_xz(r, dipole_J.r, quadl)*dipole_J.z;
    H.y = this->G_HJ_yx(r, dipole_J.r, quadl)*dipole_J.x+this->G_HJ_yy(r, dipole_J.r, quadl)*dipole_J.y+this->G_HJ_yz(r, dipole_J.r, quadl)*dipole_J.z;
    H.z = this->G_HJ_zx(r, dipole_J.r, quadl)*dipole_J.x+this->G_HJ_zy(r, dipole_J.r, quadl)*dipole_J.y+this->G_HJ_zz(r, dipole_J.r, quadl)*dipole_J.z;
    return H;
}


// plane wave

enum side_t{top, bottom};

near_field_plane_wave_t configuration_t::compute_plane_wave(const position_t r, const real_t theta_i, const real_t phi_i, 
    const incident_plane_wave_field_E_0_t E_0){
    side_t side;
    if (cos(theta_i)>=0){
        side = top;
    }else{
        side = bottom;
    }
    near_field_plane_wave_t fields;
    if (side==top){
        const real_t z_=this->layers[0].z_max;
        const complex_t eps_1=this->layers[0].eps;
        const complex_t mu_1=this->layers[0].mu;
        const complex_t eps=this->layers[find_layer(r.z)].eps;
        const complex_t mu=this->layers[find_layer(r.z)].mu;
        const complex_t k_1=this->layers[0].k;
        const complex_t k_rho_i=k_1*sin(theta_i);
        TLGF_t TLGF;
        TLGF = configuration_t::TLGF(k_rho_i, r.z, z_, v_source, sheet_I);
        const complex_t Z_e_1=this->layers[0].eta*cos(theta_i);
        const complex_t Z_h_1=this->layers[0].eta/cos(theta_i);
        vector_t<real_t> vector_rho_i, vector_phi_i, vector_z_i;
        vector_rho_i = vector_t<real_t>(cos(phi_i), sin(phi_i), 0.0);
        vector_phi_i = vector_t<real_t>(-sin(phi_i), cos(phi_i), 0.0);
        vector_z_i = vector_t<real_t>(0.0, 0.0, 1.0);
        const complex_t j=complex_t(0.0, +1.0);
        vector_t<complex_t> vector_k_rho_i=vector_t<complex_t>(k_rho_i*cos(phi_i), k_rho_i*sin(phi_i), 0.0);
        vector_t<real_t> vector_rho=vector_t<real_t>(r.x, r.y, 0.0);
        const complex_t factor=exp(+j*vector_k_rho_i*vector_rho)*exp(+j*k_1*cos(theta_i)*z_);
        vector_t<complex_t> E, H;
        E = -2.0*E_0.TM*(TLGF.V_e*cos(theta_i)*vector_rho_i-TLGF.I_e*Z_e_1*(eps_1/eps)*sin(theta_i)*vector_z_i)*factor+
            -2.0*E_0.TE*(TLGF.V_h*vector_phi_i)*factor;
        H = -2.0*E_0.TM*(TLGF.I_e*cos(theta_i)*vector_phi_i)*factor+
            -2.0*E_0.TE*(-TLGF.I_h*vector_rho_i+TLGF.V_h*(mu_1/(mu*Z_h_1))*tan(theta_i)*vector_z_i)*factor;
        fields.E.x = E.x;
        fields.E.y = E.y;
        fields.E.z = E.z;
        fields.H.x = H.x;
        fields.H.y = H.y;
        fields.H.z = H.z;
    }
    if (side==bottom){
        const real_t z_=this->layers[this->N-1].z_min;
        const complex_t eps_N=this->layers[this->N-1].eps;
        const complex_t mu_N=this->layers[this->N-1].mu;
        const complex_t eps=this->layers[find_layer(r.z)].eps;
        const complex_t mu=this->layers[find_layer(r.z)].mu;
        const complex_t k_N=this->layers[this->N-1].k;
        const complex_t k_rho_i=k_N*sin(theta_i);
        TLGF_t TLGF;
        TLGF = configuration_t::TLGF(k_rho_i, r.z, z_, v_source, sheet_I);
        const complex_t Z_e_N=this->layers[this->N-1].eta*cos(theta_i);
        const complex_t Z_h_N=this->layers[this->N-1].eta/cos(theta_i);
        vector_t<real_t> vector_rho_i, vector_phi_i, vector_z_i;
        vector_rho_i = vector_t<real_t>(cos(phi_i), sin(phi_i), 0.0);
        vector_phi_i = vector_t<real_t>(-sin(phi_i), cos(phi_i), 0.0);
        vector_z_i = vector_t<real_t>(0.0, 0.0, 1.0);
        const complex_t j=complex_t(0.0, +1.0);
        vector_t<complex_t> vector_k_rho_i=vector_t<complex_t>(k_rho_i*cos(phi_i), k_rho_i*sin(phi_i), 0.0);
        vector_t<real_t> vector_rho=vector_t<real_t>(r.x, r.y, 0.0);
        const complex_t factor=exp(+j*vector_k_rho_i*vector_rho)*exp(+j*k_N*cos(theta_i)*z_);
        vector_t<complex_t> E, H;
        E = -2.0*E_0.TM*(TLGF.V_e*cos(theta_i)*vector_rho_i-TLGF.I_e*Z_e_N*(eps_N/eps)*sin(theta_i)*vector_z_i)*factor+
            -2.0*E_0.TE*(TLGF.V_h*vector_phi_i)*factor;
        H = -2.0*E_0.TM*(TLGF.I_e*cos(theta_i)*vector_phi_i)*factor+
            -2.0*E_0.TE*(-TLGF.I_h*vector_rho_i+TLGF.V_h*(mu_N/(mu*Z_h_N))*tan(theta_i)*vector_z_i)*factor;
        fields.E.x = E.x;
        fields.E.y = E.y;
        fields.E.z = E.z;
        fields.H.x = H.x;
        fields.H.y = H.y;
        fields.H.z = H.z;
    }
    return fields;
}


