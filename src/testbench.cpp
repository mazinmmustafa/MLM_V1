//
#include "testbench.hpp"

void test_configuration(){

    const size_t N_layers=3;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    const complex_t j=complex_t(0.0, 1.0);
    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0.0*units::nm, +2000.0*units::nm, 1.0, 2.3013);
    config.add_layer(n++, -50.0*units::nm, +0.0*units::nm, 1.0, -11.753-j*1.2596);
    config.add_layer(n++, -2000.0*units::nm, -50.0*units::nm, 1.0, 1.0);

    config.log();

    const size_t Ns=1001;
    range_t theta;
    theta.set(0.0, pi/2.0, Ns);
    theta.linspace();
    file_t file;
    const complex_t k_1=config.layers[0].k;
    const complex_t k_N=config.layers[N_layers-1].k;
    const sheet_t sheet=sheet_I;
    file.open("data/GoldKretschmann/data.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        real_t theta_i=theta(i);
        complex_t k_rho;
        k_rho = k_1*sin(theta_i);
        Gamma_t Gamma_d=config.refl_d(k_rho, 0, sheet);
        k_rho = k_N*sin(theta_i);
        Gamma_t Gamma_u=config.refl_u(k_rho, N_layers-1, sheet);
        file.write("%21.14E %21.14E %21.14E %21.14E %21.14E\n", theta_i*180.0/pi, 
            abs(Gamma_d.Gamma_e)*abs(Gamma_d.Gamma_e), 
            abs(Gamma_d.Gamma_h)*abs(Gamma_d.Gamma_h),
            abs(Gamma_u.Gamma_e)*abs(Gamma_u.Gamma_e), 
            abs(Gamma_u.Gamma_h)*abs(Gamma_u.Gamma_h));
    }
    file.close();
    theta.unset();

    tau_t tau;
    tau = config.tau_u((1.2+j)*config.k_0, +1.0, +1.0, 2, 0, sheet_I);
    print(tau.tau_e);
    print(tau.tau_h);
    tau = config.tau_d((1.2+j)*config.k_0, +1.0, +1.0, 0, 2, sheet_I);
    print(tau.tau_e);
    print(tau.tau_h);

    print(config.find_layer(+2000.0*units::nm));
    print(config.find_layer(-0.0*units::nm));
    print(config.find_layer(-40.0*units::nm));
    print(config.find_layer(-2000.0*units::nm));

    //

    TLGF_t TLGF;
    TLGF = config.TLGF_r((1.2+j)*config.k_0, +20E-9, +25E-9, v_source, sheet_I);
    print(TLGF.V_e);
    print(TLGF.I_e);
    print(TLGF.V_h);
    print(TLGF.I_h);
    TLGF = config.TLGF_r((1.2+j)*config.k_0, -60E-9, -70E-9, i_source, sheet_I);
    print(TLGF.V_e);
    print(TLGF.I_e);
    print(TLGF.V_h);
    print(TLGF.I_h);
    print("\n");
    TLGF = config.TLGF_r((1.2+j)*config.k_0, +50E-9, -70E-9, i_source, sheet_I);
    print(TLGF.V_e);
    print(TLGF.I_e);
    print(TLGF.V_h);
    print(TLGF.I_h);
    TLGF = config.TLGF_r((1.2+j)*config.k_0, -60E-9, +70E-9, i_source, sheet_I);
    print(TLGF.V_e);
    print(TLGF.I_e);
    print(TLGF.V_h);
    print(TLGF.I_h);
    config.unset();
}


void test_TLGFs(){

    const size_t N_layers=4;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0500.0*units::nm, +1000.0*units::nm, 1.0, 01.0);
    config.add_layer(n++, +0000.0*units::nm, +0500.0*units::nm, 1.0, 02.0);
    config.add_layer(n++, -0500.0*units::nm, +0000.0*units::nm, 1.0, 10.0);
    config.add_layer(n++, -1000.0*units::nm, -0500.0*units::nm, 1.0, 01.0);

    config.log();

    const size_t Ns=1001;
    const real_t z_=+750*units::nm;
    range_t z;
    z.set(-1000.0*units::nm, +1000.0*units::nm, Ns);
    z.linspace();
    const complex_t k_rho=0.7*config.k_0;
    const sheet_t sheet=sheet_I;
    TLGF_t TLGF;
    file_t file_v, file_i;
    file_v.open("data/Paulus/data_v.txt", 'w');
    file_i.open("data/Paulus/data_i.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        file_v.write("%21.14E ", z(i)/units::nm);
        file_i.write("%21.14E ", z(i)/units::nm);
        // TLGF = config.TLGF(k_rho, z(i), z_, v_source, sheet);
        TLGF = config.TLGF_r(k_rho, z(i), z_, v_source, sheet);
        file_v.write("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E %21.14E %21.14E\n", 
            real(TLGF.V_e), imag(TLGF.V_e), 
            real(TLGF.I_e), imag(TLGF.I_e),
            real(TLGF.V_h), imag(TLGF.V_h),
            real(TLGF.I_h), imag(TLGF.I_h));
        // TLGF = config.TLGF(k_rho, z(i), z_, i_source, sheet);
        TLGF = config.TLGF_r(k_rho, z(i), z_, i_source, sheet);
        file_i.write("%21.14E %21.14E %21.14E %21.14E %21.14E %21.14E %21.14E %21.14E\n", 
            real(TLGF.V_e), imag(TLGF.V_e), 
            real(TLGF.I_e), imag(TLGF.I_e),
            real(TLGF.V_h), imag(TLGF.V_h),
            real(TLGF.I_h), imag(TLGF.I_h));
    }
    file_v.close();
    file_i.close();
    z.unset();

    print(config.k_min);

    config.unset();
}


struct integrand_args_t{
    real_t rho, phi, z, z_, v;
    configuration_t &config;
};

complex_t integrand_I(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    TLGF_t TLGF=args->config.TLGF(k_rho*args->config.k_0, args->z, args->z_, i_source, sheet_I);
    return (TLGF.V_e+TLGF.V_h)*besselj(args->v, k_rho*args->config.k_0*args->rho)*k_rho*args->config.k_0;
}

complex_t integrand_II(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, args->rho, abs(args->z-args->z_));
    TLGF_t TLGF=args->config.TLGF(k_rho_, args->z, args->z_, i_source, sheet_I);
    return factor*(TLGF.V_e+TLGF.V_h)*besselj(args->v, k_rho_*args->rho)*k_rho_;
}

complex_t integrand_III(const complex_t k_rho, void *args_){
    integrand_args_t *args=(integrand_args_t*)args_;
    complex_t k_rho_=k_rho;
    complex_t factor=args->config.detour(k_rho_, args->rho, abs(args->z-args->z_));
    TLGF_t TLGF=args->config.TLGF_r(k_rho_, args->z, args->z_, i_source, sheet_I);
    return factor*(TLGF.V_e+TLGF.V_h)*besselj(args->v, k_rho_*args->rho)*k_rho_;
}

void test_integrands(){

    const size_t N_layers=4;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0500.0*units::nm, +1000.0*units::nm, 1.0, 01.0);
    config.add_layer(n++, +0000.0*units::nm, +0500.0*units::nm, 1.0, 02.0);
    config.add_layer(n++, -0500.0*units::nm, +0000.0*units::nm, 1.0, 10.0);
    config.add_layer(n++, -1000.0*units::nm, -0500.0*units::nm, 1.0, 01.0);

    config.log();

    const size_t Ns=10001;
    const real_t z=+750*units::nm;
    const real_t z_=+750*units::nm;
    const real_t phi=45.0*pi/180.0;
    const real_t rho=633.0*units::nm;
    const real_t v=0.0;

    range_t k_rho;
    k_rho.set(0.0, 6.0, Ns);
    k_rho.linspace();

    file_t file;
    integrand_args_t args={rho, phi, z, z_, v, config};
    file.open("data/Paulus/data_I.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        file.write("%21.14E %21.14E ", k_rho(i), abs(integrand_I(k_rho(i), &args)));
        file.write("%21.14E %21.14E\n", abs(integrand_II(k_rho(i), &args)), abs(integrand_III(k_rho(i), &args)));
    }
    file.close();

    k_rho.unset();

    config.unset();
}

void test_DGFs_Paulus(){

    const size_t N_layers=4;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0500.0*units::nm, +10000.0*units::nm, 1.0, 1.0);
    config.add_layer(n++, +0000.0*units::nm, +0500.0*units::nm, 1.0, 2.0);
    config.add_layer(n++, -0500.0*units::nm, +0000.0*units::nm, 1.0, 10.0);
    config.add_layer(n++, -10000.0*units::nm, -0500.0*units::nm, 1.0, 1.0);

    config.log();

    //
    quadl_t quadl;
    const size_t N_quadl=32;
    const size_t k_max=15;
    const real_t tol=1.0E-4;
    quadl.set(N_quadl, k_max, tol);

    const size_t Ns=401;
    const real_t z_=+750*units::nm;
    const real_t phi=+45.0*pi/180.0;
    const real_t rho=633.0*units::nm;

    range_t z;
    const real_t z_min=-1000.0*units::nm;
    const real_t z_max=+1000.0*units::nm;
    z.set(z_min, z_max, Ns);
    z.linspace();

    file_t file;
    file.open("data/Paulus/data_EJ.txt", 'w');
    stopwatch_t timer;
    timer.set();
    for (size_t i=0; i<Ns; i++){
        file.write("%21.14E ", z(i));
        position_t r=cylindrical_t(rho, phi, z(i));
        position_t r_=cylindrical_t(0.0, 0.0, z_);
        Greens_functions_t DGF;
        DGF.xx = config.G_EJ_xx(r, r_, quadl);
        DGF.zx = config.G_EJ_zx(r, r_, quadl);
        DGF.xz = config.G_EJ_xz(r, r_, quadl);
        DGF.zz = config.G_EJ_zz(r, r_, quadl);
        file.write("%21.14E ", abs(DGF.xx));
        file.write("%21.14E ", abs(DGF.zx));
        file.write("%21.14E ", abs(DGF.xz));
        file.write("%21.14E ", abs(DGF.zz));
        file.write("\n");
    }
    timer.unset();
    file.close();

    z.unset();
    quadl.unset();
    config.unset();

}

void test_DGFs_Chew(){

    const size_t N_layers=7;
    configuration_t config;

    const real_t lambda_0=1.0*units::m;
    const real_t freq=c_0/lambda_0;

    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0.0*units::m, +10.0*units::m, 1.0, 1.0);
    config.add_layer(n++, -0.2*units::m, +0.0*units::m, 1.0, 2.6);
    config.add_layer(n++, -0.5*units::m, -0.2*units::m, 3.2, 6.5);
    config.add_layer(n++, -1.0*units::m, -0.5*units::m, 6.0, 4.2);
    config.add_layer(n++, -1.3*units::m, -1.0*units::m, 3.2, 6.5);
    config.add_layer(n++, -1.5*units::m, -1.3*units::m, 1.0, 2.6);
    config.add_layer(n++, -10.0*units::m, -1.5*units::m, 1.0, 1.0);

    const position_t r_=cartesian_t(0.0*units::m, 0.0*units::m, -1.4*units::m);
    const real_t theta_0=+20.0*pi/180.0;
    const real_t phi_0=+30.0*pi/180.0;
    const complex_t Il=+1.0;
    const complex_t Kl=+1.0;
    dipole_t J=dipole_t(r_, theta_0, phi_0, Il);
    dipole_t M=dipole_t(r_, theta_0, phi_0, Kl);
    const real_t y=+1.0*units::m;
    const real_t z=-0.3*units::m;

    config.log();

    //
    quadl_t quadl;
    const size_t N_quadl=32;
    const size_t k_max=15;
    const real_t tol=1.0E-4;
    quadl.set(N_quadl, k_max, tol);

    const size_t Ns=401;

    range_t x;
    const real_t x_min=-3.0*units::m;
    const real_t x_max=+3.0*units::m;
    x.set(x_min, x_max, Ns);
    x.linspace();

    file_t file_EJ, file_EM;
    file_t file_HJ, file_HM;
    file_EJ.open("data/Chew/data_EJ.txt", 'w');
    file_EM.open("data/Chew/data_EM.txt", 'w');
    file_HJ.open("data/Chew/data_HJ.txt", 'w');
    file_HM.open("data/Chew/data_HM.txt", 'w');
    stopwatch_t timer;
    timer.set();
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing E fields...");
        file_EJ.write("%21.14E ", x(i));
        file_EM.write("%21.14E ", x(i));
        file_HJ.write("%21.14E ", x(i));
        file_HM.write("%21.14E ", x(i));
        position_t r=cartesian_t(x(i), y, z);
        near_field_t E;
        E = config.compute_E_J_near_field(r, J, quadl);
        file_EJ.write("%21.14E ", abs(E.x));
        file_EJ.write("%21.14E ", abs(E.y));
        file_EJ.write("%21.14E ", abs(E.z));
        file_EJ.write("\n");
        E = config.compute_E_M_near_field(r, M, quadl);
        file_EM.write("%21.14E ", abs(E.x));
        file_EM.write("%21.14E ", abs(E.y));
        file_EM.write("%21.14E ", abs(E.z));
        file_EM.write("\n");
        progress_bar(i, Ns, "computing H fields...");
        near_field_t H;
        H = config.compute_H_J_near_field(r, J, quadl);
        file_HJ.write("%21.14E ", abs(H.x));
        file_HJ.write("%21.14E ", abs(H.y));
        file_HJ.write("%21.14E ", abs(H.z));
        file_HJ.write("\n");
        H = config.compute_H_M_near_field(r, M, quadl);
        file_HM.write("%21.14E ", abs(H.x));
        file_HM.write("%21.14E ", abs(H.y));
        file_HM.write("%21.14E ", abs(H.z));
        file_HM.write("\n");
    }
    timer.unset();
    file_EJ.close();
    file_EM.close();
    file_HJ.close();
    file_HM.close();

    x.unset();
    quadl.unset();
    config.unset();

}


void test_DGFs_Gold_Kretschmann(){

    const size_t N_layers=3;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    const complex_t j=complex_t(0.0, +1.0);
    config.set(N_layers, freq);
    size_t n=0;
    // Nadson
    // config.add_layer(n++, +0.0*units::nm, +2000.0*units::nm, 1.0, 3.0615);
    // config.add_layer(n++, -50.0*units::nm, +0.0*units::nm, 1.0, -11.67-j*1.35);
    // config.add_layer(n++, -2000.0*units::nm, -50.0*units::nm, 1.0, 1.0);
    // Gold Kretschmann
    config.add_layer(n++, +0.0*units::nm, +2000.0*units::nm, 1.0, 2.3013);
    config.add_layer(n++, -50.0*units::nm, +0.0*units::nm, 1.0, -11.753-j*1.2596);
    config.add_layer(n++, -2000.0*units::nm, -50.0*units::nm, 1.0, 1.0);

    const position_t r_=cartesian_t(0.0*units::nm, 0.0*units::nm, +20.0*units::nm);
    const real_t theta_0=+0.0*pi/180.0;
    const real_t phi_0=+0.0*pi/180.0;
    const complex_t Il=+1.0;
    dipole_t J=dipole_t(r_, theta_0, phi_0, Il);
    const real_t y=+0.0*units::nm;

    config.log();

    //
    quadl_t quadl;
    const size_t N_quadl=8;
    const size_t k_max=15;
    const real_t tol=1.0E-4;
    quadl.set(N_quadl, k_max, tol);

    const size_t Ns=201;
    const size_t Nx=Ns, Nz=Ns;

    range_t x, z;
    const real_t x_min=-2000.0*units::nm;
    const real_t x_max=+2000.0*units::nm;
    const real_t z_min=-2000.0*units::nm;
    const real_t z_max=+2000.0*units::nm;
    x.set(x_min, x_max, Nx);
    z.set(z_min, z_max, Nz);
    x.linspace();
    z.linspace();

    file_t file_data;
    file_data.open("data/GoldKretschmann/data_VED.txt", 'w');
    stopwatch_t timer;
    timer.set();
    for (size_t ii=0; ii<Nx; ii++){
        progress_bar(ii, Nx, "computing E fields...");
        for (size_t jj=0; jj<Nz; jj++){
            position_t r=cartesian_t(x(ii), y, z(jj));
            near_field_t E;
            E = config.compute_E_J_near_field(r, J, quadl);
            // const complex_t E_mag=sqrt(abs(E.x*E.x)+abs(E.y*E.y)+abs(E.z*E.z));
            // const real_t E_mag=sqrt(pow(real(E.z), 2.0));
            const real_t E_mag=sqrt(pow(real(E.x), 2.0)+pow(real(E.y), 2.0)+pow(real(E.z), 2.0));

            file_data.write("%21.14E %21.14E %21.14E\n", x(ii)/1.0E-6, z(jj)/1.0E-6, 20.0*log10(E_mag)); // Gnu_plot
        }
        file_data.write("\n");
    }
    timer.unset();
    file_data.close();

    x.unset();
    z.unset();


    quadl.unset();
    config.unset();

}

void test_DGFs_Gold_Kretschmann_cut(){

    const size_t N_layers=3;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    const complex_t j=complex_t(0.0, +1.0);
    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0.0*units::nm, +2500.0*units::nm, 1.0, 2.3013);
    config.add_layer(n++, -50.0*units::nm, +0.0*units::nm, 1.0, -11.753-j*1.2596);
    config.add_layer(n++, -2500.0*units::nm, -50.0*units::nm, 1.0, 1.0);
;

    const position_t r_=cartesian_t(0.0*units::nm, 0.0*units::nm, +20.0*units::nm);
    const real_t theta_0=+0.0*pi/180.0;
    const real_t phi_0=+0.0*pi/180.0;
    const complex_t Il=+1.0;
    dipole_t J=dipole_t(r_, theta_0, phi_0, Il);
    const real_t x=+0.0*units::nm;
    const real_t y=+0.0*units::nm;

    config.log();

    //
    quadl_t quadl;
    const size_t N_quadl=32;
    const size_t k_max=15;
    const real_t tol=1.0E-4;
    quadl.set(N_quadl, k_max, tol);

    const size_t Ns=1001;

    range_t z;
    const real_t z_min=-400.0*units::nm;
    const real_t z_max=+400.0*units::nm;
    z.set(z_min, z_max, Ns);
    z.linspace();

    file_t file;
    file.open("data/GoldKretschmann/data_cut.txt", 'w');
    stopwatch_t timer;
    timer.set();
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing E fields...");
        file.write("%21.14E ", z(i));
        position_t r=cartesian_t(x, y, z(i));
        near_field_t E;
        E = config.compute_E_J_near_field(r, J, quadl);
        file.write("%21.14E %21.14E %21.14E\n", abs(E.x), abs(E.y), abs(E.z));
        file.write("\n");
    }
    timer.unset();
    file.close();

    z.unset();
    
    quadl.unset();
    config.unset();

}


void test_DGFs_Paulus_near_field(){

    const size_t N_layers=4;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0500.0*units::nm, +10000.0*units::nm, 1.0, 1.0);
    config.add_layer(n++, +0000.0*units::nm, +0500.0*units::nm, 1.0, 2.0);
    config.add_layer(n++, -0500.0*units::nm, +0000.0*units::nm, 1.0, 10.0);
    config.add_layer(n++, -10000.0*units::nm, -0500.0*units::nm, 1.0, 1.0);
;

    const position_t r_=cartesian_t(0.0*units::nm, 0.0*units::nm, +200.0*units::nm);
    const real_t theta_0=+0.0*pi/180.0;
    const real_t phi_0=+0.0*pi/180.0;
    const complex_t Il=+1.0;
    dipole_t J=dipole_t(r_, theta_0, phi_0, Il);
    const real_t y=+0.0*units::nm;

    config.log();

    //
    quadl_t quadl;
    const size_t N_quadl=32;
    const size_t k_max=15;
    const real_t tol=1.0E-4;
    quadl.set(N_quadl, k_max, tol);

    const size_t Ns=201;
    const size_t Nx=Ns, Nz=Ns;

    range_t x, z;
    const real_t x_min=-2000.0*units::nm;
    const real_t x_max=+2000.0*units::nm;
    const real_t z_min=-2000.0*units::nm;
    const real_t z_max=+2000.0*units::nm;
    x.set(x_min, x_max, Nx);
    z.set(z_min, z_max, Nz);
    x.linspace();
    z.linspace();

    file_t file_x, file_z, file_data;
    file_x.open("data/Paulus/data_x.txt", 'w');
    file_z.open("data/Paulus/data_z.txt", 'w');
    file_data.open("data/Paulus/data.txt", 'w');
    stopwatch_t timer;
    timer.set();
    for (size_t ii=0; ii<Nx; ii++){
        progress_bar(ii, Nx, "computing E fields...");
        for (size_t jj=0; jj<Nz; jj++){
            file_x.write("%21.14E ", x(ii));
            file_z.write("%21.14E ", z(jj));
            position_t r=cartesian_t(x(ii), y, z(jj));
            near_field_t E;
            E = config.compute_E_J_near_field(r, J, quadl);
            // const complex_t E_mag=sqrt(abs(E.x*E.x)+abs(E.y*E.y)+abs(E.z*E.z));
            const complex_t E_mag=sqrt(pow(real(E.x), 2.0)+pow(real(E.y), 2.0)+pow(real(E.z), 2.0));
            file_data.write("%21.14E ", E_mag);
        }
        file_x.write("\n");
        file_z.write("\n");
        file_data.write("\n");
    }
    timer.unset();
    file_x.close();
    file_z.close();
    file_data.close();

    x.unset();
    z.unset();

    quadl.unset();
    config.unset();

}


void test_planar_wave_Chew_cut(){

    const size_t N_layers=7;
    configuration_t config;

    const real_t lambda_0=1.0*units::m;
    const real_t freq=c_0/lambda_0;

    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0.0*units::m, +10.0*units::m, 1.0, 1.0);
    config.add_layer(n++, -0.2*units::m, +0.0*units::m, 1.0, 2.6);
    config.add_layer(n++, -0.5*units::m, -0.2*units::m, 3.2, 6.5);
    config.add_layer(n++, -1.0*units::m, -0.5*units::m, 6.0, 4.2);
    config.add_layer(n++, -1.3*units::m, -1.0*units::m, 3.2, 6.5);
    config.add_layer(n++, -1.5*units::m, -1.3*units::m, 1.0, 2.6);
    config.add_layer(n++, -10.0*units::m, -1.5*units::m, 1.0, 1.0);

    const real_t theta_i=+30.0*pi/180.0;
    const real_t phi_i=+160.0*pi/180.0;
    const real_t psi_i=(+90.0+45.0)*pi/180.0;
    incident_plane_wave_field_E_0_t E_0;
    E_0.TM = 1.0*cos(psi_i);
    E_0.TE = 1.0*sin(psi_i);

    const real_t x=+0.0*units::m;
    const real_t y=+0.0*units::m;

    config.log();

    //
    const size_t Ns=4001;

    range_t z;
    const real_t z_min=-2.0*units::m;
    const real_t z_max=+1.0*units::m;
    z.set(z_min, z_max, Ns);
    z.linspace();

    file_t file;
    file.open("data/Chew/data_PW.txt", 'w');
    stopwatch_t timer;
    timer.set();
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing plane wave fields...");
        file.write("%21.14E ", z(i));
        position_t r=cartesian_t(x, y, z(i));
        near_field_plane_wave_t fields;
        fields = config.compute_plane_wave(r, theta_i, phi_i, E_0);
        file.write("%21.14E %21.14E %21.14E ", abs(fields.E.x), abs(fields.E.y), abs(fields.E.z));
        file.write("%21.14E %21.14E %21.14E ", abs(fields.H.x), abs(fields.H.y), abs(fields.H.z));
        file.write("\n");
    }
    timer.unset();
    file.close();

    z.unset();
    config.unset();

}


void test_Gold_Kretschmann_reflection(){

    const size_t N_layers=3;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    const complex_t j=complex_t(0.0, +1.0);
    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0.0*units::nm, +2500.0*units::nm, 1.0, 2.3013);
    config.add_layer(n++, -50.0*units::nm, +0.0*units::nm, 1.0, -11.753-j*1.2596);
    config.add_layer(n++, -2500.0*units::nm, -50.0*units::nm, 1.0, 1.0);

    config.log();

    //
    const size_t Ns=1001;

    range_t theta_i;

    const real_t theta_i_min=deg_to_rad(+0.0);
    const real_t theta_i_max=deg_to_rad(+180.0);
    theta_i.set(theta_i_min, theta_i_max, Ns);
    theta_i.linspace();

    file_t file;
    file.open("data/GoldKretschmann/data_refl.txt", 'w');
    
    stopwatch_t timer;
    timer.set();
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing reflection coefficients...");
        Gamma_t Gamma=config.compute_Gamma(theta_i(i), sheet_I);
        file.write("%21.14E %21.14E %21.14E\n", theta_i(i)*180.0/pi, 
            pow(abs(Gamma.Gamma_e), 2.0), pow(abs(Gamma.Gamma_h), 2.0));
    }
    timer.unset();
    file.close();

    theta_i.unset();
    config.unset();

}

void test_planar_wave_Gold_Kretschmann_near_field(){

    const size_t N_layers=3;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    const complex_t j=complex_t(0.0, +1.0);
    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0.0*units::nm, +2500.0*units::nm, 1.0, 2.3013);
    config.add_layer(n++, -50.0*units::nm, +0.0*units::nm, 1.0, -11.753-j*1.2596);
    config.add_layer(n++, -2500.0*units::nm, -50.0*units::nm, 1.0, 1.0);

    const real_t theta_i=deg_to_rad(+43.7);
    const real_t phi_i=deg_to_rad(+180.0);
    incident_plane_wave_field_E_0_t E_0;
    E_0.TM = 1.0;
    E_0.TE = 0.0;

    const real_t y=+0.0*units::nm;

    config.log();

    //
    const size_t Nx=1001, Nz=1001;

    range_t x, z;

    const real_t x_min=-1000.0*units::nm;
    const real_t x_max=+1000.0*units::nm;
    const real_t z_min=-1000.0*units::nm;
    const real_t z_max=+1000.0*units::nm;
    x.set(x_min, x_max, Nx);
    z.set(z_min, z_max, Nz);
    x.linspace();
    z.linspace();


    file_t file_E, file_H;
    file_E.open("data/GoldKretschmann/data_PW_E.txt", 'w');
    file_H.open("data/GoldKretschmann/data_PW_H.txt", 'w');
    stopwatch_t timer;
    timer.set();
    for (size_t ii=0; ii<Nx; ii++){
        progress_bar(ii, Nx, "computing plane wave fields...");
        for (size_t jj=0; jj<Nz; jj++){
            position_t r=cartesian_t(x(ii), y, z(jj));
            near_field_plane_wave_t fields;
            fields = config.compute_plane_wave(r, theta_i, phi_i, E_0);
            const real_t E_mag=sqrt(pow(real(fields.E.x), 2.0)+pow(real(fields.E.y), 2.0)+pow(real(fields.E.z), 2.0));
            const real_t H_mag=sqrt(pow(real(fields.H.x), 2.0)+pow(real(fields.H.y), 2.0)+pow(real(fields.H.z), 2.0));
            file_E.write("%21.14E %21.14E %21.14E\n", x(ii)/1E-6, z(jj)/1E-6, E_mag);
            file_H.write("%21.14E %21.14E %21.14E\n", x(ii)/1E-6, z(jj)/1E-6, H_mag/1.0E-3);
        }
        file_E.write("\n");
        file_H.write("\n");
    }
    timer.unset();
    file_E.close();

    x.unset();
    z.unset();
    config.unset();

}


void test_plasmonic_WG_far_field(){

    const size_t N_layers=19;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    const complex_t j=complex_t(0.0, +1.0);
    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +1523.0*units::nm, +2500.0*units::nm, 1.0, 2.3013);
    config.add_layer(n++, +1445.0*units::nm, +1523.0*units::nm, 1.0, 4.5967);
    config.add_layer(n++, +1329.0*units::nm, +1445.0*units::nm, 1.0, 2.1229);
    config.add_layer(n++, +1241.0*units::nm, +1329.0*units::nm, 1.0, 4.5967);
    config.add_layer(n++, +1115.0*units::nm, +1241.0*units::nm, 1.0, 2.1229);
    config.add_layer(n++, +1037.0*units::nm, +1115.0*units::nm, 1.0, 4.5967);
    config.add_layer(n++, +0911.0*units::nm, +1037.0*units::nm, 1.0, 2.1229);
    config.add_layer(n++, +0833.0*units::nm, +0911.0*units::nm, 1.0, 4.5967);
    config.add_layer(n++, +0707.0*units::nm, +0833.0*units::nm, 1.0, 2.1229);
    config.add_layer(n++, +0629.0*units::nm, +0707.0*units::nm, 1.0, 4.5967);
    config.add_layer(n++, +0503.0*units::nm, +0629.0*units::nm, 1.0, 2.1229);
    config.add_layer(n++, +0425.0*units::nm, +0503.0*units::nm, 1.0, 4.5967);
    config.add_layer(n++, +0299.0*units::nm, +0425.0*units::nm, 1.0, 2.1229);
    config.add_layer(n++, +0221.0*units::nm, +0299.0*units::nm, 1.0, 4.5967);
    config.add_layer(n++, +0069.0*units::nm, +0221.0*units::nm, 1.0, 2.1229);
    config.add_layer(n++, +0042.0*units::nm, +0069.0*units::nm, 1.0, 2.1229);
    config.add_layer(n++, +0000.0*units::nm, +0042.0*units::nm, 1.0, -18.3511-j*0.4331);
    config.add_layer(n++, -0027.0*units::nm, +0000.0*units::nm, 1.0, 2.1229);
    config.add_layer(n++, -0500.0*units::nm, -0027.0*units::nm, 1.0, 1.0);

    const position_t r_=cartesian_t(0.0*units::nm, 0.0*units::nm, -15.0*units::nm);
    const real_t theta_0=deg_to_rad(+135.0);
    const real_t phi_0=deg_to_rad(+0.0);
    const complex_t Il=+1.0;
    dipole_t J=dipole_t(r_, theta_0, phi_0, Il);

    real_t phi_s;

    config.log();

    //
    const size_t Ns=4001;

    range_t theta_s;

    const real_t theta_s_min=deg_to_rad(-180.0);
    const real_t theta_s_max=deg_to_rad(+180.0);
    theta_s.set(theta_s_min, theta_s_max, Ns);
    theta_s.linspace();


    file_t file;
    file.open("data/plasmonicWG/data_far_field.txt", 'w');
    stopwatch_t timer;
    timer.set();
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing far fields...");
        const real_t r=1000.0*lambda_0;
        file.write("%21.14E ", theta_s(i));
        far_field_t field;
        phi_s = deg_to_rad(0.0);
        field = config.compute_far_field_J(r, theta_s(i), phi_s, J);
        file.write("%21.14E ", 20.0*log10(abs(field.E.theta)));
        phi_s = deg_to_rad(90.0);
        field = config.compute_far_field_J(r, theta_s(i), phi_s, J);
        file.write("%21.14E ", 20.0*log10(abs(field.E.theta)));
        file.write("\n");
    }
    timer.unset();
    file.close();

    theta_s.unset();
    config.unset();

}

void test_Chew_far_field(){

    const size_t N_layers=7;
    configuration_t config;

    const real_t lambda_0=1.0*units::m;
    const real_t freq=c_0/lambda_0;

    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0.0*units::m, +10.0*units::m, 1.0, 1.0);
    config.add_layer(n++, -0.2*units::m, +0.0*units::m, 1.0, 2.6);
    config.add_layer(n++, -0.5*units::m, -0.2*units::m, 3.2, 6.5);
    config.add_layer(n++, -1.0*units::m, -0.5*units::m, 6.0, 4.2);
    config.add_layer(n++, -1.3*units::m, -1.0*units::m, 3.2, 6.5);
    config.add_layer(n++, -1.5*units::m, -1.3*units::m, 1.0, 2.6);
    config.add_layer(n++, -10.0*units::m, -1.5*units::m, 1.0, 1.0);

    const position_t r_=cartesian_t(0.0*units::m, 0.0*units::m, -1.4*units::m);
    const real_t theta_0=+20.0*pi/180.0;
    const real_t phi_0=+30.0*pi/180.0;
    const complex_t Il=+1.0;
    const complex_t Kl=+1.0;
    dipole_t J=dipole_t(r_, theta_0, phi_0, Il);
    dipole_t M=dipole_t(r_, theta_0, phi_0, Kl);

    real_t phi_s;

    config.log();

    //
    const size_t Ns=4001;

    range_t theta_s;

    const real_t theta_s_min=deg_to_rad(-180.0);
    const real_t theta_s_max=deg_to_rad(+180.0);
    theta_s.set(theta_s_min, theta_s_max, Ns);
    theta_s.linspace();


    file_t file;
    stopwatch_t timer;
    timer.set();
    file.open("data/Chew/data_far_field_J.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing far fields for J...");
        const real_t r=1000.0*lambda_0;
        file.write("%21.14E ", theta_s(i));
        far_field_t field;
        phi_s = deg_to_rad(0.0);
        field = config.compute_far_field_J(r, theta_s(i), phi_s, J);
        file.write("%21.14E ", 20.0*log10(abs(field.E.theta)));
        phi_s = deg_to_rad(90.0);
        field = config.compute_far_field_J(r, theta_s(i), phi_s, J);
        file.write("%21.14E ", 20.0*log10(abs(field.E.theta)));
        phi_s = deg_to_rad(0.0);
        field = config.compute_far_field_J(r, theta_s(i), phi_s, J);
        file.write("%21.14E ", 20.0*log10(abs(field.E.phi)));
        phi_s = deg_to_rad(90.0);
        field = config.compute_far_field_J(r, theta_s(i), phi_s, J);
        file.write("%21.14E ", 20.0*log10(abs(field.E.phi)));
        file.write("\n");
    }
    file.close();

    file.open("data/Chew/data_far_field_M.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing far fields for M...");
        const real_t r=1.0*lambda_0;
        file.write("%21.14E ", theta_s(i));
        far_field_t field;
        phi_s = deg_to_rad(0.0);
        field = config.compute_far_field_M(r, theta_s(i), phi_s, M);
        file.write("%21.14E ", 20.0*log10(abs(field.E.theta)));
        phi_s = deg_to_rad(90.0);
        field = config.compute_far_field_M(r, theta_s(i), phi_s, M);
        file.write("%21.14E ", 20.0*log10(abs(field.E.theta)));
        phi_s = deg_to_rad(0.0);
        field = config.compute_far_field_M(r, theta_s(i), phi_s, M);
        file.write("%21.14E ", 20.0*log10(abs(field.E.phi)));
        phi_s = deg_to_rad(90.0);
        field = config.compute_far_field_M(r, theta_s(i), phi_s, M);
        file.write("%21.14E ", 20.0*log10(abs(field.E.phi)));
        file.write("\n");
    }
    file.close();

    timer.unset();

    theta_s.unset();
    config.unset();

}


void test_Chew_far_field_2_elements(){

    const size_t N_layers=7;
    configuration_t config;

    const real_t lambda_0=1.0*units::m;
    const real_t freq=c_0/lambda_0;

    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0.0*units::m, +10.0*units::m, 1.0, 1.0);
    config.add_layer(n++, -0.2*units::m, +0.0*units::m, 1.0, 2.6);
    config.add_layer(n++, -0.5*units::m, -0.2*units::m, 3.2, 6.5);
    config.add_layer(n++, -1.0*units::m, -0.5*units::m, 6.0, 4.2);
    config.add_layer(n++, -1.3*units::m, -1.0*units::m, 3.2, 6.5);
    config.add_layer(n++, -1.5*units::m, -1.3*units::m, 1.0, 2.6);
    config.add_layer(n++, -10.0*units::m, -1.5*units::m, 1.0, 1.0);

    const position_t r_1=cartesian_t(+0.5*units::m, 0.0*units::m, -0.6*units::m);
    const position_t r_2=cartesian_t(-0.5*units::m, 0.0*units::m, -0.6*units::m);
    const complex_t Il=+1.0;
    dipole_t J_1=dipole_t(r_1, deg_to_rad(+0.0), deg_to_rad(+0.0), Il);
    dipole_t J_2=dipole_t(r_2, deg_to_rad(+0.0), deg_to_rad(+0.0), Il);

    real_t phi_s;

    config.log();

    //
    const size_t Ns=4001;

    range_t theta_s;

    const real_t theta_s_min=deg_to_rad(-180.0);
    const real_t theta_s_max=deg_to_rad(+180.0);
    theta_s.set(theta_s_min, theta_s_max, Ns);
    theta_s.linspace();


    file_t file;
    stopwatch_t timer;
    timer.set();
    file.open("data/Chew/data_far_field_J_2.txt", 'w');
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing far fields for J...");
        const real_t r=1000.0*lambda_0;
        file.write("%21.14E ", theta_s(i));
        far_field_t field_1, field_2;
        phi_s = deg_to_rad(0.0);
        field_1 = config.compute_far_field_J(r, theta_s(i), phi_s, J_1);
        field_2 = config.compute_far_field_J(r, theta_s(i), phi_s, J_2);
        file.write("%21.14E ", 20.0*log10(abs(field_1.E.theta+field_2.E.theta)));
        phi_s = deg_to_rad(90.0);
        field_1 = config.compute_far_field_J(r, theta_s(i), phi_s, J_1);
        field_2 = config.compute_far_field_J(r, theta_s(i), phi_s, J_2);
        file.write("%21.14E ", 20.0*log10(abs(field_1.E.theta+field_2.E.theta)));
        file.write("\n");
    }
    file.close();

    timer.unset();

    theta_s.unset();
    config.unset();

}


void test_DGFs_Gold_Kretschmann_near_far_fields(){

    const size_t N_layers=3;
    configuration_t config;

    const real_t lambda_0=633.0*units::nm;
    const real_t freq=c_0/lambda_0;

    const complex_t j=complex_t(0.0, +1.0);
    config.set(N_layers, freq);
    size_t n=0;
    // Gold Kretschmann
    config.add_layer(n++, +0.0*units::nm, +2000.0*units::nm, 1.0, 2.3013);
    config.add_layer(n++, -50.0*units::nm, +0.0*units::nm, 1.0, -11.753-j*1.2596);
    config.add_layer(n++, -2000.0*units::nm, -50.0*units::nm, 1.0, 1.0);

    const position_t r_=cartesian_t(0.0*units::nm, 0.0*units::nm, +20.0*units::nm);
    const real_t theta_0=+0.0*pi/180.0;
    const real_t phi_0=+0.0*pi/180.0;
    const complex_t Il=+1.0;
    dipole_t J=dipole_t(r_, theta_0, phi_0, Il);
    const real_t r=+100.0*lambda_0;
    
    config.log();

    //
    quadl_t quadl;
    const size_t N_quadl=32;
    const size_t k_max=15;
    const real_t tol=1.0E-4;
    quadl.set(N_quadl, k_max, tol);

    const size_t Ns=1001;

    range_t theta;
    const real_t theta_min=deg_to_rad(-180.0);
    const real_t theta_max=deg_to_rad(+180.0);
    theta.set(theta_min, theta_max, Ns);
    theta.linspace();

    file_t file_J, file_M;
    file_J.open("data/GoldKretschmann/data_near_far_J.txt", 'w');
    file_M.open("data/GoldKretschmann/data_near_far_M.txt", 'w');
    stopwatch_t timer;
    timer.set();
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "computing E fields...");
        position_t r_position;
        real_t phi;
        phi = deg_to_rad(0.0);
        r_position = cartesian_t(r*sin(theta(i))*cos(phi), 
            r*sin(theta(i))*sin(phi), r*cos(theta(i)));
        near_field_t E;
        E = config.compute_E_J_near_field(r_position, J, quadl);
        const real_t E_mag_J=sqrt(abs(E.x*E.x)+abs(E.y*E.y)+abs(E.z*E.z));
        E = config.compute_E_M_near_field(r_position, J, quadl);
        const real_t E_mag_M=sqrt(abs(E.x*E.x)+abs(E.y*E.y)+abs(E.z*E.z));

        far_field_t fields_J, fields_M;
        fields_J = config.compute_far_field_J(r, theta(i), phi , J);
        fields_M = config.compute_far_field_M(r, theta(i), phi , J);
        const real_t E_mag_J_=sqrt(abs(fields_J.E.theta*fields_J.E.theta)+abs(fields_J.E.phi*fields_J.E.phi));
        const real_t E_mag_M_=sqrt(abs(fields_M.E.theta*fields_M.E.theta)+abs(fields_M.E.phi*fields_M.E.phi));
        file_J.write("%21.14E %21.14E %21.14E\n", theta(i), 20.0*log10(E_mag_J), 20.0*log10(E_mag_J_));
        file_M.write("%21.14E %21.14E %21.14E\n", theta(i), 20.0*log10(E_mag_M), 20.0*log10(E_mag_M_)); 
    }
    timer.unset();
    file_J.close();
    file_M.close();

    theta.unset();

    quadl.unset();
    config.unset();

}


complex_t func_modal(const complex_t z, void *args_){
    function_args_t *args=(function_args_t*)args_;
    return z*z*(z-2.0)*(z-2.0)*(exp(2.*z)*cos(z)+z*z*z-1.0-sin(z)) + args->dummy;
}

void test_modal_analysis(){

    configuration_t config;
    function_args_t args;
    const complex_t j=complex_t(0.0, +1.0);
    const complex_t z=+1.6-j*3.7;
    complex_t ans;
    ans = function_derivative(func_modal, z, &args);

    contour_t contour={-3.0, -3.0, +3.0, +3.0};
    real_t k=0.0;
    args.func = func_modal;
    args.args = null;

    quadl_t quadl;
    const size_t N_quadl=32;
    const size_t k_max=15;
    const real_t tol=1.0E-6;
    quadl.set(N_quadl, k_max, tol);

    complex_t mu_k;

    mu_k = compute_mu_k(func_modal, &args, contour, k, quadl);
    // print(mu_k);

    // const size_t N=4;
    // real_t C[N]={+1.0, -2.0, +3.0, -1.0};
    // complex_t *roots=(complex_t*)calloc(N, sizeof(complex_t));
    // find_polynomial_roots(N, C, roots);
    // for (size_t i=0; i<N; i++){
    //     print(roots[i]);
    // }
    // free(roots);

    #define N_MAX 10
    complex_t *zeros=(complex_t*)calloc(N_MAX, sizeof(complex_t));
    evaluate_Delves_Lyness(func_modal, &args, contour, quadl, zeros);
    for (size_t i=0; i<N_MAX; i++){
        print(zeros[i]);
    }
    print("\n");
    for (size_t i=0; i<N_MAX; i++){
        complex_t zr=zeros[i];
        const real_t h=0.1;
        const real_t eps=1.0E-12;
        const size_t maxit=100;
        polish_Muller_methed(func_modal, &args, zr, h, eps, maxit);
        print(zr);
    }

    free(zeros);

    quadl.unset();
}


void test_CIM(){

    configuration_t config;
    function_args_t args;

    args.func = func_modal;
    args.args = null;

    quadl_t quadl;
    const size_t N_quadl=32;
    const size_t k_max=15;
    const real_t tol=1.0E-6;
    quadl.set(N_quadl, k_max, tol);

    contour_t contour={-3.0, -3.0, +3.0, +3.0};
    CIM_t CIM(contour, quadl);
    CIM.compute_zeros(func_modal, &args, func_modal, contour);
    print("\n");
    CIM.disp();

    CIM.clear();
    quadl.unset();
    
}

void test_CIM_layered_media(){

    const size_t N_layers=5;
    configuration_t config;

    const real_t lambda_0=1.3*units::um;
    const real_t freq=c_0/lambda_0;

    const complex_t j=complex_t(0.0, +1.0);
    config.set(N_layers, freq);
    size_t n=0;
    config.add_layer(n++, +0.0*units::um, +0.0*units::um, 1.0, 1.0);
    config.add_layer(n++, -0.6*units::um, +0.0*units::um, 1.0, pow(3.4-j*0.002, 2));
    config.add_layer(n++, -1.0*units::um, -0.6*units::um, 1.0, pow(3.6+j*0.010, 2));
    config.add_layer(n++, -1.6*units::um, -1.0*units::um, 1.0, pow(3.4-j*0.002, 2));
    config.add_layer(n++, -1.6*units::um, -1.6*units::um, 1.0, 1.0);

    config.log();

    termination_t top=INF, bottom=INF;
    polarization_t pol=e;
    modal_args_t args={top, bottom, pol, sheet_I, config};
    print(dispersion_function_on_sheet(1.0+j*0.01, &args));
    print("\n");
    // exit(0);

    quadl_t quadl;
    const size_t N_quadl=32;
    const size_t k_max=15;
    const real_t tol=1.0E-6;
    quadl.set(N_quadl, k_max, tol);

    contour_t contour={+0.5, -+0.0, 4.0, +0.1};
    CIM_t CIM(contour, quadl);
    CIM.compute_zeros(dispersion_function_all_sheets, &args, dispersion_function_on_sheet, contour);
    CIM.disp();

    quadl.unset();
}