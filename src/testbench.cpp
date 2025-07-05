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
    config.unset();
}
