//
#include "testbench_gold_Kretschmann.h"


void test_gold_Kretschmann_reflection_coefficients(){

    const size_t N_layers=3;
    const real_t lamda=633.0*_nm;
    const real_t freq=_c_0/lamda;

    config_t config=default_config_t;
    init_config_t(&config, N_layers, 0.0, 0.0, freq);

    size_t n=0;
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, 2.3103-_1j*0.0, +0.0*_nm, +2000.0*_nm);
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, -11.753-_1j*1.2596, -50.0*_nm, +0.0*_nm);
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, 1.0-_1j*0.0, -2000.0*_nm, -50.0*_nm);
    check_config_t(&config);
    print_config_t(&config);

    const size_t Ns=1001;
    FILE *file=fopen("data/gold_kretschmann/Refl.txt", "w");
    assert(file!=null);
    complex_t Gamma_e_u, Gamma_h_u;
    complex_t Gamma_e_d, Gamma_h_d;
    const complex_t k_1=config.layers[0].k;
    const complex_t k_N=config.layers[config.N-1].k;
    const real_t theta_min=0.0;
    const real_t theta_max=_pi/2.0;
    real_t theta;
    complex_t k_rho;
    for (size_t i=0; i<Ns; i++){
        theta = linspace(theta_min, theta_max, Ns, i);
        fprintf(file, "%21.14E ", theta*180.0/_pi);
        k_rho = k_1*sin(theta);
        Refl(k_rho, &config, 0, sheet_I, &Gamma_e_u, &Gamma_h_u, &Gamma_e_d, &Gamma_h_d);
        fprintf(file, "%21.14E %21.14E ", pow(cabs(Gamma_e_d), 2.0), pow(cabs(Gamma_h_d), 2.0));
        k_rho = k_N*sin(theta);
        Refl(k_rho, &config, config.N-1, sheet_I, &Gamma_e_u, &Gamma_h_u, &Gamma_e_d, &Gamma_h_d);
        fprintf(file, "%21.14E %21.14E ", pow(cabs(Gamma_e_u), 2.0), pow(cabs(Gamma_h_u), 2.0));
        fprintf(file, "\n");
    }
    fclose(file);

    clear_config_t(&config);
}