//
#include "testbench.h"

void test_basic_definitions(){

    printf("pi = %21.14E\n", _pi);
    printf("c_0 = %21.14E m/s\n", _c_0);
    printf("mu_0 = %21.14E H/m\n", _mu_0);
    printf("eps_0 = %21.14E F/m\n", _eps_0);
    printf("eta_0 = %21.14E ohm\n", _eta_0);
    printf("c_0 estimate = %21.14E m/s\n", 1.0/sqrt(_mu_0*_eps_0));
    printc(3.7-_1j*2.4);
}

void test_utilities(){

    const size_t Ns=1000;
    for (size_t i=0; i<Ns; i++){
        progress_bar(i, Ns, "doing something");
    }

    assert_error(1==2, "an error was captured correctly!");

}

void test_layers(){

    const size_t N=4;
    const real_t lamda=633.0*_nm;
    const real_t freq=_c_0/lamda;

    config_t config=default_config_t;
    init_config_t(&config, N, 0.0, 0.0, freq);

    size_t n=0;
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, 1.0-_1j*0.0, +500.0*_nm, +1000.0*_nm);
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, 2.0-_1j*0.0, +0.0*_nm, +500.0*_nm);
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, 10.0-_1j*0.0, -500.0*_nm, +0.0*_nm);
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, 1.0-_1j*0.0, -1000.0*_nm, -500.0*_nm);
    check_config_t(&config);

    print_config_t(&config);

    clear_config_t(&config);
}

void test_bessel(){

    const real_t range=4.0;
    const size_t Ns=100+1;
    real_t x, y;
    complex_t z, ans;

    FILE *file_x=fopen("data/tests/data_x.txt", "w"); assert(file_x!=null);
    FILE *file_y=fopen("data/tests/data_y.txt", "w"); assert(file_y!=null);
    FILE *file_m=fopen("data/tests/data_besselh_m.txt", "w"); assert(file_m!=null);
    FILE *file_p=fopen("data/tests/data_besselh_p.txt", "w"); assert(file_p!=null);
    for (size_t i=0; i<Ns; i++){
        x = linspace(-range, +range, Ns, i);
        for (size_t j=0; j<Ns; j++){    
            y = linspace(-range, +range, Ns, j);
            fprintf(file_x, "%21.14E ", x);
            fprintf(file_y, "%21.14E ", y);
            z = x+_1j*y;
            ans = besselh(2, 0.0, z);
            fprintf(file_m, "%21.14E ", cabs(ans));
            fprintf(file_p, "%21.14E ", carg(ans));
        }
        fprintf(file_x, "\n");
        fprintf(file_y, "\n");
        fprintf(file_m, "\n");
        fprintf(file_p, "\n");
    }
    fclose(file_x);
    fclose(file_y);
    fclose(file_m);
    fclose(file_p);

}

typedef struct func_arg{
    real_t a;
}func_arg;

complex_t func_test(const complex_t z, void *args_){
    func_arg *args=(func_arg*)args_;
    const real_t a=args->a;
    return ccos(a*z+_1j);
}

void test_quad(){

    const size_t N_quad=32;

    real_t *x_quad=(real_t*)calloc(N_quad, sizeof(real_t));
    real_t *w_quad=(real_t*)calloc(N_quad, sizeof(real_t));
    get_quad_rule(N_quad, x_quad, w_quad);     

    printf("Rule %zu:\n", N_quad);
    for (size_t i=0; i<N_quad; i++){
        printf("%2zu: %21.14E %21.14E\n", i, x_quad[i], w_quad[i]);
    }

    complex_t ans=0.0;

    const real_t a=-1.6, b=+2.4;
    const real_t omega=40000.0;
    func_arg args={omega};
    ans = quad_rule(func_test, a, b, &args, x_quad, w_quad, N_quad);
    printc(ans);

    const real_t tol=1.0E-8;
    const size_t k_max=30;
    size_t flag;
    printf("\n");
    ans = quadL(func_test, a, b, &args, x_quad, w_quad, N_quad, tol, k_max, &flag, 0.0, 0);
    printc(ans);
    printc(_1j*csinh(1.0-omega*_1j*b)/omega-_1j*csinh(1.0-omega*_1j*a)/omega);
    if (flag==true){printf("no convergence!\n");}
    
    free(x_quad);
    free(w_quad);

}

void test_transmission_line(){

    const size_t N_layers=4;
    const real_t lamda=633.0*_nm;
    const real_t freq=_c_0/lamda;

    config_t config=default_config_t;
    init_config_t(&config, N_layers, 0.0, 0.0, freq);

    size_t n=0;
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, 1.0-_1j*0.0, +500.0*_nm, +1000.0*_nm);
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, 2.0-_1j*0.0, +0.0*_nm, +500.0*_nm);
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, 10.0-_1j*0.0, -500.0*_nm, +0.0*_nm);
    set_layer_config_t(&config, n++, 1.0-_1j*0.0, 1.0-_1j*0.0, -1000.0*_nm, -500.0*_nm);
    check_config_t(&config);

    print_config_t(&config);

    clear_config_t(&config);

    //

    const real_t range=4.0;
    const size_t Ns=100+1;
    real_t x, y;
    complex_t z, ans;

    FILE *file_x=fopen("data/tests/data_x.txt", "w"); assert(file_x!=null);
    FILE *file_y=fopen("data/tests/data_y.txt", "w"); assert(file_y!=null);
    FILE *file_r=fopen("data/tests/data_sqrt_r.txt", "w"); assert(file_r!=null);
    FILE *file_i=fopen("data/tests/data_sqrt_i.txt", "w"); assert(file_i!=null);
    const size_t sheet=0;
    const complex_t k=2.0-_1j*0.5;
    for (size_t i=0; i<Ns; i++){
        x = linspace(-range, +range, Ns, i);
        for (size_t j=0; j<Ns; j++){    
            y = linspace(-range, +range, Ns, j);
            fprintf(file_x, "%21.14E ", x);
            fprintf(file_y, "%21.14E ", y);
            z = x+_1j*y;
            ans = sqrt_Riemann(k*k-z*z, sheet);
            fprintf(file_r, "%21.14E ", cabs(ans));
            fprintf(file_i, "%21.14E ", carg(ans));
        }
        fprintf(file_x, "\n");
        fprintf(file_y, "\n");
        fprintf(file_r, "\n");
        fprintf(file_i, "\n");
    }
    fclose(file_x);
    fclose(file_y);
    fclose(file_r);
    fclose(file_i);

}