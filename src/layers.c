//
#include "layers.h"

void init_layer_t(layer_t *layer, const complex_t mu, const complex_t eps,
    const real_t z_d, const real_t z_u){
    //
    assert(layer!=null);
    layer->mu = mu;
    layer->eps = eps;
    layer->k = csqrt(mu*eps);
    layer->lambda = 1.0/layer->k;
    layer->eta = csqrt(mu/eps);
    layer->z_u = z_u;
    layer->z_d = z_d;
    layer->d = z_u-z_d;
    assert_error(z_u>=z_d, "invalid layer description!");
}


void init_config_t(config_t *config, const size_t N, const complex_t Gamma_u, const complex_t Gamma_d, 
    const real_t freq){
    assert(config!=null);
    clear_config_t(config);
    assert_error(N>0, "invalid number of layers");
    assert_error(freq>0, "invalid operating frequency");
    config->N = N;
    config->Gamma_d = Gamma_d;
    config->Gamma_u = Gamma_u;
    config->freq = freq;
    config->omega = 2.0*_pi*freq;
    config->lambda_0 = _c_0/freq;
    config->k_0 = 2.0*_pi/config->lambda_0;
    config->layers = (layer_t*)calloc(N, sizeof(layer_t));
    assert(config->layers!=null);
    config->is_allocated = true;
}

void clear_config_t(config_t *config){
    assert(config!=null);
    config->freq = 0.0;
    config->omega = 0.0;
    config->Gamma_d = 0.0;
    config->Gamma_u = 0.0;
    config->k_0 = 0.0;
    config->lambda_0 = 0.0;
    config->N = 0;
    if (config->is_allocated==true){
        free(config->layers);
        config->is_allocated = false;
    }
    config->is_allocated = false;
}

void set_layer_config_t(config_t *config, const size_t n, const complex_t mu, const complex_t eps,
    const real_t z_d, const real_t z_u){
    assert(config!=null);
    assert_error(config->is_allocated==true, "configuration is not yet allocated");
    assert_error(n<config->N, "layer index is out of range");
    layer_t layer;
    init_layer_t(&layer, mu, eps, z_d, z_u);
    config->layers[n] = layer;
}

void check_config_t(config_t *config){
    assert(config!=null);
    assert_error(config->is_allocated==true, "configuration is not yet allocated");
    for (size_t n=0; n<config->N-1; n++){
        char msg[100];
        sprintf(msg, "invalid layer %zu description", n);
        assert_error(config->layers[n].z_d==config->layers[n+1].z_u, msg);
    }
    for (size_t n=0; n<config->N; n++){
        config->layers[n].lambda*=config->lambda_0;
        config->layers[n].k*=config->k_0;
        config->layers[n].eta*=_eta_0;
    }
}

void print_config_t(config_t *config){
    assert(config!=null);
    assert_error(config->is_allocated==true, "configuration is not yet allocated");
    FILE *file=fopen("data/configuration.txt", "w");
    fprintf(file, "configuration details:\n\n");
    fprintf(file, "number of layers is %zu\n", config->N);
    fprintf(file, "frequency is %21.14E Hz\n", config->freq);
    fprintf(file, "radial frequency is %21.14E rad/s\n", config->omega);
    fprintf(file, "free-space wavelength is %21.14E m\n", config->lambda_0);
    fprintf(file, "free-space wave constant is %21.14E rad/m\n", config->k_0);
    fprintf(file, "\nlayers details:\n");
    for (size_t n=0; n<config->N; n++){
        fprintf(file, "\nlayer %zu:\n", n);
        fprintf(file, "range from %21.14E m to %21.14E m\n", config->layers[n].z_d, config->layers[n].z_u);
        fprintf(file, "thickness is %21.14E m\n", config->layers[n].d);
        fprintf(file, "magnetic constant is %21.14E +j %21.14E\n", creal(config->layers[n].mu), cimag(config->layers[n].mu));
        fprintf(file, "electric constant is %21.14E +j %21.14E\n", creal(config->layers[n].eps), cimag(config->layers[n].eps));
        fprintf(file, "intrinsic impedance is %21.14E +j %21.14E ohm\n", creal(config->layers[n].eta), cimag(config->layers[n].eta));
        fprintf(file, "complex wavelength is %21.14E +j %21.14E m\n", creal(config->layers[n].lambda), cimag(config->layers[n].lambda));
        fprintf(file, "complex wave constant is %21.14E +j %21.14E rad/m\n", creal(config->layers[n].k), cimag(config->layers[n].k));
    }
    fclose(file);
}
