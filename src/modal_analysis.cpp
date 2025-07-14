//
#include "modal_analysis.hpp"

complex_t function_derivative(complex_t func(const complex_t, void*), const complex_t z, void *args){
    const real_t h=1.0E-6;
    const complex_t u_x=(real(func(z+h/2.0, args))-real(func(z-h/2.0, args)))/h;    
    const complex_t v_x=(imag(func(z+h/2.0, args))-imag(func(z-h/2.0, args)))/h;
    const complex_t j=complex_t(0.0, +1.0);
    return u_x+j*v_x;
}

// complex_t get_mu_k_integrand(const complex_t z, void *args_){
//     function_args_t *args=(function_args_t*)args_;
//     const complex_t z=args->z;
//     const real_t k=args->k;
//     const complex_t j=complex_t(0.0, +1.0);
//     return (function_derivative(args->func, z, args)/args->func(z, args))/(2.0*pi*j);
// }

complex_t mu_k_integrand_1(const complex_t z, void *args_){
    function_args_t *args=(function_args_t*)args_;
    const real_t k=args->k;
    const complex_t j=complex_t(0.0, +1.0);
    const real_t t=real(z);
    const real_t x=t;
    const real_t y=args->contour.y_1;
    const complex_t z_=x+j*y;
    const complex_t factor=+1.0;
    return factor*pow(z_, k)*(function_derivative(args->func, z_, args->args)/args->func(z_, args->args))/(2.0*pi*j);
}

complex_t mu_k_integrand_2(const complex_t z, void *args_){
    function_args_t *args=(function_args_t*)args_;
    const real_t k=args->k;
    const complex_t j=complex_t(0.0, +1.0);
    const real_t t=real(z);
    const real_t x=args->contour.x_2;
    const real_t y=t;
    const complex_t z_=x+j*y;
    const complex_t factor=+j;
    return factor*pow(z_, k)*(function_derivative(args->func, z_, args->args)/args->func(z_, args->args))/(2.0*pi*j);
}

complex_t mu_k_integrand_3(const complex_t z, void *args_){
    function_args_t *args=(function_args_t*)args_;
    const real_t k=args->k;
    const complex_t j=complex_t(0.0, +1.0);
    const real_t t=real(z);
    const real_t x=t;
    const real_t y=args->contour.y_2;
    const complex_t z_=x+j*y;
    const complex_t factor=-1.0;
    return factor*pow(z_, k)*(function_derivative(args->func, z_, args->args)/args->func(z_, args->args))/(2.0*pi*j);
}

complex_t mu_k_integrand_4(const complex_t z, void *args_){
    function_args_t *args=(function_args_t*)args_;
    const real_t k=args->k;
    const complex_t j=complex_t(0.0, +1.0);
    const real_t t=real(z);
    const real_t x=args->contour.x_1;
    const real_t y=t;
    const complex_t z_=x+j*y;
    const complex_t factor=-j;
    return factor*pow(z_, k)*(function_derivative(args->func, z_, args->args)/args->func(z_, args->args))/(2.0*pi*j);
}

complex_t compute_mu_k(complex_t func(const complex_t, void*), void *args_, 
    contour_t contour, const real_t k, quadl_t quadl){
    complex_t I_1, I_2, I_3, I_4;
    function_args_t args={func, contour, k, args_};
    size_t flag;
    I_1 = quadl.integral_1d(mu_k_integrand_1, &args, contour.x_1, contour.x_2, flag);
    I_2 = quadl.integral_1d(mu_k_integrand_2, &args, contour.y_1, contour.y_2, flag);
    I_3 = quadl.integral_1d(mu_k_integrand_3, &args, contour.x_1, contour.x_2, flag);
    I_4 = quadl.integral_1d(mu_k_integrand_4, &args, contour.y_1, contour.y_2, flag);
    return I_1+I_2+I_3+I_4;
}
