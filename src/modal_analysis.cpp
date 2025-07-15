//
#include "modal_analysis.hpp"

extern "C" void n_eigen(double *_a, int n, double *wr, double *wi);

complex_t function_derivative(complex_t func(const complex_t, void*), const complex_t z, void *args){
    const real_t h=1.0E-6;
    const complex_t u_x=(real(func(z+h/2.0, args))-real(func(z-h/2.0, args)))/h;    
    const complex_t v_x=(imag(func(z+h/2.0, args))-imag(func(z-h/2.0, args)))/h;
    const complex_t j=complex_t(0.0, +1.0);
    return u_x+j*v_x;
}

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

void find_polynomial_roots(const size_t N, const real_t *C, complex_t *roots){
    real_t *a=(real_t*)calloc(N*N, sizeof(real_t));
    real_t *wr=(real_t*)calloc(N, sizeof(real_t));
    real_t *wi=(real_t*)calloc(N, sizeof(real_t));
    size_t counter=0;
    for (size_t i=0; i<(N*N); i++){
        a[i] = 0.0;
    }
    for (size_t i=0; i<N; i++){
        for (size_t j=0; j<N; j++){
            if (i-1==j){
                a[counter] = 1.0;  
            }
            if (j==N-1){
                a[counter] = -C[i]; 
            }
            counter++;
        }
    }
    n_eigen(a, N, wr, wi);
    const complex_t j=complex_t(0.0, +1.0);
    for (size_t i=0; i<N; i++){
        roots[i] = wr[i]+j*wi[i];
    }
    free(a);
    free(wr);
    free(wi);
}

#define __min_tol 1.0E-2

void evaluate_Delves_Lyness(complex_t func(const complex_t, void*), void *args, 
    contour_t contour, quadl_t quadl, complex_t *zeros){
    complex_t N_estimate=compute_mu_k(func, args, contour, 0.0, quadl);
    if (abs(imag(N_estimate))<__min_tol && real(N_estimate)>=0.0){
        const size_t N=real(N_estimate);
        real_t *mu_k=(real_t*)calloc(N+1, sizeof(real_t));
        for (size_t k=0; k<=N; k++){
            mu_k[k] = real(compute_mu_k(func, args, contour, k, quadl));
        }
        real_t *sigma=(real_t*)calloc(N+1, sizeof(real_t));
        sigma[0] = 1.0;
        for (size_t k=1; k<=N; k++){
            real_t sum=0.0;
            for (size_t j=1; j<=(k-1); j++){
                sum = sum+sigma[k-j]*mu_k[j];
            }
            sigma[k] = (-1.0/k)*(mu_k[k]+sum);
        }
        real_t *C=(real_t*)calloc(N, sizeof(real_t));
        for (size_t i=0; i<N; i++){
            C[N-i-1] = sigma[i+1];
        }
        find_polynomial_roots(N, C, zeros);
        free(mu_k);
        free(sigma);
        free(C);
    }
}

void polish_Muller_methed(complex_t func(const complex_t, void*), void *args, 
    complex_t &xr, const real_t h, const real_t eps, const size_t maxit){
    complex_t x2=xr;
    complex_t x1=xr+h*xr;
    complex_t x0=xr-h*xr; 
    for (size_t iter=0; iter<maxit; iter++){
        complex_t h0=x1-x0;
        complex_t h1=x2-x1;
        complex_t d0=(func(x1, args)-func(x0, args))/h0;
        complex_t d1=(func(x2, args)-func(x1, args))/h1;
        complex_t a=(d1-d0)/(h1+h0);
        complex_t b=a*h1+d1;
        complex_t c=func(x2, args);
        complex_t rad=sqrt(b*b-4.0*a*c);
        complex_t den;
        if (abs(b+rad)>abs(b-rad)){
            den = b+rad;
        }else{
            den = b-rad;
        }
        complex_t dxr=-2.0*c/den;
        xr = x2+dxr;
        if (abs(dxr)<eps*abs(xr)){
            return;
        }else{
            x0 = x1;
            x1 = x2;
            x2 = xr;
        }
    }
}

//

void CIM_t::compute_zeros_internal(complex_t func(const complex_t, void*), void *args, contour_t contour){
    const size_t trial_max=10;
    complex_t N_estimate=0.0;
    size_t is_box_ok=false;
    if (this->counter==0){
        contour = this->contour_0;
        this->counter++;
    }
    for (size_t i=0; i<trial_max; i++){
        N_estimate = compute_mu_k(func, args, contour, 0.0, this->quadl);
        if (abs(imag(N_estimate))<__min_tol && real(N_estimate)>=0.0){
            is_box_ok = true;
            break;
        }else{
            contour.x_1 -= __min_tol;
            contour.x_2 += __min_tol;
            contour.y_1 -= __min_tol;
            contour.y_2 += __min_tol;
        }
    }
    assert_error(is_box_ok==true, "unable to find a usable box");
    const size_t N=real(N_estimate);
    if (N>this->N_max){
        const real_t x_m=(contour.x_1+contour.x_2)/2.0;
        const real_t y_m=(contour.y_1+contour.y_2)/2.0;
        contour_t new_contour_1={contour.x_1, contour.y_1, x_m, y_m};
        contour_t new_contour_2={x_m, contour.y_1, contour.x_2, y_m};
        contour_t new_contour_3={contour.x_1, y_m, x_m, contour.y_2};
        contour_t new_contour_4={x_m, y_m, contour.x_2, contour.y_2};
        this->compute_zeros_internal(func, args, new_contour_1);
        this->compute_zeros_internal(func, args, new_contour_2);
        this->compute_zeros_internal(func, args, new_contour_3);
        this->compute_zeros_internal(func, args, new_contour_4);
        
    }
    this->N_zeros_temp+=N;
    complex_t *new_zeros=(complex_t*)calloc(N, sizeof(complex_t));
    evaluate_Delves_Lyness(func, args, contour, this->quadl, new_zeros);
    for (size_t i=0; i<N; i++){
        polish_Muller_methed(func, args, new_zeros[i], this->h, this->eps, this->maxit);
        this->zeros_temp.push(new_zeros[i]);
    }
    free(new_zeros);
}

complex_t round_n(const complex_t z, const real_t eps){
    real_t x=real(z), y=imag(z);
    const complex_t j=complex_t(0.0, +1.0);
    const int_t n=log10(eps);
    const real_t base=pow(10.0, n);
    x = base*round(x/base);
    y = base*round(y/base);
    return x+j*y;
}

void CIM_t::compute_zeros(complex_t func(const complex_t, void*), void *args, contour_t contour){
    this->compute_zeros_internal(func, args, contour);
    struct data_temp_t{
        complex_t data=0;
        size_t flag=false;
    };
    data_temp_t *data_temp=(data_temp_t*)calloc(this->N_zeros_temp, sizeof(data_temp_t));

    for (size_t i=0; i<this->N_zeros_temp; i++){
        data_temp[i].data = round_n(this->zeros_temp.top(), this->eps*10.0);
        this->zeros_temp.pop();
    }
    for (size_t i=0; i<this->N_zeros_temp; i++){
        for (size_t k=i+1; k<this->N_zeros_temp; k++){
            if (abs(data_temp[i].data-data_temp[k].data)<this->eps){
                data_temp[k].flag = true;
            }
        }
    }
    for (size_t i=0; i<this->N_zeros_temp; i++){
        if (data_temp[i].flag==false){
            this->zeros_temp.push(data_temp[i].data);
            this->N_zeros++;
        }
    }
    free(data_temp);
    this->zeros = (complex_t*)calloc(this->N_zeros, sizeof(complex_t));
    assert(this->zeros!=null);
    this->is_allocated = true;
    for (size_t i=0; i<this->N_zeros; i++){
        this->zeros[i] = this->zeros_temp.top();
        this->zeros_temp.pop();
    }
}

void CIM_t::display(){
    for (size_t i=0; i<this->N_zeros; i++){
        print(this->zeros[i]);
    }
}