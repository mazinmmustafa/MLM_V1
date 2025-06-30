//
#include "quad.h"

extern void cgqf_f77_(int *rule, int *order, double *x, double *w);

void get_quad_rule(const size_t N, real_t *x, real_t *w){
	int rule=1, order=N;
	cgqf_f77_(&rule, &order, x, w);
}


complex_t quad_rule(complex_t func(complex_t, void*), const real_t a, const real_t b, void *args,
	real_t *x_quad, real_t *w_quad, const size_t N){
	const real_t h_p=(b+a)/2.0;
    const real_t h_m=(b-a)/2.0;
    complex_t sum=0.0;
    for (size_t i=0; i<N; i++){
        sum+=h_m*w_quad[i]*func(h_m*x_quad[i]+h_p, args);
    }
    return sum;
}

complex_t quadL(complex_t func(complex_t, void*), const real_t a, const real_t b, void *args,
	real_t *x_quad, real_t *w_quad, const size_t N, const real_t tol, const size_t k_max, 
	size_t *flag, const complex_t I_p, size_t k){
	complex_t I_1, I_2, I_n;
	const real_t m=(a+b)/2.0;
	I_1 = quad_rule(func, a, m, args, x_quad, w_quad, N);
	I_2 = quad_rule(func, m, b, args, x_quad, w_quad, N);
	I_n = I_1+I_2;
	const real_t error=cabs(I_n-I_p);
	(*flag) = false;
	if (k<k_max && error>tol*cabs(I_p)){
		I_1 = quadL(func, a, m, args, x_quad, w_quad, N, tol, k_max, flag, I_1, ++k);
		I_2 = quadL(func, m, b, args, x_quad, w_quad, N, tol, k_max, flag, I_2, ++k);
		I_n = I_1+I_2;
		// printf("%zu: %21.14E, %21.14E\n", k, creal(I_n), cimag(I_n));
		return I_n;
	}else{
		if (k>=k_max){(*flag) = true;}
		return I_n;
	}
}

