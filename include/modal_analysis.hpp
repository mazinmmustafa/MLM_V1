#ifndef __MODAL_ANALYSIS_HPP__
#define __MODAL_ANALYSIS_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "special_functions.hpp"
#include "math_utilities.hpp"
#include "quadl.hpp"
#include "vector.hpp"

// Definitions
struct contour_t{
    real_t x_1=0.0, y_1=0.0, x_2=0.0, y_2=0.0;
};

class CIM_t{
    private:
        const real_t h=0.1;
        const real_t eps=1.0E-12;
        const size_t maxit=100;
        size_t counter=0;
    public:
        complex_t *zeros=null;
        size_t N_zeros=0;
        quadl_t quadl;
        size_t N_max=4;
        contour_t contour_0;
        CIM_t(){}
        ~CIM_t(){}
        void clear(){
            free(this->zeros);
            this->N_zeros = 0;
        }
        CIM_t(const contour_t contour, quadl_t quadl){
            this->contour_0 = contour;
            this->quadl = quadl;
        }
        void compute_zeros(complex_t func(const complex_t, void*), void *args, contour_t contour);
        void display();
};

struct function_args_t{
    complex_t (*func)(const complex_t, void*);
    contour_t contour;
    real_t k=0.0;
    void *args=null;
    complex_t dummy=0.0;
};

// Functions
complex_t function_derivative(complex_t func(const complex_t, void*), const complex_t z, void *args);
complex_t compute_mu_k(complex_t func(const complex_t, void*), void *args_, 
    contour_t contour, const real_t k, quadl_t quadl);
void find_polynomial_roots(const size_t N, const real_t *C, complex_t *roots);
void evaluate_Delves_Lyness(complex_t func(const complex_t, void*), void *args, 
    contour_t contour, quadl_t quadl, complex_t *zeros);
void polish_Muller_methed(complex_t func(const complex_t, void*), void *args, 
    complex_t &xr, const real_t h, const real_t eps, const size_t maxit);

#endif