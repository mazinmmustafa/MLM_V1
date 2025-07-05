#ifndef __QUADL_HPP__
#define __QUADL_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "vector.hpp"

// Definitions
class quadl_t{
    private:
        size_t N=0;
        size_t k_max=0;
        real_t *x=null, *w=null;
        size_t is_allocated=false;
        real_t tol=1.0E-4;
        //
        complex_t quadl_1d(complex_t (*func)(const complex_t, void*), 
            void *args, const real_t a, const real_t b);
        complex_t quadl_1d_(complex_t (*func)(const complex_t, void*), 
            void *args, const real_t a, const real_t b, size_t &k, const complex_t I_p);
        complex_t quadl_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y);
        complex_t quadl_2d_(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y, 
            size_t &k, const complex_t I_p);
        complex_t quadl_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y,
            const real_t a_z, const real_t b_z);
        complex_t quadl_3d_(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y, 
            const real_t a_z, const real_t b_z, size_t &k, const complex_t I_p);
    public:
        quadl_t();
        ~quadl_t();
        void set(const size_t N, const size_t k_max, const real_t tol);
        void unset();
        void disp();
        complex_t integral_1d(complex_t (*func)(const complex_t, void*), 
            void *args, const real_t a, const real_t b, size_t &flag);
        complex_t integral_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y, size_t &flag);
        complex_t integral_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y,
            const real_t a_z, const real_t b_z, size_t &flag);
};


#endif
