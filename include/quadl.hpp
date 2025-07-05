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

struct line_domain_t{ // domain x
    vector_t<real_t> v1, v2;
    real_t length(){
        return abs(v2.x-v1.x);
    }
};

struct triangle_domain_t{ // domain x, y
    vector_t<real_t> v1, v2, v3;
    real_t area(){
        return mag((v2-v1)^(v3-v1))/2.0;
    }
};

struct tetrahedron_domain_t{ // domain x, y, z
    vector_t<real_t> v1, v2, v3, v4;
    real_t volume(){
        return ((v2-v1)^(v3-v1))*(v4-v1)/6.0;
    }
};
 	
class quadl_domain_t{
    private:
        // 1d
        size_t k_max_1d=40;
        real_t tol_1d=1.0E-4;
        const size_t N_1d=3;
        const real_t x_1d[3]={(1.0-sqrt(0.6))/2.0, (1.0+0.0)/2.0, (1.0+sqrt(0.6))/2.0};
        const real_t w_1d[3]={5.0/18.0, 4.0/9.0, 5.0/18.0};
        // 2d
        size_t k_max_2d=20;
        real_t tol_2d=1.0E-4;
        const size_t N_2d=4;
        const real_t x_2d[4]={1.0/3.0, 0.6, 0.2, 0.2};
        const real_t y_2d[4]={1.0/3.0, 0.2, 0.6, 0.2};
        const real_t w_2d[4]={-27.0/48.0, 25.0/48.0, 25.0/48.0, 25.0/48.0};
        // 3d
        size_t k_max_3d=10;
        real_t tol_3d=1.0E-4;
        const size_t N_3d=5;
        const real_t x_3d[5]={1.0/4.0, 1.0/2.0, 1.0/6.0, 1.0/6.0, 1.0/6.0};
        const real_t y_3d[5]={1.0/4.0, 1.0/6.0, 1.0/2.0, 1.0/6.0, 1.0/6.0};
        const real_t z_3d[5]={1.0/4.0, 1.0/6.0, 1.0/6.0, 1.0/2.0, 1.0/6.0};
        const real_t w_3d[5]={-4.0/5.0, 9.0/20.0, 9.0/20.0, 9.0/20.0, 9.0/20.0};
        //
        complex_t quadl_1d(complex_t (*func)(const complex_t, void*), 
            void *args, const line_domain_t line);
        complex_t quadl_1d_(complex_t (*func)(const complex_t, void*), 
            void *args, line_domain_t line, size_t &k, const complex_t I_p);
        complex_t quadl_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, triangle_domain_t triangle);
        complex_t quadl_2d_(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, const triangle_domain_t triangle, size_t &k, const complex_t I_p);
        complex_t quadl_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, tetrahedron_domain_t tetrahedron);
        complex_t quadl_3d_(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, const tetrahedron_domain_t tetrahedron, size_t &k, const complex_t I_p);
    public:
        quadl_domain_t(){}
        ~quadl_domain_t(){}
        void set_1d(const size_t k_max, const real_t tol);
        void set_2d(const size_t k_max, const real_t tol);
        void set_3d(const size_t k_max, const real_t tol);
        complex_t integral_1d(complex_t (*func)(const complex_t, void*), 
            void *args, const line_domain_t line, size_t &flag);
        complex_t integral_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
            void *args, const triangle_domain_t triangle, size_t &flag);
        complex_t integral_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
            void *args, const tetrahedron_domain_t tetrahedron, size_t &flag);
};

#endif
