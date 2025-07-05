#ifndef __ENGINE_HPP__
#define __ENGINE_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"

// Definitions
template <typename type_t> 
class vector_t{
    private:
    public:
        type_t x=(type_t)0.0;
        type_t y=(type_t)0.0;
        type_t z=(type_t)0.0;
        vector_t(){}
        ~vector_t(){}
        vector_t(const type_t x, const type_t y, const type_t z){
            this->x = x;
            this->y = y;
            this->z = z;
        }
};

// Functions
void print(const vector_t<complex_t> A);
void print(const vector_t<real_t> A);

vector_t<complex_t> operator + (const vector_t<complex_t> A, const vector_t<complex_t> B);
vector_t<real_t> operator + (const vector_t<real_t> A, const vector_t<real_t> B);

vector_t<complex_t> operator - (const vector_t<complex_t> A, const vector_t<complex_t> B);
vector_t<real_t> operator - (const vector_t<real_t> A, const vector_t<real_t> B);

complex_t operator * (const vector_t<complex_t> A, const vector_t<complex_t> B);
real_t operator * (const vector_t<real_t> A, const vector_t<real_t> B);

vector_t<complex_t> operator ^ (const vector_t<complex_t> A, const vector_t<complex_t> B);
vector_t<real_t> operator ^ (const vector_t<real_t> A, const vector_t<real_t> B);

vector_t<complex_t> operator * (const complex_t a, const vector_t<complex_t> A);
vector_t<complex_t> operator * (const real_t a, const vector_t<complex_t> A);
vector_t<real_t> operator * (const real_t a, const vector_t<real_t> A);

vector_t<complex_t> operator * (const vector_t<complex_t> A, const complex_t a);
vector_t<complex_t> operator * (const vector_t<complex_t> A, const real_t a);
vector_t<real_t> operator * (const vector_t<real_t> A, const real_t a);

vector_t<complex_t> operator / (const vector_t<complex_t> A, const complex_t a);
vector_t<complex_t> operator / (const vector_t<complex_t> A, const real_t a);
vector_t<real_t> operator / (const vector_t<real_t> A, const real_t a);

real_t mag(const vector_t<real_t> A);
real_t mag(const vector_t<complex_t> A);

vector_t<real_t> unit(const vector_t<real_t> A);

int is_equal(const vector_t<real_t> A, const vector_t<real_t> B, const real_t tol);

vector_t<real_t> real_v(const vector_t<complex_t> A);
vector_t<real_t> imag_v(const vector_t<complex_t> A);

vector_t<complex_t> operator * (const vector_t<real_t> A, const complex_t a);
vector_t<complex_t> operator * (const complex_t a, const vector_t<real_t> A);
vector_t<complex_t> operator / (const vector_t<real_t> A, const complex_t a);

//

complex_t operator * (const vector_t<complex_t> A, const vector_t<real_t> B);
complex_t operator * (const vector_t<real_t> A, const vector_t<complex_t> B);
vector_t<complex_t> operator ^ (const vector_t<complex_t> A, const vector_t<real_t> B);
vector_t<complex_t> operator ^ (const vector_t<real_t> A, const vector_t<complex_t> B);

//
vector_t<real_t> get_spherical(const real_t r, const real_t theta, const real_t phi);

#endif