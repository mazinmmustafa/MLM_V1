//
#include "bessel.h"

// Fortran Functions
extern void zbesj_f77_(double *fnu, double *zr, 
	double *zi, double *cyr, double *cyi);
	
extern void zbesy_f77_(double *fnu, double *zr, 
	double *zi, double *cyr, double *cyi);
	
extern void zbesi_f77_(double *fnu, double *zr, 
	double *zi, double *cyr, double *cyi);
	
extern void zbesk_f77_(double *fnu, double *zr, 
	double *zi, double *cyr, double *cyi);
	
extern void zbesh_f77_(int *k, double *fnu, double *zr, 
	double *zi, double *cyr, double *cyi);
	
complex_t besselj(real_t n, complex_t z){
	real_t zr=creal(z);
	real_t zi=cimag(z);
	real_t cyr, cyi;
	complex_t j=csqrt(-1.0);
	zbesj_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex_t bessely(real_t n, complex_t z){
	real_t zr=creal(z);
	real_t zi=cimag(z);
	real_t cyr, cyi;
	complex_t j=csqrt(-1.0);
	zbesy_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex_t besseli(real_t n, complex_t z){
	real_t zr=creal(z);
	real_t zi=cimag(z);
	real_t cyr, cyi;
	complex_t j=csqrt(-1.0);
	zbesi_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex_t besselk(real_t n, complex_t z){
	real_t zr=creal(z);
	real_t zi=cimag(z);
	real_t cyr, cyi;
	complex_t j=csqrt(-1.0);
	zbesk_f77_(&n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}

complex_t besselh(int_t k, real_t n, complex_t z){
	assert(k==1||k==2);
	real_t zr=creal(z);
	real_t zi=cimag(z);
	real_t cyr, cyi;
	complex_t j=csqrt(-1.0);
	zbesh_f77_(&k, &n, &zr, &zi, &cyr, &cyi);
	return cyr+j*cyi;
}