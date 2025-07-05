//
#include "special_functions.hpp"

// Amos

complex_t besselj(const real_t n, const complex_t z){
	const int_t kode=1, N=1;
	int_t ierr;
	complex_t cy[16];
	xsf::amos::besj(z, n, kode, N, cy, &ierr);
	return cy[0];
}

complex_t bessely(const real_t n, const complex_t z){
	const int_t kode=1, N=1;
	int_t ierr;
	complex_t cy[16];
	xsf::amos::besy(z, n, kode, N, cy, &ierr);
	return cy[0];
}

complex_t besseli(const real_t n, const complex_t z){
	const int_t kode=1, N=1;
	int_t ierr;
	complex_t cy[16];
	xsf::amos::besi(z, n, kode, N, cy, &ierr);
	return cy[0];
}

complex_t besselk(const real_t n, const complex_t z){
	const int_t kode=1, N=1;
	int_t ierr;
	complex_t cy[16];
	xsf::amos::besk(z, n, kode, N, cy, &ierr);
	return cy[0];
}

complex_t besselh1(const real_t n, const complex_t z){
	const int_t kode=1, N=1;
	int_t ierr;
	complex_t cy[16];
	xsf::amos::besh(z, n, kode, 1, N, cy, &ierr);
	return cy[0];
}

complex_t besselh2(const real_t n, const complex_t z){
	const int_t kode=1, N=1;
	int_t ierr;
	complex_t cy[16];
	xsf::amos::besh(z, n, kode, 2, N, cy, &ierr);
	return cy[0];
}

// not propper!

complex_t sphj(const real_t n, const complex_t z){
	return sqrt(pi/(2.0*z))*besselj(n+0.5, z);
}

complex_t sphy(const real_t n, const complex_t z){
	return sqrt(pi/(2.0*z))*bessely(n+0.5, z);
}

complex_t sphh1(const real_t n, const complex_t z){
	const complex_t j=complex_t(0.0, 1.0);
	return sphj(n, z)+j*sphy(n, z);
}

complex_t sphh2(const real_t n, const complex_t z){
	const complex_t j=complex_t(0.0, 1.0);
	return sphj(n, z)-j*sphy(n, z);
}



