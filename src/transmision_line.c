//
#include "transmission_line.h"

complex_t sqrt_Riemann(const complex_t z, const size_t sheet){
	complex_t w=csqrt(z);
	if (sheet==0){
		return cimag(w)<=0.0 ? +w : -w;
	}
	if (sheet==1){
		return cimag(w)<=0.0 ? -w : +w;
	}
	return 0.0;
}

// void spectral_components(){

// }

// complex_t gamma_TL(){

// }
