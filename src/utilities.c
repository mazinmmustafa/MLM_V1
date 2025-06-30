//
#include "utilities.h"

void __assert_error(const size_t condition, const char *msg, const char *file, const size_t line){
    if (condition==false){
        printf("%s:%zu: error: %s\n", file, line, msg);
        exit(0);
    }
}

void progress_bar(const size_t n, const size_t Ns, const char *msg){
    assert(msg!=null);
    const size_t N_bars=10;
    printf("\rprogress: [");
    for (size_t i=0; i<(n+1)*N_bars/Ns; i++){
        printf("#");
    }
    for (size_t i=(n+1)*N_bars/Ns; i<N_bars; i++){
        printf("-");
    }
    printf("] %3zu%%, %s", 100*(n+1)/Ns, msg);
    if ((n+1)==Ns){
        printf(" ,done!\n");
    }
}

void printc(const complex_t z){
    printf("(%21.14E, %21.14E)\n", creal(z), cimag(z));
}

real_t linspace(const real_t x_min, const real_t x_max, const size_t Ns, const size_t i){
    const real_t dx=(x_max-x_min)/(Ns-1);
    return x_min+i*dx;
}

real_t logspace(const real_t x_min, const real_t x_max, const size_t Ns, const size_t i){
    const real_t X_min=log10(x_min);
    const real_t X_max=log10(x_max);
    const real_t dX=(X_max-X_min)/(Ns-1);
    return pow(10.0, X_min+i*dX);
}