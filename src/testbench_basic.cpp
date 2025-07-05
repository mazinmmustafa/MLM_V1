//
#include "testbench_basic.hpp"

void test_bessel(){

    const size_t Ns=401;
    const real_t R=2.0;

    range_t x, y;
    x.set(-R, +R, Ns);
    y.set(-R, +R, Ns);

    x.linspace();
    y.linspace();

    file_t file_x, file_y, file_data;
    file_x.open("data/x_bessel.dat", 'w');
    file_y.open("data/y_bessel.dat", 'w');
    file_data.open("data/data_bessel.dat", 'w');

    complex_t z;
    const complex_t j=complex_t(0.0, 1.0);

    for (size_t m=0; m<Ns; m++){
        file_x.write("%21.14E\n", x(m));
        file_y.write("%21.14E\n", y(m));
        for (size_t n=0; n<Ns; n++){
            z = x(m)+j*y(n);
            file_data.write("%21.14E ", abs(besselh2(0, z)));
        }
        file_data.write("\n");
    }
    
    file_data.close();
    file_x.close();
    file_y.close();

    x.unset();
    y.unset();

}

void test_matrix(){

    size_t N=3;
    matrix_t<complex_t> A, b, x;
    A.set(N, N);
    b.set(N, 1);
    x.set(N, 1);
    //
    A(0, 0) = +1.0;
    A(0, 1) = +2.0;
    A(0, 2) = +3.0;

    A(1, 0) = +3.0;
    A(1, 1) = +2.0;
    A(1, 2) = +1.0;

    A(2, 0) = +1.0;
    A(2, 1) = +3.0;
    A(2, 2) = +4.0;

    A.lup();

    b(0, 0) = +1.0;
    b(1, 0) = +2.0;
    b(2, 0) = +3.0;

    A.solve(b, x);

    x.disp(N);

    A.unset();
    b.unset();
    x.unset();

}

struct func_args{
    real_t a=0.0;
};

complex_t func_1d(const complex_t x, void *args_){
    func_args *args=(func_args*)args_;
    const complex_t j=complex_t(0.0, +1.0);
    return cos(args->a*x+j);
}

void test_quadl(){

    quadl_t quadl;
    const size_t N_quadl=32;
    const size_t k_max=15;

    const real_t Omega=1000;
    func_args args={Omega};
    const real_t a=-1.6, b=+3.2;

    const real_t tol=1.0E-6;
    quadl.set(N_quadl, k_max, tol);

    complex_t ans=0.0;
    size_t flag;

    ans = quadl.integral_1d(func_1d, &args, a, b, flag);
    if (flag){print("no convergence!\n");}else{print("converged!\n");}
    print(ans);
    const complex_t j=complex_t(0.0, +1.0);
    print("correct answer:\n");
    print((sin(Omega*b+j)-sin(Omega*a+j))/Omega);

}
