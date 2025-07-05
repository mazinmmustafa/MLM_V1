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
    real_t a=0.0, b=0.0, c=0.0;
};

complex_t func_1d(const complex_t x, void *args_){
    func_args *args=(func_args*)args_;
    return args->a+args->b*x-args->c*x*x+exp(-x);
    // return args->a*log(1.0+x);
}

complex_t func_2d(const complex_t x, const complex_t y, void *args_){
    func_args *args=(func_args*)args_;
    return args->a+args->b*(x+y)-args->c*x*y+args->c*x*x+exp(-x-y);
    // return args->a*log(1.0+x+y);
}

complex_t func_3d(const complex_t x, const complex_t y, const complex_t z, void *args_){
    func_args *args=(func_args*)args_;
    return args->a+args->b*(x+y)-args->c*x*y+args->c*z*z+exp(-x-y-z);
    // return args->a*log(2.0+x+y+z);
}

void test_quadl(){

    func_args args={+1.0, +1.0, +1.0};

    vector_t<real_t> v1, v2, v3, v4;
    v1 = vector_t<real_t>(+0.0, +0.0, +0.0);
    v2 = vector_t<real_t>(+1.0, +0.0, +0.0);
    v3 = vector_t<real_t>(+0.0, +1.0, +0.0);
    v4 = vector_t<real_t>(+0.0, +0.0, +1.0);

    quadl_domain_t quadl;
    const real_t tol=1.0E-6;
    quadl.set_1d(20, tol);
    quadl.set_2d(20, tol);
    quadl.set_3d(20, tol);

    line_domain_t line;
    line.v1 = v1;
    line.v2 = v2;
    
    triangle_domain_t triangle;
    triangle.v1 = v1;
    triangle.v2 = v2;
    triangle.v3 = v3;

    tetrahedron_domain_t tetrahedron;
    tetrahedron.v1 = v1;
    tetrahedron.v2 = v2;
    tetrahedron.v3 = v3;
    tetrahedron.v4 = v4;

    complex_t ans=0.0;
    size_t flag;

    // +3.86294E-01
    ans = quadl.integral_1d(func_1d, &args, line, flag);
    if (flag){print("no convergence!\n");}else{print("converged!\n");}
    print(ans);

    // +2.50000E-01
    ans = quadl.integral_2d(func_2d, &args, triangle, flag);
    if (flag){print("no convergence!\n");}else{print("converged!\n");}
    print(ans);

    // +1.68167E-01
    ans = quadl.integral_3d(func_3d, &args, tetrahedron, flag);
    if (flag){print("no convergence!\n");}else{print("converged!\n");}
    print(ans);

}
