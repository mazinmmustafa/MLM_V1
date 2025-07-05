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

void test_mesh(){

    const real_t L=1.0*units::cm;
    const real_t port_length=0.2*units::cm;
    const real_t theta=deg_to_rad(30.0);
    const real_t phi=deg_to_rad(45.0);
    const real_t clmax=0.1*units::cm;
    const vector_t<real_t> r0=vector_t<real_t>(0.4*units::cm, 0.0*units::cm, 0.0*units::cm);

    mesh_t mesh;
    mesh.open_geometry();
    mesh.add_wire_dipole(L, -1.0*r0, +theta, +phi, clmax, 4, 1, port_length);
    mesh.add_wire_dipole(L, +0.0*r0, +theta, -phi, clmax, 5, 2, port_length);
    mesh.add_wire_dipole(L, +1.0*r0, +theta, +phi, clmax, 6, 3, port_length);
    mesh.close_geometry();
    mesh.mesh(clmax);

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

void test_shape(){

    const size_t N_ports=1;
    const real_t L=25.0*units::mm;
    const real_t port_length=2.0*units::mm;
    const real_t theta=deg_to_rad(0.0);
    const real_t phi=deg_to_rad(0.0);
    const real_t clmax=1.0*units::mm;
    const vector_t<real_t> r0=vector_t<real_t>(0.0*units::mm, 0.0*units::mm, 0.0*units::mm);
    const real_t freq=2.45*units::GHz;

    shape_t shape;
    shape.set_parameters(freq, units::mm, N_ports);
    shape.mesh.open_geometry();
    shape.mesh.add_wire_dipole(L/units::mm, r0/units::mm, 
        +theta, +phi, clmax/units::mm, 2, 1, port_length/units::mm);
    shape.mesh.close_geometry();
    shape.mesh.mesh(clmax/units::mm);
    shape.get_mesh();
    
    shape.unset();

}

struct func_I_1d_args{
    basis_1d_t bn;
    real_t a;
    vector_t<real_t> r;
    vector_t<real_t> dot;
};

complex_t func_I_1_1d(const complex_t alpha, void *args_){
    func_I_1d_args *args=(func_I_1d_args*)args_;
    complex_t ans=0.0;
    vector_t<real_t> rho;
    real_t R;
    // term-
    rho = real(alpha)*args->bn.L_m[0];
    R = mag(args->r-(args->bn.v_m+rho));
    R = sqrt(R*R+args->a*args->a);
    ans += args->bn.l_m[0]*1.0/R;
    // term+
    rho = real(alpha)*args->bn.L_p[0];
    R = mag(args->r-(args->bn.v_p+rho));
    R = sqrt(R*R+args->a*args->a);
    ans += args->bn.l_p[0]*1.0/R;
    //
    return ans;
}

complex_t func_I_2_1d(const complex_t alpha, void *args_){
    func_I_1d_args *args=(func_I_1d_args*)args_;
    complex_t ans=0.0;
    vector_t<real_t> rho;
    real_t R;
    projection_para_1d_t para;
    // term-
    rho = real(alpha)*args->bn.L_m[0];
    para = projection_1d(args->bn.v_m, args->bn.e[0], args->r, args->a);
    R = mag(args->r-(args->bn.v_m+rho));
    R = sqrt(R*R+args->a*args->a);
    ans += args->bn.l_m[0]*args->dot*(+1.0*rho+para.l_m*para.l_unit)/R;
    // term+
    rho = real(alpha)*args->bn.L_p[0];
    para = projection_1d(args->bn.e[0], args->bn.v_p, args->r, args->a);
    R = mag(args->r-(args->bn.v_p+rho));
    R = sqrt(R*R+args->a*args->a);
    ans += args->bn.l_p[0]*args->dot*(+1.0*rho+para.l_p*para.l_unit)/R;
    //
    return ans;
}

complex_t func_I_3_1d(const complex_t alpha, void *args_){
    func_I_1d_args *args=(func_I_1d_args*)args_;
    complex_t ans=0.0;
    vector_t<real_t> rho;
    real_t R;
    projection_para_1d_t para;
    // term-
    rho = real(alpha)*args->bn.L_m[0];
    para = projection_1d(args->bn.v_m, args->bn.e[0], args->r, args->a);
    R = mag(args->r-(args->bn.v_m+rho));
    R = sqrt(R*R+args->a*args->a);
    ans += -args->bn.l_m[0]*(args->dot)*(+1.0*rho+para.l_m*para.l_unit)/(R*R*R);
    ans += -args->bn.l_m[0]*(args->dot)*(para.P0*para.P0_unit)/(R*R*R);
    // term+
    rho = real(alpha)*args->bn.L_p[0];
    para = projection_1d(args->bn.e[0], args->bn.v_p, args->r, args->a);
    R = mag(args->r-(args->bn.v_p+rho));
    R = sqrt(R*R+args->a*args->a);
    ans += -args->bn.l_p[0]*(args->dot)*(+1.0*rho+para.l_p*para.l_unit)/(R*R*R);
    ans += -args->bn.l_p[0]*(args->dot)*(para.P0*para.P0_unit)/(R*R*R);
    //
    return ans;
}

void test_engine_projection_1d(){

    // test single points
    {
        quadl_domain_t quadl;
        complex_t ans;
        size_t flag;

        const real_t a=1.0E-2;
        basis_1d_t bn;
        bn.v_m = vector_t<real_t>(+1.0, +0.0, -2.0);
        bn.v_p = vector_t<real_t>(+1.0, +0.0, +1.0);
        bn.e[0] = vector_t<real_t>(+0.0, +1.0, +0.0);
        bn.get_parameters();
        const vector_t<real_t> dot=unit(vector_t<real_t>(+1.0, -1.0, +1.0));

        const vector_t<real_t> r=vector_t<real_t>(+1.0, +1.0, -2.0);

        func_I_1d_args args_1d={bn, a, r, dot};

        line_domain_t line;
        line.v1 = vector_t<real_t>(+0.0, +0.0, +0.0);
        line.v2 = vector_t<real_t>(+1.0, +0.0, +0.0);
        
        print(I_1_1d(bn, r, a));
        ans = quadl.integral_1d(func_I_1_1d, &args_1d, line, flag);
        if (flag){print("no convergence!\n");}else{print("converged!\n");}
        print(ans);

        print(I_2_1d(bn, r, a)*dot);
        ans = quadl.integral_1d(func_I_2_1d, &args_1d, line, flag);
        if (flag){print("no convergence!\n");}else{print("converged!\n");}
        print(ans);

        print(I_3_1d(bn, r, a)*dot);
        ans = quadl.integral_1d(func_I_3_1d, &args_1d, line, flag);
        if (flag){print("no convergence!\n");}else{print("converged!\n");}
        print(ans);
    }

    // test range
    // {
    //     quadl_domain_t quadl;
    //     complex_t ans;
    //     size_t flag;
    //     quadl.set_1d(15, 1.0E-4);

    //     const size_t Ns=401;
    //     const real_t R=+2.0;
    //     range_t x, y;
    //     x.set(-R, +R, Ns);
    //     y.set(-R, +R, Ns);
    //     x.linspace();
    //     y.linspace();

    //     const real_t a=1.0E-4;
    //     basis_1d_t bn;
    //     bn.v_m = vector_t<real_t>(+1.0, -1.0, +0.0);
    //     bn.v_p = vector_t<real_t>(+1.0, +1.0, +0.0);
    //     bn.e[0] = vector_t<real_t>(+0.0, +0.0, +0.0);
    //     bn.get_parameters();
    //     const vector_t<real_t> dot=unit(vector_t<real_t>(+1.0, +1.0, +0.0));

    //     vector_t<real_t> r;

    //     file_t file_x, file_y, file_data;
    //     file_x.open("data/data_x.dat", 'w');
    //     file_y.open("data/data_y.dat", 'w');
    //     for (size_t i=0; i<Ns; i++){
    //         file_x.write("%21.14E\n", x(i));
    //         file_y.write("%21.14E\n", y(i));
    //     }
    //     file_x.close();
    //     file_y.close();

    //     file_data.open("data/data.dat", 'w');

    //     line_domain_t line;
    //     line.v1 = vector_t<real_t>(+0.0, +0.0, +0.0);
    //     line.v2 = vector_t<real_t>(+1.0, +0.0, +0.0);

    //     for (size_t i=0; i<Ns; i++){
    //         for (size_t j=0; j<Ns; j++){
    //             r = vector_t<real_t>(x(i), y(j), +0.0);
    //             func_I_1d_args args_1d={bn, a, r, dot};
    //             //
    //             // complex_t ans_1=I_1_1d(bn, r, a);
    //             // complex_t ans_2=quadl.integral_1d(func_I_1_1d, &args_1d, line, flag);
    //             //
    //             // complex_t ans_1=dot*I_2_1d(bn, r, a);
    //             // complex_t ans_2=quadl.integral_1d(func_I_2_1d, &args_1d, line, flag);
    //             //
    //             complex_t ans_1=dot*I_3_1d(bn, r, a);
    //             complex_t ans_2=quadl.integral_1d(func_I_3_1d, &args_1d, line, flag);
    //             //
    //             real_t err=abs(ans_1-ans_2)/abs(ans_1);
    //             file_data.write("%21.14E ", err);
    //         }
    //         file_data.write("\n");
    //     }
        
    //     file_data.close();
    //     x.unset();
    //     y.unset();
    // }
}

void test_engine_1d(){

    quadl_domain_t quadl;
    complex_t ans;
    const real_t tol=1.0E-6;
    const size_t k_max=6;
    quadl.set_1d(k_max, tol);

    const real_t lambda=1.0;
    const real_t k=2.0*pi/lambda;
    const real_t eta_0=sqrt(mu_0/eps_0);

    const real_t a=1.0E-3;
    basis_1d_t bm, bn;
    //
    bm.v_m = vector_t<real_t>(+0.2, +0.2, +0.0);
    bm.v_p = vector_t<real_t>(+0.2, -0.2, +0.0);
    bm.e[0] = vector_t<real_t>(+0.0, +0.0, +0.0);
    //
    bn.v_m = vector_t<real_t>(+0.2, +0.2, +0.0);
    bn.v_p = vector_t<real_t>(+0.2, -0.2, +0.0);
    bn.e[0] = vector_t<real_t>(+0.0, +0.0, +0.0);
    bn.get_parameters();
    
    ans = Z_mn_1d(bm, bn, k, lambda, quadl, a, eta_0);
    print(ans);

}