//
#include "engine_1d.hpp"

real_t I_1_1d(const basis_1d_t bn, const vector_t<real_t> r, const real_t a){
    projection_para_1d_t para;
    real_t ans=0.0;
    // term-
    para = projection_1d(bn.v_m, bn.e[0], r, a);
    ans += log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    // term+
    para = projection_1d(bn.e[0], bn.v_p, r, a);
    ans += log((para.P_p+para.l_p)/(para.P_m+para.l_m));
    //
    return ans;
}

vector_t<real_t> I_2_1d(const basis_1d_t bn, const vector_t<real_t> r, const real_t a){
    projection_para_1d_t para;
    vector_t<real_t> ans=vector_t<real_t>(0.0, 0.0, 0.0);
    // term-
    para = projection_1d(bn.v_m, bn.e[0], r, a);
    ans = ans+(para.P_p-para.P_m)*para.l_unit;
    // term+
    para = projection_1d(bn.e[0], bn.v_p, r, a);
    ans = ans+(para.P_p-para.P_m)*para.l_unit;
    //
    return ans;
}

vector_t<real_t> I_3_1d(const basis_1d_t bn, const vector_t<real_t> r, const real_t a){
    projection_para_1d_t para;
    vector_t<real_t> ans=vector_t<real_t>(0.0, 0.0, 0.0);
    // term-
    para = projection_1d(bn.v_m, bn.e[0], r, a);
    ans = ans+(1.0/para.P_p-1.0/para.P_m)*para.l_unit;
    ans = ans-(para.l_p/para.P_p-para.l_m/para.P_m)*para.P0_unit/para.P0;
    // term+
    para = projection_1d(bn.e[0], bn.v_p, r, a);
    ans = ans+(1.0/para.P_p-1.0/para.P_m)*para.l_unit;
    ans = ans-(para.l_p/para.P_p-para.l_m/para.P_m)*para.P0_unit/para.P0;
    //
    return ans;
}

//

struct __integrand_1d_args_t{
    basis_1d_t bm;
    basis_1d_t bn;
    real_t k;
    real_t lambda;
    quadl_domain_t quadl;
    real_t a;
    complex_t alpha;
};

complex_t psi_integrand_1_internal(const complex_t alpha_, void *args_){
    __integrand_1d_args_t *args=(__integrand_1d_args_t*)args_;
    complex_t alpha=args->alpha;
    basis_1d_t bm=args->bm, bn=args->bn;
    real_t k=args->k;
    real_t a=args->a;
    //
    const complex_t j=complex_t(0.0, 1.0);
    vector_t<real_t> R_vec;
    real_t R_vec_, R;
    complex_t A, B;
    complex_t ans=0.0;
    // term--
    R_vec = bm.v_m-bn.v_m+real(alpha)*bm.L_m[0]-real(alpha_)*bn.L_m[0];
    R_vec_ = mag(R_vec);
    R = sqrt(R_vec_*R_vec_+a*a);
    A = -j*k*exp(-j*k*R/2.0)/(4.0*pi);
    B = sinc(k*R/2.0);
    ans += (+1.0)*alpha*alpha_*A*B;
    // term-+
    R_vec = bm.v_m-bn.v_p+real(alpha)*bm.L_m[0]-real(alpha_)*bn.L_p[0];
    R_vec_ = mag(R_vec);
    R = sqrt(R_vec_*R_vec_+a*a);
    A = -j*k*exp(-j*k*R/2.0)/(4.0*pi);
    B = sinc(k*R/2.0);
    ans += (-1.0)*alpha*alpha_*A*B;
    // term+-
    R_vec = bm.v_p-bn.v_m+real(alpha)*bm.L_p[0]-real(alpha_)*bn.L_m[0];
    R_vec_ = mag(R_vec);
    R = sqrt(R_vec_*R_vec_+a*a);
    A = -j*k*exp(-j*k*R/2.0)/(4.0*pi);
    B = sinc(k*R/2.0);
    ans += (-1.0)*alpha*alpha_*A*B;
    // term++
    R_vec = bm.v_p-bn.v_p+real(alpha)*bm.L_p[0]-real(alpha_)*bn.L_p[0];
    R_vec_ = mag(R_vec);
    R = sqrt(R_vec_*R_vec_+a*a);
    A = -j*k*exp(-j*k*R/2.0)/(4.0*pi);
    B = sinc(k*R/2.0);
    ans += (+1.0)*alpha*alpha_*A*B;
    //
    return ans;
}

complex_t psi_integrand_1_external(const complex_t alpha, void *args_){
    __integrand_1d_args_t *args=(__integrand_1d_args_t*)args_;
    quadl_domain_t quadl=args->quadl;
    args->alpha = alpha;
    size_t flag;
    line_domain_t line;
    line.v1 = vector_t<real_t>(+0.0, +0.0, +0.0);
    line.v2 = vector_t<real_t>(+1.0, +0.0, +0.0);
    complex_t ans=quadl.integral_1d(psi_integrand_1_internal, args, line, flag);
    return ans;
}

complex_t Z_mn_1d(const basis_1d_t bm, const basis_1d_t bn, const real_t k, 
    const real_t lambda, quadl_domain_t quadl, const real_t a, const real_t eta){
    //
    __integrand_1d_args_t args={bm, bn, k, lambda, quadl, a, 0.0};
    complex_t I1, I2, I3, I4, I5;
    size_t flag;
    line_domain_t line;
    line.v1 = vector_t<real_t>(+0.0, +0.0, +0.0);
    line.v2 = vector_t<real_t>(+1.0, +0.0, +0.0);
    const complex_t j=complex_t(0.0, 1.0);
    //
    I1 = quadl.integral_1d(psi_integrand_1_external, &args, line, flag);
    if (flag){print("warning: no convergence!\n");}
    //
    return j*k*eta*(1.0*I1+0.0*I2+0.0*I3)-j*(eta/k)*(0.0*I4+0.0*I5);
}
