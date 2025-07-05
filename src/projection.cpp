//
#include "projection.hpp"

void project_1d(const vector_t<real_t> v1, const vector_t<real_t> v2, const vector_t<real_t> p, 
    vector_t<real_t> &p0){
    vector_t<real_t> v21=v2-v1;
    real_t alpha=v21*(p-v1)/(v21*v21);
    p0 = v1+alpha*v21;
}

projection_para_1d_t projection_1d(const vector_t<real_t> v1, const vector_t<real_t> v2, const vector_t<real_t> p, 
    const real_t a){
    projection_para_1d_t para;
    vector_t<real_t> p0;
    project_1d(v1, v2, p, p0);
    vector_t<real_t> l_m, l_p;
    para.l_unit = unit(v2-v1);
    l_m = v1-p0;
    l_p = v2-p0;
    para.l_m = l_m*para.l_unit;
    para.l_p = l_p*para.l_unit;
    para.P0 = mag(p0-p);
    para.P0 = sqrt(para.P0*para.P0+a*a);
    para.P0_unit = unit(p0-p);
    para.P_m = sqrt(para.l_m*para.l_m+para.P0*para.P0);
    para.P_p = sqrt(para.l_p*para.l_p+para.P0*para.P0);
    return para;
}