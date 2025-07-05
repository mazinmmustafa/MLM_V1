//
#include "quadl.hpp"

size_t max_size_t(const size_t a, const size_t b){
    return a>b ? a : b;
}

size_t min_size_t(const size_t a, const size_t b){
    return a<b ? a : b;
}

extern "C" void cgqf_f77_(size_t *rule, size_t *order, double *x, double *w);

quadl_t::quadl_t(){
}

quadl_t::~quadl_t(){
}

void quadl_t::set(const size_t N, const size_t k_max, const real_t tol){
    assert_error(this->is_allocated==false, "quad is already set");
    assert_error(tol>0.0, "invalid tolerance");
    this->N = N;
    this->k_max = k_max;
    this->x = (real_t*)calloc(this->N, sizeof(real_t));
    assert(this->x!=null);
    this->w = (real_t*)calloc(this->N, sizeof(real_t));
    assert(this->w!=null);
    size_t rule=1;
    size_t N_=this->N;
    double *x_, *w_;
    x_ = (double*)calloc(this->N, sizeof(double));
    assert(this->x!=null);
    w_ = (double*)calloc(this->N, sizeof(double));
    assert(this->w!=null);
    cgqf_f77_(&rule, &N_, x_, w_);
    this->x = (real_t*)(x_);
    this->w = (real_t*)(w_);
    this->is_allocated = true;
}

void quadl_t::unset(){
    if (this->is_allocated){
        free(x);
        free(w);
        this->N = 0;
        this->k_max = 0;
        this->is_allocated = false;
    }
}

void quadl_t::disp(){
    assert_error(this->is_allocated, "quadrature rule is unset");
    for (size_t i=0; i<this->N; i++){
        print("%2zu: %21.14E, %21.14E\n", i, this->x[i], this->w[i]);
    }
}

// 1d

complex_t quadl_t::quadl_1d(complex_t (*func)(const complex_t, void*), 
    void *args, const real_t a, const real_t b){
    real_t h_p=(b+a)/2.0;
    real_t h_m=(b-a)/2.0;
    complex_t sum=0.0;
    real_t w_n;
    real_t x_n;
    for (size_t n=0; n<this->N; n++){
        w_n = this->w[n];
        x_n = this->x[n];
        sum+=h_m*w_n*func(h_m*x_n+h_p, args);
    }
    return sum;
}

complex_t quadl_t::quadl_1d_(complex_t (*func)(const complex_t, void*), 
    void *args, const real_t a, const real_t b, size_t &k, const complex_t I_p){
    real_t m=(a+b)/2.0;
    complex_t I1=quadl_t::quadl_1d(func, args, a, m);
    complex_t I2=quadl_t::quadl_1d(func, args, m, b);
    complex_t I_n=I1+I2;
    real_t error=abs(I_n-I_p);
    if (error>this->tol*abs(I_n)&&k<this->k_max&&error>0.0){
        size_t k_1, k_2;
        k_1 = k_2 = k;
        k_1++;
        k_2++;
        I1 = quadl_t::quadl_1d_(func, args, a, m, k_1, I1);
        I2 = quadl_t::quadl_1d_(func, args, m, b, k_2, I2);
        I_n = I1+I2;
        k = max_size_t(k_1, k_2);
    }
    return I_n;
}

complex_t quadl_t::integral_1d(complex_t (*func)(const complex_t, void*), 
    void *args, const real_t a, const real_t b, size_t &flag){
    assert_error(this->is_allocated, "quadrature rule is unset");
    flag = false;
    size_t k=0;
    complex_t ans=quadl_t::quadl_1d_(func, args, a, b, k, 0.0);
    if (k>=this->k_max){flag = true;}
    return ans;
}

// 2d

complex_t quadl_t::quadl_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y){
    real_t h_x_p=(b_x+a_x)/2.0;
    real_t h_x_m=(b_x-a_x)/2.0;
    real_t h_y_p=(b_y+a_y)/2.0;
    real_t h_y_m=(b_y-a_y)/2.0;
    complex_t sum=0.0;
    real_t w_m;
    real_t x_m;
    real_t w_n;
    real_t x_n;
    for (size_t m=0; m<this->N; m++){
        w_m = this->w[m];
        x_m = this->x[m];
        for (size_t n=0; n<this->N; n++){
            w_n = this->w[n];
            x_n = this->x[n];
            sum+=h_x_m*h_y_m*w_m*w_n*func(h_x_m*x_m+h_x_p, h_y_m*x_n+h_y_p, args);
        }
    }
    return sum;
}

complex_t quadl_t::quadl_2d_(complex_t (*func)(const complex_t, const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y, 
    size_t &k, const complex_t I_p){
    real_t m_x=(a_x+b_x)/2.0;
    real_t m_y=(a_y+b_y)/2.0;
    complex_t I1=quadl_t::quadl_2d(func, args, a_x, m_x, a_y, m_y);
    complex_t I2=quadl_t::quadl_2d(func, args, a_x, m_x, m_y, b_y);
    complex_t I3=quadl_t::quadl_2d(func, args, m_x, b_x, a_y, m_y);
    complex_t I4=quadl_t::quadl_2d(func, args, m_x, b_x, m_y, b_y);
    complex_t I_n=I1+I2+I3+I4;
    real_t error=abs(I_n-I_p);
    if (error>this->tol*abs(I_n)&&k<this->k_max&&error>0.0){
        size_t k_1, k_2, k_3, k_4;
        k_1 = k_2 = k_3 = k_4 = k;
        k_1++;
        k_2++;
        k_3++;
        k_4++;
        I1 = quadl_t::quadl_2d_(func, args, a_x, m_x, a_y, m_y, k_1, I1);
        I2 = quadl_t::quadl_2d_(func, args, a_x, m_x, m_y, b_y, k_2, I2);
        I3 = quadl_t::quadl_2d_(func, args, m_x, b_x, a_y, m_y, k_3, I3);
        I4 = quadl_t::quadl_2d_(func, args, m_x, b_x, m_y, b_y, k_4, I4);
        I_n = I1+I2+I3+I4; 
        k = max_size_t(k_1, max_size_t(k_2, max_size_t(k_3, k_4)));
    }
    return I_n;
}

complex_t quadl_t::integral_2d(complex_t (*func)(const complex_t, const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y,
    size_t &flag){
    assert_error(this->is_allocated, "quadrature rule is unset");
    flag = false;
    size_t k=0;
    complex_t ans=quadl_t::quadl_2d_(func, args, a_x, b_x, a_y, b_y, k, 0.0);
    if (k>=this->k_max){flag = true;}
    return ans;
}

// 3d

complex_t quadl_t::quadl_3d(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y,
    const real_t a_z, const real_t b_z){
    real_t h_x_p=(b_x+a_x)/2.0;
    real_t h_x_m=(b_x-a_x)/2.0;
    real_t h_y_p=(b_y+a_y)/2.0;
    real_t h_y_m=(b_y-a_y)/2.0;
    real_t h_z_p=(b_z+a_z)/2.0;
    real_t h_z_m=(b_z-a_z)/2.0;
    complex_t sum=0.0;
    real_t w_p;
    real_t x_p;
    real_t w_m;
    real_t x_m;
    real_t w_n;
    real_t x_n;
    for (size_t p=0; p<this->N; p++){
        w_p = this->w[p];
        x_p = this->x[p];
        for (size_t m=0; m<this->N; m++){
            w_m = this->w[m];
            x_m = this->x[m];
            for (size_t n=0; n<this->N; n++){
                w_n = this->w[n];
                x_n = this->x[n];
                sum+=h_x_m*h_y_m*h_z_m*w_p*w_m*w_n*func(
                    h_x_m*x_m+h_x_p, h_y_m*x_n+h_y_p, h_z_p*x_p+h_z_p, args);
            }
        }
    }
    return sum;
}

complex_t quadl_t::quadl_3d_(complex_t (*func)(const complex_t, const complex_t, const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y, 
    const real_t a_z, const real_t b_z, size_t &k, const complex_t I_p){
    real_t m_x=(a_x+b_x)/2.0;
    real_t m_y=(a_y+b_y)/2.0;
    real_t m_z=(a_z+b_z)/2.0;
    complex_t I1=quadl_t::quadl_3d(func, args, a_x, m_x, a_y, m_y, a_z, m_z);
    complex_t I2=quadl_t::quadl_3d(func, args, a_x, m_x, a_y, m_y, m_z, b_z);
    complex_t I3=quadl_t::quadl_3d(func, args, a_x, m_x, m_y, b_y, a_z, m_z);
    complex_t I4=quadl_t::quadl_3d(func, args, a_x, m_x, m_y, b_y, m_z, b_z);
    complex_t I5=quadl_t::quadl_3d(func, args, m_x, b_x, a_y, m_y, a_z, m_z);
    complex_t I6=quadl_t::quadl_3d(func, args, m_x, b_x, a_y, m_y, m_z, b_z);
    complex_t I7=quadl_t::quadl_3d(func, args, m_x, b_x, m_y, b_y, a_z, m_z);
    complex_t I8=quadl_t::quadl_3d(func, args, m_x, b_x, m_y, b_y, m_z, b_z);
    complex_t I_n=I1+I2+I3+I4+I5+I6+I7+I8;
    real_t error=abs(I_n-I_p);
    if (error>this->tol*abs(I_n)&&k<this->k_max&&error>0.0){
        size_t k_1, k_2, k_3, k_4, k_5, k_6, k_7, k_8;
        k_1 = k_2 = k_3 = k_4 = k_5 = k_6 = k_7 = k_8 = k;
        k_1++;
        k_2++;
        k_3++;
        k_4++;
        k_5++;
        k_6++;
        k_7++;
        k_8++;
        I1 = quadl_t::quadl_3d_(func, args, a_x, m_x, a_y, m_y, a_z, m_z, k_1, I1);
        I2 = quadl_t::quadl_3d_(func, args, a_x, m_x, a_y, m_y, m_z, b_z, k_2, I2);
        I3 = quadl_t::quadl_3d_(func, args, a_x, m_x, m_y, b_y, a_z, m_z, k_3, I3);
        I4 = quadl_t::quadl_3d_(func, args, a_x, m_x, m_y, b_y, m_z, b_z, k_4, I4);
        I5 = quadl_t::quadl_3d_(func, args, m_x, b_x, a_y, m_y, a_z, m_z, k_5, I5);
        I6 = quadl_t::quadl_3d_(func, args, m_x, b_x, a_y, m_y, m_z, b_z, k_6, I6);
        I7 = quadl_t::quadl_3d_(func, args, m_x, b_x, m_y, b_y, a_z, m_z, k_7, I7);
        I8 = quadl_t::quadl_3d_(func, args, m_x, b_x, m_y, b_y, m_z, b_z, k_8, I8);
        I_n = I1+I2+I3+I4+I5+I6+I7+I8; 
        k = max_size_t(k_1, max_size_t(k_2, max_size_t(k_3, max_size_t(k_4, 
            max_size_t(k_5, max_size_t(k_6, max_size_t(k_7, k_8)))))));
    }
    return I_n;
}

complex_t quadl_t::integral_3d(complex_t (*func)(const complex_t, const complex_t,const complex_t, void*), 
    void *args, const real_t a_x, const real_t b_x, const real_t a_y, const real_t b_y,
    const real_t a_z, const real_t b_z, size_t &flag){
    assert_error(this->is_allocated, "quadrature rule is unset");
    flag = false;
    size_t k=0;
    complex_t ans=quadl_t::quadl_3d_(func, args, a_x, b_x, a_y, b_y, a_z, b_z, k, 0.0);
    if (k>=this->k_max){flag = true;}
    return ans;
}
