#ifndef __SHAPE_HPP__
#define __SHAPE_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "math_utilities.hpp"
#include "vector.hpp"
#include "mesh.hpp"

// Definitions
#define _Z0 50.0
#define __tol_mesh 1.0E-6

struct port_t{
    size_t port_index=0, port_pg=0;
    complex_t V=0.0, Z=0.0;
};

struct line_t{
    vector_t<real_t> v[2];
    size_t pg;
};

struct  basis_1d_t{
    vector_t<real_t> v_m, v_p, e[1];
    //
    vector_t<real_t> L_m[1], L_p[1];
    real_t l_m[1], l_p[1];
    size_t pg_m=0, pg_p=0;
    size_t index[2];
    const size_t N_d=1;
    void get_parameters(){
        for (size_t i=0; i<this->N_d; i++){
            this-> L_m[i] = e[i]-v_m;
            this-> L_p[i] = e[i]-v_p;
            this->l_m[i] = mag(this-> L_m[i]);
            this->l_p[i] = mag(this-> L_p[i]);
        }
    }
};

class shape_t{
    private:
        real_t freq=0.0, lambda=0.0;
        real_t metric_unit=1.0;
        size_t N_ports=0;
        port_t *ports_list=null;
        size_t N_points=0; 
        size_t N_lines=0; 
        size_t N_triangles=0; 
        size_t N_tetrahedrons=0; 
        size_t N_basis_0d=0;
        size_t N_basis_1d=0;
        size_t N_basis_2d=0;
        size_t N_basis_3d=0;
        basis_1d_t *basis_1d_list;
        line_t *lines_list=null;
        //
        size_t is_ports_list_allocated=false;
        size_t is_parameters_set=false;
        size_t is_basis_1d_list_allocated=false;
    public:
    mesh_t mesh;
    shape_t(){}
    ~shape_t(){}
    void set_parameters(const real_t freq, const real_t metric_unit, const size_t N_ports);
    void get_mesh();
    void unset();
    void get_info_1d(size_t &N_lines, size_t &N_basis_1d);
    void assign_RF_port(const size_t port_index, const size_t port_pg, 
        const complex_t V, const complex_t Z);
};

// Functions


#endif