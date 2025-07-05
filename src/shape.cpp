//
#include "shape.hpp"

void shape_t::set_parameters(const real_t freq, const real_t metric_unit, const size_t N_ports){
    assert_error(freq>0.0, "invalid frequency");
    assert_error(metric_unit>0.0, "invalid metric unit");
    this->N_ports = N_ports;
    this->metric_unit = metric_unit;
    this->freq = freq;
    this->lambda = c_0/freq;
    this->ports_list = (port_t*)calloc(this->N_ports, sizeof(port_t));
    this->is_parameters_set = true;
}

void shape_t::get_mesh(){
    assert_error(this->is_parameters_set, "parameters are not set yet");
    file_t file;
    file.open("mesh/mesh/info.txt", 'r');
    file.read("%zu %zu %zu %zu\n", 
        &this->N_points, &this->N_lines, &this->N_triangles, &this->N_tetrahedrons);
    file.close();    
    file_t file_read, file_write;
    // basis 1d
    this->lines_list = (line_t*)calloc(this->N_lines, sizeof(line_t));
    this->is_basis_1d_list_allocated = true;
    file_read.open("mesh/mesh/elements_1d.txt", 'r');
    for (size_t i=0; i<this->N_lines; i++){
        line_t new_line;
        file_read.read("%lf %lf %lf ", &new_line.v[0].x, &new_line.v[0].y, &new_line.v[0].z);
        file_read.read("%lf %lf %lf ", &new_line.v[1].x, &new_line.v[1].y, &new_line.v[1].z);
        file_read.read("%zu\n", &new_line.pg);
        new_line.v[0] = new_line.v[0]*this->metric_unit;
        new_line.v[1] = new_line.v[1]*this->metric_unit;
        this->lines_list[i] = new_line;
    }
    file_read.close();
    file_write.open("mesh/basis/basis_1d.txt", 'w');
    line_t line_s, line_d;
    this->N_basis_1d = 0;
    for (size_t i=0; i<this->N_lines; i++){
        line_s = this->lines_list[i];
        for (size_t j=i+1; j<this->N_lines; j++){
            line_d = this->lines_list[j];
            // checking possible basis
            size_t counter=0;
            int_t indices_p[2*2];
            int_t indices_q[2*2];
            for (size_t i=0; i<2*2; i++){
                indices_p[i] = -1;
                indices_q[i] = -1;
            }
            for (size_t p=0; p<2; p++){
                for (size_t q=0; q<2; q++){
                    if (is_equal(line_s.v[p], line_d.v[q], __tol_mesh*this->lambda)){
                        indices_p[counter] = p;
                        indices_q[counter] = q;
                        counter++;
                    }
                }
            }
            assert_error(counter<2, "invalid mesh");
            if (counter==1){
                size_t check;
                check = 0;
                for (int_t k=0; k<2; k++){
                    for (int_t p=0; p<2; p++){
                        if (k!=indices_p[p]&&indices_p[p]>-1){
                            file_write.write("%21.14E %21.14E %21.14E ", line_s.v[k].x, line_s.v[k].y, line_s.v[k].z);
                            check++;
                        }
                    }
                }
                assert_error(check==1, "invalid mesh");
                check = 0;
                for (int_t k=0; k<2; k++){
                    for (int_t q=0; q<2; q++){
                        if (k!=indices_q[q]&&indices_q[q]>-1){
                            file_write.write("%21.14E %21.14E %21.14E ", line_d.v[k].x, line_d.v[k].y, line_d.v[k].z);
                            check++;
                        }
                    }
                }
                assert_error(check==1, "invalid mesh");
                check = 0;
                for (int_t k=0; k<2; k++){
                    for (int_t q=0; q<2; q++){
                        if (k==indices_q[q]&&indices_q[q]>-1){
                            file_write.write("%21.14E %21.14E %21.14E ", line_d.v[k].x, line_d.v[k].y, line_d.v[k].z);
                            check++;
                        }
                    }
                }
                assert_error(check==1, "invalid mesh");
                file_write.write("%zu %zu\n", line_s.pg, line_d.pg);
                this->N_basis_1d++;
            }
        }
    }
    file_write.close();
    this->basis_1d_list = (basis_1d_t*)calloc(this->N_basis_1d, sizeof(basis_1d_t));
    file_read.open("mesh/basis/basis_1d.txt", 'r');
    for (size_t i=0; i<this->N_basis_1d; i++){
        basis_1d_t new_basis;
        file_read.read("%lf %lf %lf ", &new_basis.v_m.x, &new_basis.v_m.y, &new_basis.v_m.z);
        file_read.read("%lf %lf %lf ", &new_basis.v_p.x, &new_basis.v_p.y, &new_basis.v_p.z);
        for (size_t j=0; j<1; j++){
            file_read.read("%lf %lf %lf ", &new_basis.e[j].x, &new_basis.e[j].y, &new_basis.e[j].z);
        }
        file_read.read("%zu %zu\n", &new_basis.pg_m, &new_basis.pg_p);
    }
    file_read.close();
    assert(this->basis_1d_list!=null);
    //
}

void shape_t::get_info_1d(size_t &N_lines, size_t &N_basis_1d){
    N_lines = this->N_lines;
    N_basis_1d = this->N_basis_1d;
}


void shape_t::assign_RF_port(const size_t port_index, const size_t port_pg, 
    const complex_t V, const complex_t Z){
    assert_error(this->is_parameters_set, "parameters are not set");
    port_t new_port={port_index, port_pg, V, Z};
    this->ports_list[port_index] = new_port;
}

void shape_t::unset(){
    if (this->is_parameters_set){
        free(this->ports_list);
        this->is_parameters_set = false;
    }
    if (this->is_basis_1d_list_allocated){
        free(this->lines_list);
        free(this->basis_1d_list);
        this->is_basis_1d_list_allocated = false;
    }
}