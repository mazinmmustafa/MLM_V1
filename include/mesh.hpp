#ifndef __MESH_HPP__
#define __MESH_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "math_utilities.hpp"
#include "vector.hpp"

// Definitions
#define __gmsh_filename "mesh/shape.geo"

class mesh_t{
    private:
        file_t file;
        size_t is_file_open=false;
        size_t point_index, line_index;
        size_t curve_index;
    public:
        mesh_t(){}
        ~mesh_t(){}
        void open_geometry();
        void close_geometry();
        void mesh(const real_t clmax);
        // shapes
        void add_wire_dipole(const real_t length, const vector_t<real_t> center, 
            const real_t theta, const real_t phi, const real_t clmax, 
            const size_t wire_pg, const size_t port_pg, const real_t port_length);
};

// Functions


#endif