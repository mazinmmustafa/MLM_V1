//
#include "mesh.hpp"

void mesh_t::open_geometry(){
    assert_error(this->is_file_open==false, "file is already open");
    this->file.open(__gmsh_filename, 'w');
    this->file.write("SetFactory(\"OpenCASCADE\");\n");
    this->file.close();
    this->is_file_open = true;
    this->point_index = 0;
    this->line_index = 0;
    this->curve_index = 0;
}

void mesh_t::close_geometry(){
    assert_error(this->is_file_open==true, "file is not open yet");
    this->is_file_open = false;
    this->point_index = 0;
    this->line_index = 0;
    this->curve_index = 0;
}

//

void mesh_t::add_wire_dipole(const real_t length, const vector_t<real_t> center, 
    const real_t theta, const real_t phi, const real_t clmax, 
    const size_t wire_pg, const size_t port_pg, const real_t port_length){
    assert_error(length>0, "invalid length");
    assert_error(clmax>0, "invalid clmax");
    assert_error(port_length<=length && port_length>0, "invalid port length");
    vector_t<real_t> v1, v2, v3, v4, v5;
    v1 = center-get_spherical(length/2.0, theta, phi);
    v2 = center-get_spherical(port_length/2.0, theta, phi);
    v3 = center;
    v4 = center+get_spherical(port_length/2.0, theta, phi);
    v5 = center+get_spherical(length/2.0, theta, phi);
    file_t file;
    file.open(__gmsh_filename, 'a');
    file.write("//\n");
    file.write("Point(%zu) = {%21.14E, %21.14E, %21.14E};\n", ++this->point_index, v1.x, v1.y, v1.z);
    file.write("Point(%zu) = {%21.14E, %21.14E, %21.14E};\n", ++this->point_index, v2.x, v2.y, v2.z);
    file.write("Point(%zu) = {%21.14E, %21.14E, %21.14E};\n", ++this->point_index, v3.x, v3.y, v3.z);
    file.write("Point(%zu) = {%21.14E, %21.14E, %21.14E};\n", ++this->point_index, v4.x, v4.y, v4.z);
    file.write("Point(%zu) = {%21.14E, %21.14E, %21.14E};\n", ++this->point_index, v5.x, v5.y, v5.z);
    file.write("MeshSize {%zu, %zu, %zu, %zu, %zu} = %21.14E;\n", 
    this->point_index-4, this->point_index-3, this->point_index-2, this->point_index-1, this->point_index-0, 
    clmax);
    file.write("Line(%zu) = {%zu, %zu};\n", ++this->line_index, this->point_index-4, this->point_index-3);
    file.write("Line(%zu) = {%zu, %zu};\n", ++this->line_index, this->point_index-3, this->point_index-2);
    file.write("Line(%zu) = {%zu, %zu};\n", ++this->line_index, this->point_index-2, this->point_index-1);
    file.write("Line(%zu) = {%zu, %zu};\n", ++this->line_index, this->point_index-1, this->point_index-0);
    file.write("Physical Curve(\"Wire%zu\", %zu) = {%zu, %zu};\n", 
        ++this->curve_index, wire_pg, this->line_index-3, this->line_index-0);
    file.write("Physical Curve(\"Port%zu\", %zu) = {%zu, %zu};\n", 
        ++this->curve_index, port_pg, this->line_index-2, this->line_index-1);
    file.close();
}

void mesh_t::mesh(const real_t clmax){
    //
    assert_error(clmax>0.0, "invalid mesh tolerance");
    int_t max_length=200;
    char *cmd=(char*)calloc(max_length, sizeof(char));
    print("calling gmsh...");
    sprintf(cmd, "gmsh %s -3 -clmax %0.4f -format vtk -save_all -o mesh/shape.vtk > mesh/shape_log.txt", __gmsh_filename, clmax);
    assert_error(!system(cmd), "unable to mesh geometry");
    print(", done!\n");
    #ifdef __windows__
    sprintf(cmd, "python mesh/read_vtk.py");
    #endif
    #ifdef __linux__
    sprintf(cmd, "python3 mesh/read_vtk.py");
    #endif
    assert_error(!system(cmd), "unable to generate mesh");
    free(cmd);
    //
    
}
