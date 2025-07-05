#ifndef __ENGINE_HPP__
#define __ENGINE_HPP__

// Libraries
#include "basic_definitions.hpp"
#include "utilities.hpp"
#include "special_functions.hpp"
#include "math_utilities.hpp"
#include "quadl.hpp"

// Definitions
class configuration_t{
    private:
        size_t N=1;
        real_t freq=0.0, omega=0.0;

    public:
        configuration_t(){}
        ~configuration_t(){}
};

// Functions

#endif