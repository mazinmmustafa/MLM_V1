#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__

// Libraries
#include "basic_definitions.hpp"

// Definitions
#define assert_error(check, msg) __assert_error((check), (msg), __FILE__, __LINE__);

class range_t{
    private:
        real_t *data=null;
        real_t x_min=0.0, x_max=0.0;
        size_t Ns=0;
        int is_allocated=false;
    public:
        range_t();
        ~range_t();
        void linspace();
        void logspace();
        real_t operator() (const size_t index) const;
        void set(const real_t x_min, const real_t x_max, const size_t Ns);
        void unset();
        void get_info(real_t &x_min, real_t &x_max, size_t &Ns);
};

class stopwatch_t{
private:
    int is_set=false;
    #ifdef __windows__
    time_t start, stop;
    #endif
    #ifdef __linux__
    struct timespec start, stop;
    #endif
    double elapsed=0.0;
public:
    stopwatch_t();
    ~stopwatch_t();
    void set();
    void unset();
    void unset_silent();
    double get_elapsed();
};

class file_t{
    protected:
        FILE *file_ptr=NULL;
        char mode='0';
        const size_t max_length=200;
        char *filename=NULL;
        int is_open=false;
    public:
        file_t();
        ~file_t();
        void open(const char *filename, const char mode);
        void close();
        void write(const char *format, ...);
        int_t read(const char *format, ...);
};

class binary_file_t : public file_t{
    private:
    public:
        binary_file_t(){}
        ~binary_file_t(){
            binary_file_t::close();
        }
        void open(const char *filename, const char mode);
        template <typename type_t>
        int write(const type_t *data){
            return fwrite(data, sizeof(*data), 1, this->file_ptr);
        }
        template <typename type_t>
        int read(type_t *data){
            return fread(data, sizeof(*data), 1, this->file_ptr);
        }
};

// Functions
void __assert_error(const int_t condition, const char *error_msg, 
    const char* filename, const size_t line);

void print(const char *format, ...);
void print(const int_t n);
void print(const size_t n);
void print(const real_t x);
void print(const complex_t z);

void progress_bar(const size_t i, const size_t N, const char *msg);


#endif