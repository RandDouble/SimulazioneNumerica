#ifndef  __OPTIONS__
#define __OPTIONS__

#include <fstream>

#include "vector3D.h"

struct Options
{
    std::size_t n_block{1};
    std::size_t n_step{10};
    double delta{1.};
    Vector3D start_pos{1., 0., 0.};
    int first_initializer_row{0}, second_initializer_row{0};
};

void file_parser(std::ifstream &in, Options &opt);


#endif // __OPTIONS__
