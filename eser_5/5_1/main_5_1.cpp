#include <iomanip>
#include <iostream>

#include <vector>

#include "metropolis.h"
#include "utilities.h"
#include "wave_function.h"
#include "options.h"

// int argc, char** argv

// Remember in metropolis.h there is a TEST Macro for enabling output

int main()
{
    std::ifstream config_file("config.dat");
    Options opt;
    file_parser(config_file, opt);
    config_file.close();

    Metropolis metr;
    initializer(*(metr.get_rng_ptr()), opt.first_initializer_row);
    metr.set_n_step(opt.n_step);

    const std::size_t n_block = opt.n_block;

    Random rng;
    initializer(rng, opt.second_initializer_row);

    const double delta = opt.delta;

    const Vector3D start_pos = opt.start_pos;

    for (size_t i = 0; i < n_block; i++)
    {
        std::cout << i << '\t'
                  << metr.generate<Vector3D>(start_pos, ground_state, [&]
                                             { return Vector3D::generate_unif(&rng, delta); })
                  << '\n';
    }

    return 0;
}

