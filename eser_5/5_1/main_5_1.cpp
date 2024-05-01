#include <iomanip>
#include <iostream>

#include <vector>

#include "metropolis.h"
#include "options.h"
#include "utilities.h"
#include "wave_function.h"

// Remember in metropolis.h there is a TEST Macro for enabling output

int main()
{
    std::ifstream config_file("config.dat");
    Options opt;
    file_parser(config_file, opt);
    config_file.close();
    std::ofstream output;

    Metropolis metr;
    initializer(*(metr.get_rng_ptr()), opt.first_initializer_row);
    metr.set_n_step(opt.n_step);

    const std::size_t n_block = opt.n_block;

    Random rng;
    initializer(rng, opt.second_initializer_row);

    const double delta = opt.delta;

    const Vector3D start_pos = opt.start_pos;

    if (opt.output_file)
    {
        output.open("output.csv");
    }

    for (size_t i = 0; i < n_block; i++)
    {
        auto current_pos = metr.generate<Vector3D>(start_pos, opt.convert[opt.func], [&]
                                                   { return Vector3D::generate_unif(&rng, delta); });
        if (opt.output_video)
        {
            std::cout << i << '\t'
                      << current_pos
                      << '\n';
        }

        output << current_pos << '\n';
    }

    if (output.is_open())
    {
        output.close();
    }

    return 0;
}
