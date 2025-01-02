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
    // Loading Simulation Options
    std::ifstream config_file("config.jsonc");
    Options opt;
    file_parser(config_file, opt);
    config_file.close();

    // Sampler
    generator_func generator = opt.generator[opt.new_pos_generator];

    // Initialize Metropolis Random Number Generator
    Metropolis metr;
    initializer(*(metr.get_rng_ptr()), opt.first_initializer_row);

    Random rng;
    initializer(rng, opt.second_initializer_row);

    // Movement Parameters
    const double delta = opt.delta;
    const Vector3D start_pos = opt.start_pos;

    // Initialing Blocks
    const std::size_t n_block = opt.n_block;
    std::vector v_mean_r(n_block, 0.);
    std::vector<values> v_radius(n_block, {0., 0.});

    std::cout << "Nr Blocks : " << opt.n_block << '\n'
              << "Nr Block Steps : " << opt.n_block_step << '\n'
              << "Nr Particle Step : " << opt.n_particle_step << '\n'
              << "Nr Thermalization Steps : " << opt.n_thermalization_step << '\n';

    // Output Report Files
    std::ofstream output, output_info;
    if (opt.output_file)
    {
        output.open("output.csv");
    }

    if (opt.output_info)
    {
        output_info.open("info.txt");
        output_info << "Block,Acceptance\n";
    }

    // Thermalization
    auto current_pos = start_pos;
    metr.set_n_step(opt.n_thermalization_step);
    current_pos = metr.generate<Vector3D>(current_pos, opt.convert[opt.func], [&] { return generator(&rng, delta); });

    // Main Loop
    metr.set_n_step(opt.n_particle_step);
    for (size_t i = 0; i < n_block; i++)
    {
        current_pos = (opt.reset) ? start_pos : current_pos;

        std::vector v_instant_pos(opt.n_block_step, 0.);

        for (auto &instant_pos : v_instant_pos)
        {
            current_pos =
                metr.generate<Vector3D>(current_pos, opt.convert[opt.func], [&] { return generator(&rng, delta); });
            instant_pos = current_pos.distance();
        }

        if (opt.output_video)
        {
            std::cout << i << '\t' << current_pos << '\n';
        }

        if (opt.output_info)
        {
            output_info << i << "," << metr.get_acceptance() << '\n';
        }

        v_mean_r[i] = calc_mean(v_instant_pos);
        output << current_pos << '\n';
    }

    auto begin = v_mean_r.begin();

    // Data Blocking and Moving Average
    for (size_t i = 0; i < v_mean_r.size(); i++)
    {
        v_radius[i] = {calc_mean(begin, begin + i), calc_std(begin, begin + i)};
    }

    if (opt.output_mean_radius)
    {
        std::ofstream mean_radius_stream("mean_radius.csv");

        if (!mean_radius_stream)
        {
            std::cerr << "Could Not open mean_radius.dat\nExiting\n";
            exit(-1);
        }
        mean_radius_stream << "idx,instant_r,mean_r,std_r\n";

        for (std::size_t i = 0; i < v_radius.size(); i++)
        {
            mean_radius_stream << i << ',' << v_mean_r[i] << ',' << v_radius[i] << "\n";
        }

        mean_radius_stream.close();
    }

    // Cleaning
    if (output.is_open())
    {
        output.close();
    }
    if (output_info.is_open())
    {
        output_info.close();
    }

    return 0;
}
