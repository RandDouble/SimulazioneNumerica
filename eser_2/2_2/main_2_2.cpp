#include <algorithm>
#include <array>
#include <execution>
#include <fstream>
#include <iostream>
#include <vector>

#include "random.h"
#include "random_walk.h"
#include "utilities.h"

/*
    1. on a cubic lattice with lattice constant $a=1$;
    at each discrete time the walker makes a forward or backward step of length equal to $a$
    in one of the 3 principal directions of the lattice: $x$, $y$ or $z$

    2. in the continuum; at each discrete time the walker makes a step of length equal to $a(=1)$
    along a **random direction** obtained by sampling **uniformly**
    the solid angle: $\theta \in [0,\pi]$ and $\phi \in [0,2\pi]$
 */

std::vector<values> walking_the_walker_discrete(Random_Walk &dog, const std::size_t time_steps, const std::size_t n_throws, const std::size_t n_blocks);
std::vector<values> walking_the_walker_continuos(Random_Walk &dog, const std::size_t time_steps, const std::size_t n_throws, const std::size_t n_blocks);

int main()
{
    Random rng;
    initializer(rng);
    Random_Walk walker(1., {0., 0., 0.}, &rng);

    constexpr std::size_t n_time_steps = 100;
    constexpr std::size_t n_runs = 10000;
    constexpr std::size_t n_blocks = n_runs / n_time_steps;

    // Point 1

    // Vogliamo visualizzazioni fighe per vedere se tutto funziona
    std::ofstream f_off("random_walk_lattice.csv");

    f_off << "x,y,z\n";
    for (size_t i = 0; i < n_time_steps; i++)
    {
        walker.discrete_single_walk();
        f_off << walker.get_actual_pos() << "\n";
    }

    f_off.close();
    walker.reset();

    // Vogliamo visualizzazioni fighe per vedere se tutto funziona
    f_off.open("random_walk_continuos.csv");
    for (size_t i = 0; i < n_time_steps; i++)
    {
        walker.continuos_single_walk();
        f_off << walker.get_actual_pos() << "\n";
    }
    f_off.close();

    // Adesso facciamo il vero punto.
    // Conservo anche la posizione con 0 step eseguiti
    std::cout << "Inizio Discreto\n";
    auto res = walking_the_walker_discrete(walker, n_time_steps, n_runs, n_blocks);

    f_off.open("result_lattice.csv");
    print_file(f_off, res);
    f_off.close();

    // Continuos Walker
    std::cout << "Inizio Continuo\n";
    res = walking_the_walker_continuos(walker, n_time_steps, n_runs, n_blocks);
    f_off.open("result_continuos.csv");
    print_file(f_off, res);
    f_off.close();

    return 0;
}

std::vector<values> walking_the_walker_discrete(Random_Walk &dog, const std::size_t n_time_steps, const std::size_t n_runs, const std::size_t n_blocks)
{
    const unsigned int throw_per_block = n_runs / n_blocks;
    std::vector res(n_time_steps + 1, values({0., 0.}));

    for (std::size_t time_step = 0; time_step <= n_time_steps; time_step++)
    {
        std::vector positions(n_blocks, 0.);
        for (auto &pos : positions)
        {
            std::vector block_positions(throw_per_block, 0.);
            for (auto &block_value : block_positions)
            {
                dog.reset();
                dog.discrete_walk(time_step);
                block_value = dog.get_distance_squared_from_origin();
            }
            pos = std::sqrt(calc_mean(block_positions));
        }
        res[time_step] = {calc_mean(positions), calc_std(positions)};
    }
    return res;
}

std::vector<values> walking_the_walker_continuos(Random_Walk &dog, const std::size_t n_time_steps, const std::size_t n_runs, const std::size_t n_blocks)
{
    const unsigned int throw_per_block = n_runs / n_blocks;

    std::vector res(n_time_steps + 1, values({0., 0.}));

    for (std::size_t time_step = 0; time_step <= n_time_steps; time_step++)
    {
        std::vector positions(n_blocks, 0.);

        for (auto &pos : positions)
        {
            double position_accumulator = 0.;
            for (std::size_t extraction = 0; extraction < throw_per_block; extraction++)
            {
                dog.reset();
                dog.continuos_walk(time_step);
                position_accumulator += dog.get_distance_squared_from_origin();
            }
            position_accumulator /= throw_per_block;
            pos = sqrt(position_accumulator);
        }
        res[time_step] = {calc_mean(positions), calc_std(positions)};
    }
    return res;
}
