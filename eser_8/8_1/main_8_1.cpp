#include <iostream>
#include <vector>

#include <metropolis.h>
#include <utilities.h>

#include "wave_function.h"

int main()
{
    constexpr unsigned int n_step = 1000;
    constexpr unsigned int n_blocks = 1000;
    constexpr double delta = 1.;
    constexpr double mu = 0.;
    constexpr double sigma = 1.;
    constexpr double x_start = 0.;

    Metropolis metro;
    Random rng;
    initializer(*metro.get_rng_ptr());
    initializer(rng, 10); // Differente lista di numeri da cui generare

    metro.set_n_step(1);

    std::vector block_result(n_blocks, 0.);
    std::vector result(n_blocks, values{0., 0.});

    auto sampler = [&]()
    { return rng.Rannyu(-delta, delta); };

    auto PDF = [&](double x)
    { return wave_function_prob(x, mu, sigma); };

    double x_i = x_start;

    for (auto &&block : block_result)
    {
        for (unsigned int i = 0; i < n_step; i++)
        {
            x_i = metro.generate<double>(x_i, PDF, sampler);
            block += hamiltonian(x_i, mu, sigma);
        }
        block /= n_step;
    }

    auto begin = block_result.begin();
    for (auto end = block_result.begin(); end != block_result.end(); end++)
    {
        auto idx = std::distance(begin, end);
        result[idx].error = calc_std(begin, end);
        result[idx].value = calc_mean(begin, end);
    }

    std::ofstream fout("res/energy.csv");
    fout << "average,error\n";
    if (fout)
        print_file(fout, result);
    else
    {
        std::cerr << "Fail State\n";
        fout.close();
        return -1;
    }

    return 0;
}
