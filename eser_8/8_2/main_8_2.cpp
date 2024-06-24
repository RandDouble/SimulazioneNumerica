#include <iostream>
#include <vector>

#include <metropolis.h>
#include <utilities.h>

#include "wave_function.h"

void update_wave_param(WaveFunction &wave, Random &rng);
double update_temp(const double old_beta);
double calc_expexted_energy(Metropolis &metro, Random &rng, const double delta, const WaveFunction &wave, const double x_start);

constexpr double starting_temperature = 10000;
constexpr double ending_temperature = 1e-5;
constexpr std::size_t n_temperature_step = 600;
constexpr std::size_t same_temperature_step = 10;

constexpr double d_temperature = (starting_temperature - ending_temperature) / n_temperature_step; // Increase in temperature
constexpr unsigned int n_blocks = 1000;                                                            // Number of blocks for data blocking
constexpr unsigned int n_step = 1000;                                                              // Number of step for each block
// constexpr double delta_mu = 0.1;
// constexpr double delta_sigma = 0.1;

int main()
{
    constexpr double starting_beta = 1. / starting_temperature;
    constexpr double delta = 1.;   // Movement of metropolis algorithm
    constexpr double x_start = 0.; // Starting Position

    double actual_beta = starting_beta;
    int accepted_temp_step = 0;

    PsiParam initial_param{.mu = 0., .sigma = 1.};
    WaveFunction wave(initial_param);

    // Initializing metropolis
    Metropolis metro;
    Random rng;
    initializer(*metro.get_rng_ptr());
    initializer(rng, 10); // Differente lista di numeri da cui generare

    // Set how many step to do for metropolis sampler
    metro.set_n_step(1);

    double old_energy = calc_expexted_energy(metro, rng, delta, wave, x_start);

    for (size_t i = 0; i < n_temperature_step; i++)
    {
        for (size_t j = 0; j < same_temperature_step; j++)
        {
            // I'll just write explicitly this montecarlo step
            const double new_energy = calc_expexted_energy(metro, rng, delta, wave, x_start);
            double delta_e = new_energy - old_energy;
            if (rng.Rannyu() < boltzman_weight(actual_beta, delta_e))
            {
                old_energy = new_energy;
                accepted_temp_step++;
            }
        }
        std::cout << "Actual Temp : " << std::setw(5) << (1. / actual_beta)
                  << "\nAccepted Step : " << std::setw(4) << accepted_temp_step << '\n';
        accepted_temp_step = 0;
        actual_beta = update_temp(actual_beta);
    }

    return 0;
}

void update_wave_param(WaveFunction &wave, Random &rng)
{
    wave.mu(wave.mu() + rng.Rannyu());
    wave.sigma(wave.sigma() + rng.Rannyu());
}

double update_temp(const double old_beta)
{
    // I want to update temperature Linearly...
    // new_t = (1 / old_beta) + delta_t
    // 1/new_beat = (1 / old_beta) + delta_t
    // 1 / new_beta = (1 + old_beta * delta_t) / old_beta
    // new_beta = old_beta / (1 + old_beta * delta_t)

    const double result = old_beta / (1. + old_beta * d_temperature);
    return result;
}

double calc_expexted_energy(Metropolis &metro, Random &rng, const double delta, const WaveFunction &wave, const double x_start)
{
    std::vector block_result(n_blocks, 0.);
    const double weight = 1. / n_step;

    auto sampler = [&]()
    { return rng.Rannyu(-delta, delta); };

    auto PDF = [&](double x)
    { return wave.PDF(x); };

    double x_i = x_start;

    for (auto &&block : block_result)
    {
        for (unsigned int i = 0; i < n_step; i++)
        {
            x_i = metro.generate<double>(x_i, PDF, sampler);
            block += wave.hamiltonian(x_i);
        }
        block *= weight;
    }

    return calc_mean(block_result);
}