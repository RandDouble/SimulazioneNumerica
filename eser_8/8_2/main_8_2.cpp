#include <iostream>
#include <vector>

#include <metropolis.h>
#include <utilities.h>

#include "simulated_annealing.h"
#include "wave_function.h"

// constexpr double delta_mu = 0.1;
// constexpr double delta_sigma = 0.1;

int main()
{
    constexpr double starting_temperature = 1000;
    constexpr double ending_temperature = 1e-3;
    constexpr std::size_t n_temperature_step = 100;
    constexpr std::size_t same_temperature_step = 10;

    constexpr unsigned int n_blocks = 10000; // Number of blocks for data blocking
    constexpr unsigned int n_step = 100;
    constexpr PsiParam delta_params
        = {.mu = 0.1, .sigma = 0.1};         // Number of step for each block
    constexpr double integration_delta = 1.; // Movement of metropolis algorithm
    constexpr PsiParam initial_param{.mu = 1., .sigma = 0.5};
    constexpr PsiParam parameter_reduction_rate
        = {.mu = delta_params.mu / (n_temperature_step + 1),
           .sigma = delta_params.sigma / (n_temperature_step + 1)};

    SimulatedAnnealing simulator(starting_temperature,
                                 ending_temperature,
                                 n_temperature_step,
                                 same_temperature_step,
                                 integration_delta,
                                 delta_params,
                                 parameter_reduction_rate,
                                 initial_param);
    simulator.initialize_rng();

    simulator.set_metropolis_step(n_step);
    simulator.set_metropolis_block(n_blocks);

    simulator.temperature_cycle();
    simulator.write_results("./output.csv");
    simulator.save_seed();

    return 0;
}
