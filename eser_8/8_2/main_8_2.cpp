#include <array>
#include <iostream>
#include <vector>

#include <metropolis.h>
#include <utilities.h>

#include "simulated_annealing.h"
#include "wave_function.h"

int main()
{
#pragma region CONSTANTS
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
#pragma endregion CONSTANTS

#pragma region SIMULATION

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

#pragma endregion SIMULATIONs

    // Calulating Energy
#pragma region ENERGY_CALCULATION

    WaveFunction final_wave = simulator.get_params();

    std::cout << "Final Parameters : " << std::setw(10) << final_wave.mu()
              << std::setw(10) << final_wave.sigma() << '\n';

    MetropolisMemory<double> final_sampler(0.);

    auto PDF = [&final_wave](double x) {
        // std::cout << "Final Wave Params : " << final_wave.mu() << " " <<
        // final_wave.sigma() << '\n';
        return final_wave.PDF(x);
    };
    auto sampler = [&final_sampler, &final_wave]() {
        // std::cout << "Final Param Sigma : " << final_param.sigma << '\n';
        return final_sampler.get_rng_ptr()->Rannyu(-final_wave.sigma(),
                                                   final_wave.sigma());
    };

    final_sampler.get_rng_ptr()->Initializer();
    final_sampler.set_n_step(1);

    constexpr std::size_t energy_block = 1000;
    constexpr std::size_t energy_step = 100;

    std::vector<double> v_energies(energy_block, 0.);

    for (auto&& energy : v_energies)
    {
        std::vector block_energy(energy_step, 0.);
        for (auto&& energy_step : block_energy)
        {
            energy_step
                = final_wave.hamiltonian(final_sampler.generate(0., PDF, sampler));
        }
        energy = calc_mean(block_energy);
    }

    // Data Blocking

    auto begin = v_energies.begin();

    std::ofstream fout("./energy_evaluation.csv");

    fout << "Block,Energy,Error\n";

    for (auto ptr = begin; ptr != v_energies.end(); ptr++)
    {

        values val = {calc_mean(begin, ptr), calc_std(begin, ptr)};
        fout << std::distance(begin, ptr) << ", " << val << "\n";
    }

    std::cout << "Energy : " << calc_mean(v_energies) << " +- " << calc_std(v_energies)
              << '\n';

    fout.close();

#pragma endregion ENERGY_CALCULATION

#pragma region SAMPLING_PSI2

    // Plotting Psi^2  for Debugging purposes
    // std::cout << "Plotting Psi^2\n"
    //           << "Mu : " << final_wave.mu() << '\n'
    //           << "Sigma : " << final_wave.sigma() << '\n'
    //           << "1/Sigma^2 : " << final_wave.get_inv_sigma_squared() << '\n';
    fout.open("./wave_function.csv");

    for (double x = -3.; x <= 3.; x += 0.01)
    {
        fout << x << ", " << final_wave.normalized_PDF(x) << "\n";
    }
    fout.close();

    // Sampling Psi^2

    constexpr std::size_t n_sample = 100000;
    std::array<double, n_sample> sampled_position;
    double last_sampled_position = 0.;

    final_sampler.set_n_step(10);

    for (auto& sample : sampled_position)
    {
        sample = final_sampler.generate(last_sampled_position, PDF, sampler);
        last_sampled_position = sample;
    }

    fout.open("./sampled.csv");
    fout << "x\n";
    for (const auto& sample : sampled_position)
    {
        fout << sample << "\n";
    }
    fout.close();

#pragma endregion SAMPLING_PSI2
    return 0;
}
