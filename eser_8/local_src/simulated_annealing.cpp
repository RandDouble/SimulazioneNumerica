#include "simulated_annealing.h"
#define CHECK_MONTECARLO_MOVEMENT

#pragma region IO_FUNCTIONS

void SimulatedAnnealing::initialize_rng(std::size_t rows_to_skip)
{
    m_sampler.get_rng_ptr()->Initializer("seed.in", "Primes", rows_to_skip);
    m_sampler.set_n_step(1); // We want montecarlo to generate a step at time
}

void SimulatedAnnealing::write_results(const std::string &filename)
{
    std::ofstream file(filename);
    file << param_history.rdbuf();
    file.close();
}

void SimulatedAnnealing::save_seed()
{
    m_sampler.get_rng_ptr()->SaveSeed("montecarlo_seed.out");
}

#pragma endregion IO_FUNCTIONS

#pragma region TEMPERING

size_t SimulatedAnnealing::tempering_cycle()
{
    size_t accepted_temp_step = 0;
    double starting_position = 0.;
    values old_energy = evaluate_energy(starting_position);
    PsiParam old_param = m_wave.param();

    for (size_t j = 0; j < n_tempering_step; j++)
    {
        assert(!std::isinf(old_energy.value) && "Old energy is Infinite");
        assert(!std::isnan(old_energy.value) && "Old energy is NaN");
        m_wave.param(propose_new_param()); // Update wave function parameters
        const values new_energy = evaluate_energy(starting_position);
        assert(!std::isnan(new_energy.value) && "New energy is NaN");
        assert(!std::isinf(new_energy.value) && "New energy is Infinite");

        const double delta_e = new_energy.value - old_energy.value; // Calculate energy difference

        assert(!std::isnan(delta_e) && "Delta energy is NaN");
        if (m_sampler.get_rng_ptr()->Rannyu() < boltzmann_weight(delta_e))
        {
#ifndef CHECK_TEMPERING_MOVEMENT
            std::cout << "Step Accepted, current delta energy : " << delta_e << '\n';
#endif
            old_energy = new_energy;
            old_param = m_wave.param();
            accepted_temp_step++;
        }
        else
        {
#ifndef CHECK_TEMPERING_MOVEMENT
            std::cout << "Step Rejected, current delta energy : " << delta_e << '\n';
#endif
            m_wave.param(old_param); // Restore old parameters
        }

#ifndef CHECK_TEMPERING_MOVEMENT
        std::cout << "Current Energy" << std::setw(8) << old_energy.value << " +- " << std::setw(8) << old_energy.error
                  << '\n';
#endif // CHECK_TEMPERING_MOVEMENT

        param_history << get_temperature() << ',' << m_wave.mu() << ',' << m_wave.sigma() << ',' << old_energy << '\n';
    }
    return accepted_temp_step;
}

void SimulatedAnnealing::temperature_cycle()
{
    param_history.clear();
    param_history.precision(10);
    param_history << "temperature,mu,sigma,energy,err_energy\n";

    for (unsigned int i = 0; i <= n_temp_step; i++)
    {
        std::cout << "Actual Temp : " << std::setw(5) << get_temperature() << " Actual Beta : " << std::setw(5)
                  << get_beta() << '\n';
        size_t accepted_temp_step = tempering_cycle();
        std::cout << "\nAccepted Step : " << std::setw(4) << accepted_temp_step << '\n';
        if (i == n_temp_step)
            break;
        update_temperature();
    }
}

#pragma endregion TEMPERING

#pragma region ENERGY

values SimulatedAnnealing::evaluate_energy(const double x_start)
{
    std::vector block_result(monte_carlo_block, 0.);
#ifndef CHECK_MONTECARLO_MOVEMENT
    const double weight = 1. / monte_carlo_step;
    double acceptance = 0.;
#endif
    m_sampler.set_n_step(1);

    auto translation = [&]() {
        return m_sampler.get_rng_ptr()->Rannyu(-scaling_movement_law(), scaling_movement_law());
    };

    auto PDF = [&](double x) { return m_wave.PDF(x); };

    m_sampler.set_current_position(x_start);
    double x_i = x_start;

    for (auto &&block : block_result)
    {
        std::vector partial_result(monte_carlo_step, 0.);
        for (auto &partial : partial_result)
        {
#ifndef CHECK_MONTECARLO_MOVEMENT
            double x_old = x_i;
#endif
            x_i = m_sampler.generate(PDF, translation);

            partial = m_wave.hamiltonian(x_i);

#ifndef CHECK_MONTECARLO_MOVEMENT
            acceptance += (x_i != x_old);
#endif
        }
        block = calc_mean(partial_result);
        assert(!std::isnan(block) && "Block is NaN");
        assert(!std::isinf(block) && "Block is Infinite");
    }
#ifndef CHECK_MONTECARLO_MOVEMENT
    std::cout << "Montecarlo acceptance : " << acceptance * weight / monte_carlo_block << '\n';
#endif // CHECK_MONTECARLO_MOVEMENT

    values result{.value = calc_mean(block_result), .error = calc_std(block_result)};

    return result;
}

double SimulatedAnnealing::boltzmann_weight(const double delta_energy)
{
    return std::exp(-get_beta() * delta_energy);
}

#pragma endregion ENERGY

#pragma region UPDATE

void SimulatedAnnealing::update_temperature()
{ // I want to update temperature Linearly...
    // new_t = (1 / old_beta) + delta_t
    // 1/new_beta = (1 / old_beta) + delta_t
    // 1 / nesw_beta = (1 + old_beta * delta_t) / old_beta
    // new_beta = old_beta / (1 + old_beta * delta_t)
    // the rate is saved as a positive number so it need a - sign
    // new_beta = old_beta / (1 - old_beta * delta_t)

    // actual_beta = actual_beta / (1. - actual_beta * cooling_rate);

    // Found empirically that previous formula is not working well enough
    // Found this article in lecterature
    // https://www.sciencedirect.com/science/article/pii/S1568494617305768

    // t_new = cooling_rate * t_old
    // 1 / beta_new = cooling_rate / beta_old
    // beta_new = beta_old / cooling_rate

    m_actual_beta /= m_cooling_rate;

    // Adding parameter movement reduction

    delta_move.mu -= m_reduction_rate.mu;
    delta_move.sigma -= m_reduction_rate.sigma;
    std::cout << "New delta mu : " << delta_move.mu << " New delta sigma : " << delta_move.sigma << '\n';
}

PsiParam SimulatedAnnealing::propose_new_param()
{

    PsiParam new_param = m_wave.param();

    new_param.mu += m_sampler.get_rng_ptr()->Rannyu(-delta_move.mu, delta_move.mu);
    new_param.sigma += m_sampler.get_rng_ptr()->Rannyu(-delta_move.sigma, delta_move.sigma);

    assert(!std::isinf(new_param.mu) && "mu is infinite");
    assert(!std::isinf(new_param.sigma) && "sigma is infinite");

    assert(!std::isnan(new_param.mu) && "mu is NaN");
    assert(!std::isnan(new_param.sigma) && "sigma is NaN");

    return new_param;
}

double SimulatedAnnealing::scaling_movement_law() const
{
    // const double gamma = 0.6; // Experimentally found
    // return integration_delta * (1 - gamma * std::exp(-get_beta()));
    return integration_delta;
}

#pragma endregion UPDATE
