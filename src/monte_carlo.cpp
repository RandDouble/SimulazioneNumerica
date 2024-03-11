#include "monte_carlo.h"

values MonteCarlo::calculate_uniform(const double a, const double b, std::function<double(double)> integrand)
{
    std::vector<double> v_single_block(m_n_block);

    for (auto &single_block : v_single_block)
    {
        for (unsigned int i = 0; i < m_launch_size; i++)
        {
            single_block += (integrand(m_rng->Rannyu(a, b)) / m_launch_size); // Calcolo Media singolo blocco
        }
        single_block *= (b - a); // Calcolo Integrale, <f(x)>_a^b * (b-a)
    }

    values val{calc_mean(v_single_block), calc_std(v_single_block)};
    return val;
}

values MonteCarlo::calculate_uniform(const double a, const double b, const unsigned int total_throws, std::function<double(double)> integrand)
{
    std::vector<double> v_single_block(m_n_block);

    unsigned int launch_size = total_throws / m_n_block;

    for (auto &single_block : v_single_block)
    {
        for (unsigned int i = 0; i < launch_size; i++)
        {
            single_block += (integrand(m_rng->Rannyu(a, b)) / launch_size); // Calcolo Media singolo blocco
        }
        single_block *= (b - a); // Calcolo Integrale, <f(x)>_a^b * (b-a)
    }

    values val{calc_mean(v_single_block), calc_std(v_single_block)};

    return val;
}

std::vector<values> MonteCarlo::calculate_uniform_step(const double a, const double b, const unsigned int block_size, unsigned int total_throws, std::function<double(double)> integrand)
{
    if (total_throws % block_size != 0)
    {
        std::cout << "Rescaling total_throws to be a multiple of block_size\n";
        total_throws -= total_throws % block_size;
    }

    m_n_block = total_throws / block_size;
    // Container for integral values computed in a single block
    std::vector<double> v_single_integral(m_n_block);

    // Container for integral values computed with associated error
    std::vector<values> v_single_block_error(m_n_block);

    for (auto &single_integral : v_single_integral)
    {
        for (std::size_t i = 1; i <= block_size; i++)
        {
            single_integral = (single_integral * (i - 1) + (integrand(m_rng->Rannyu(a, b)))) / i; // Incremental Mean to reduce truncation error
        }
        single_integral *= (b - a); // Now we have an integral
    }

    auto single_integral_beginning = v_single_integral.begin();
    for (std::size_t i = 0; i < m_n_block; ++i)
    {
        v_single_block_error[i] = {calc_mean(single_integral_beginning, single_integral_beginning + i), calc_std(single_integral_beginning, single_integral_beginning + i)};
    }
    return v_single_block_error;
}
