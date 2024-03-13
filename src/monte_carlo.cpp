#include "monte_carlo.h"

double MonteCarlo::calculate_uniform(const double a, const double b, std::function<double(double)> integrand)
{
    // Element of single block
    double single_block = 0.;

    for (unsigned int i = 1; i <= m_launch_size; i++)
    {
        double new_element = integrand(m_rng->Rannyu(a, b));
        single_block = (single_block * (i - 1) + new_element) / i; // Incremental Average to reduce error due to overflow
    }

    single_block *= (b - a); // Calcolo Integrale, <f(x)>_a^b * (b-a)

    return single_block;
}

double MonteCarlo::calculate_uniform(const double a, const double b, const unsigned int throws, std::function<double(double)> integrand)
{
    double single_block = 0.;

    for (unsigned int i = 1; i <= throws; i++)
    {
        double new_element = integrand(m_rng->Rannyu(a, b));
        single_block = (single_block * (i - 1) + new_element) / i; // Incremental Average to reduce error due to overflow
    }
    single_block *= (b - a); // Calcolo Integrale, <f(x)>_a^b * (b-a)

    return single_block;
}

std::vector<values> MonteCarlo::calculate_uniform_blocks(const double a, const double b, const unsigned int block_size, std::function<double(double)> integrand)
{
    std::vector v_values(m_n_block, 0.);
    std::vector v_blocks(m_n_block, values{0., 0.});

    for (auto &value : v_values)
    {
        value = calculate_uniform(a, b, block_size, integrand);
    }

    auto begin = v_values.begin();

    for (unsigned int i = 0; i < m_n_block; i++)
    {
        v_blocks[i] = values{calc_mean(begin, begin + i), calc_std(begin, begin + i)};
    }

    return v_blocks;
}

double MonteCarlo::calculate_distribution_convert(const double a, const double b, const unsigned int throws, std::function<double(double)> integrand, std::function<double(double)> inverse_cumulative)
{
    double single_block = 0.;

    for (unsigned int i = 1; i <= throws; i++)
    {
        double new_element = integrand(inverse_cumulative(m_rng->Rannyu(a, b)));
        single_block = (single_block * (i - 1) + new_element) / i; // Incremental Average to reduce error due to overflow

    }
    single_block *= (b - a); // Calcolo Integrale, <f(x)>_a^b * (b-a)

    return single_block;
    ;
}

double MonteCarlo::calculate_distribution_rng(const double a, const double b, const unsigned int throws, std::function<double(double)> integrand, std::function<double(void)> ext_rng)
{
    double single_block = 0.;

    for (unsigned int i = 1; i <= throws; i++)
    {
        double new_element = integrand(ext_rng());
        single_block = (single_block * (i - 1) + new_element) / i; // Incremental Average to reduce error due to overflow
    }
    single_block *= (b - a); // Calcolo Integrale, <f(x)>_a^b * (b-a)

    return single_block;
    ;
}

std::vector<values> MonteCarlo::calculate_dist_conv_block(const double a, const double b, const unsigned int block_size, std::function<double(double)> integrand, std::function<double(double)> inverse_cumulative)
{
    std::vector v_values(m_n_block, 0.);
    std::vector v_blocks(m_n_block, values{0., 0.});

    for (auto &value : v_values)
    {
        value = calculate_distribution_convert(a, b, block_size, integrand, inverse_cumulative);
    }

    auto begin = v_values.begin();

    for (unsigned int i = 0; i < m_n_block; i++)
    {
        v_blocks[i] = values{calc_mean(begin, begin + i), calc_std(begin, begin + i)};
    }

    return v_blocks;
}

std::vector<values> MonteCarlo::calculate_dist_rng_block(const double a, const double b, const unsigned int block_size, std::function<double(double)> integrand, std::function<double(void)> ext_rng)
{
    std::vector v_values(m_n_block, 0.);
    std::vector v_blocks(m_n_block, values{0., 0.});

    for (auto &value : v_values)
    {
        value = calculate_distribution_rng(a, b, block_size, integrand, ext_rng);
    }

    auto begin = v_values.begin();

    for (unsigned int i = 0; i < m_n_block; i++)
    {
        v_blocks[i] = values{calc_mean(begin, begin + i), calc_std(begin, begin + i)};
    }

    return v_blocks;
}
