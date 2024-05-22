#include "wave_function.h"

double WaveFunction::operator()(const double x) const
{
    const double inv_sigma_squared = 1. / (sigma() * sigma());
    const double exponent = -0.5 * (x * x + mu() * mu()) * inv_sigma_squared;
    const double cosh_arg = mu() * x * inv_sigma_squared;

    return 2. * std::exp(exponent) * std::cosh(cosh_arg);
}

double WaveFunction::second_der(const double x) const
{
    const double inv_sigma_squared = 1. / (sigma() * sigma());
    const double a_2 = (x - mu()) * (x - mu()) * inv_sigma_squared;
    const double b_2 = (x + mu()) * (x + mu()) * inv_sigma_squared;
    const double exp_a = std::exp(-0.5 * a_2);
    const double exp_b = std::exp(-0.5 * b_2);
    double result = exp_a * (a_2 - 1.) + exp_b * (b_2 - 1.);
    result *= inv_sigma_squared;
    return result;
}

double WaveFunction::PDF(const double x) const
{
    const double inv_sigma_squared = 1. / (sigma() * sigma());
    const double exponent = -(x * x + mu() * mu()) * inv_sigma_squared;
    const double cosh_arg = 2. * mu() * x * inv_sigma_squared;

    return 2. * std::exp(exponent) * (std::cosh(cosh_arg) + 1.);
}

double WaveFunction::hamiltonian(const double x) const
{
    double kinetic = -0.5 * second_der(x);
    double pot = potential(x);
    double psi = operator()(x);

    return kinetic / psi + pot;
}

double potential(const double x)
{
    const double x_square = x * x;
    return (x_square - 2.5) * x_square;
}

double wave_function(const double x, const double mu, const double sigma)
{
    const double inv_sigma_squared = 1. / (sigma * sigma);
    const double exponent = -0.5 * (x * x + mu * mu) * inv_sigma_squared;
    const double cosh_arg = mu * x * inv_sigma_squared;

    return 2. * std::exp(exponent) * std::cosh(cosh_arg);
}

double wave_function_prob(const double x, const double mu, const double sigma)
{
    const double inv_sigma_squared = 1. / (sigma * sigma);
    const double exponent = -(x * x + mu * mu) * inv_sigma_squared;
    const double cosh_arg = 2. * mu * x * inv_sigma_squared;

    return 2. * std::exp(exponent) * (std::cosh(cosh_arg) + 1.);
}

double wave_function_second(const double x, const double mu, const double sigma)
{
    const double inv_sigma_squared = 1. / (sigma * sigma);
    const double a_2 = (x - mu) * (x - mu) * inv_sigma_squared;
    const double b_2 = (x + mu) * (x + mu) * inv_sigma_squared;
    const double exp_a = std::exp(-0.5 * a_2);
    const double exp_b = std::exp(-0.5 * b_2);
    double result = exp_a * (a_2 - 1.) + exp_b * (b_2 - 1.);
    result *= inv_sigma_squared;
    return result;
}

double hamiltonian(const double x, const double mu, const double sigma)
{
    double kinetic = -0.5 * wave_function_second(x, mu, sigma);
    double pot = potential(x);
    double psi = wave_function(x, mu, sigma);

    return kinetic / psi + pot;
}

double boltzman_weight(const double beta, const double energy)
{
    return std::exp(-beta * energy);
}
