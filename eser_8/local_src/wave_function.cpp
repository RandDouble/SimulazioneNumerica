#include "wave_function.h"

double WaveFunction::operator()(const double x) const
{
    const double exponent = -0.5 * (x * x + mu() * mu()) * m_inv_sigma_squared;
    const double cosh_arg = mu() * x * m_inv_sigma_squared;

    if (exponent < -708.4) // exp would underflow
        return 0.;

    assert(!std::isnan(std::exp(exponent)) && "Wavefunction exp is NaN");
    assert(!std::isinf(std::exp(exponent)) && "Wavefunction exp is Infinite");
    assert(!std::isnan(std::cosh(cosh_arg)) && "Wavefunction cosh is NaN");
    assert(!std::isinf(std::cosh(cosh_arg)) && "Wavefunction cosh is Infinite");
    assert(!std::isnan(2. * std::exp(exponent) * std::cosh(cosh_arg)) && "Wavefunction is NaN");
    assert(!std::isinf(2. * std::exp(exponent) * std::cosh(cosh_arg)) && "Wavefunction is Infinite");

    return 2. * std::exp(exponent) * std::cosh(cosh_arg);
}

double WaveFunction::second_der(const double x) const
{
    const double a_2 = (x - mu()) * (x - mu()) * m_inv_sigma_squared;
    const double b_2 = (x + mu()) * (x + mu()) * m_inv_sigma_squared;
    const double exp_a = std::exp(-0.5 * a_2);
    const double exp_b = std::exp(-0.5 * b_2);
    double result = exp_a * (a_2 - 1.) + exp_b * (b_2 - 1.);
    result *= m_inv_sigma_squared;

    assert(!std::isnan(result) && "Second derivative is NaN");
    assert(!std::isinf(result) && "Second derivative is Infinite");
    return result;
}

double WaveFunction::kinetic(const double x) const
{
    // Obtained from \frac{\Psi''}{\Psi}, the result is -0.5 * \frac{\Psi''}{\Psi}

    const double first_term = (x * x + mu() * mu()) * m_inv_sigma_squared - 1.;
    const double tanh_arg = mu() * x * m_inv_sigma_squared;
    const double res = m_inv_sigma_squared * (first_term - 2. * tanh_arg * std::tanh(tanh_arg));

    assert(!std::isnan(res) && "Kinetic energy is NaN");
    assert(!std::isinf(res) && "Kinetic energy is Infinite");

    return -0.5 * res;
}

double WaveFunction::PDF(const double x) const
{
    const double exponent = -(x * x + mu() * mu()) * m_inv_sigma_squared;
    const double cosh_arg = 2. * mu() * x * m_inv_sigma_squared;
    if (exponent < -708.4) // Exp would underflow, so cut your losses
        return 0.;

    return std::exp(exponent) * (std::cosh(cosh_arg) + 1.);
}

double WaveFunction::normalized_PDF(const double x) const
{
    const double norm = 1. / (sigma() * std::sqrt(M_PI) * (1.0 + std::exp(-(mu() * mu()) * m_inv_sigma_squared)));
    const double exponent = -(x * x + mu() * mu()) * m_inv_sigma_squared;
    const double cosh_arg = 2. * mu() * x * m_inv_sigma_squared;
    if (exponent < -708.4) // Exp would underflow, so cut your losses
        return 0.;

    return norm * std::exp(exponent) * (std::cosh(cosh_arg) + 1.);
}

double WaveFunction::hamiltonian(const double x) const
{
    double kin = kinetic(x);
    double pot = potential(x);

    assert(!std::isnan(kin + pot) && "Energy is NaN");
    assert(!std::isinf(kin + pot) && "Energy is Infinite");

    return kin + pot;
}

double WaveFunction::potential(const double x) const
{
    const double x_square = x * x;
    return (x_square - 2.5) * x_square;
}
