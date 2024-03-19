#include "finanza.h"


double FinancialOption::call_option_price_conv(double ending_price) const
{
    const double discount = std::exp(-m_r * m_T);
    return discount * std::max(0., ending_price - m_K);
}

double FinancialOption::put_option_price_conv(double ending_price) const
{
    const double discount = std::exp(-m_r * m_T);
    return discount * std::max(0., m_K - ending_price);
}

double FinancialOption::call_option_final_mean()
{
    double mean = 0.;
    for (unsigned int i = 1; i <= m_extraction; i++)
    {
        double S = ending_price_final();
        double C = call_option_price_conv(S);
        mean = (mean * (i - 1) + C) / i;
    }

    return mean;
}

double FinancialOption::put_option_final_mean()
{
    double mean = 0.;
    for (unsigned int i = 1; i <= m_extraction; i++)
    {
        double S = ending_price_final();
        double P = put_option_price_conv(S);
        mean = (mean * (i - 1) + P) / i;
    }

    return mean;
}
double FinancialOption::call_option_step_mean(const unsigned int n_step)
{
    double mean = 0.;
    for (unsigned int i = 1; i <= m_extraction; i++)
    {
        double S = ending_price_step(n_step);
        double C = call_option_price_conv(S);
        mean = (mean * (i - 1) + C) / i;
    }

    return mean;
}

double FinancialOption::put_option_step_mean(const unsigned int n_step)
{
    double mean = 0.;
    for (unsigned int i = 1; i <= m_extraction; i++)
    {
        double S = ending_price_step(n_step);
        double P = put_option_price_conv(S);
        mean = (mean * (i - 1) + P) / i;
    }

    return mean;
}

double FinancialOption::ending_price_final()
{
    const double W = Wiener(m_rng, m_T);
    const double prefix = m_r - 0.5 * m_sigma * m_sigma;
    return S_0 * std::exp(prefix * m_T + m_sigma * W);
}

double FinancialOption::ending_price_step(const unsigned int n_step)
{
    double S = S_0;
    double d_t = m_T / static_cast<double>(n_step);
    const double prefix = m_r - 0.5 * m_sigma * m_sigma;

    for (unsigned int i = 1; i <= n_step; i++)
    {
#ifndef NDEBUG
        std::cout << "step: " << i << "\t" << S << "\n";
#endif // NDEBUG
        double Z = Wiener_step(m_rng);
        S = S * std::exp(prefix * d_t + m_sigma * Z * std::sqrt(d_t));
    }

    return S;
}
