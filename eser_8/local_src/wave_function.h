#include <cassert>
#include <cmath>
#include <numeric>

#ifndef __WAVE_FUNCTION__
#define __WAVE_FUNCTION__

struct PsiParam
{
    double mu, sigma;
};

class WaveFunction
{
private:
    PsiParam m_param{.mu = 0., .sigma = 1.};
    double inv_sigma_squared = 1.;

public:
    WaveFunction() = default;
    WaveFunction(const PsiParam par) : m_param(par) {}
    WaveFunction(const double mu, const double sigma) :
        m_param({.mu = mu, .sigma = sigma}),
        inv_sigma_squared(1. / (sigma * sigma))
    {}

    double operator()(const double x) const;
    double second_der(const double x) const;
    double PDF(const double x) const;
    double hamiltonian(const double x) const;
    double kinetic(const double x) const;
    double potential(const double x) const;

    // Setters
    void mu(const double mu) { m_param.mu = mu; }
    void sigma(const double sigma)
    {
        assert(sigma != 0. && "Sigma cannot be zero");
        m_param.sigma = sigma;
        inv_sigma_squared = 1. / (sigma * sigma);
    }
    void param(const PsiParam& par)
    {
        mu(par.mu);
        sigma(par.sigma);
    }

    // Getters
    double mu() const { return m_param.mu; }
    double sigma() const { return m_param.sigma; }
    PsiParam param() const { return m_param; }
};

#endif // __WAVE_FUNCTION__
