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
    double m_inv_sigma_squared{1.};

public:
    WaveFunction() = default;
    WaveFunction(const PsiParam& par) :
        m_param(par),
        m_inv_sigma_squared(1. / (par.sigma * par.sigma))
    {}

    WaveFunction(PsiParam&& par) : m_param(std::move(par))
    {
        m_inv_sigma_squared = 1. / (m_param.sigma * m_param.sigma);
    }

    WaveFunction(const double mu, const double sigma) :
        m_param({.mu = mu, .sigma = sigma}),
        m_inv_sigma_squared(1. / (sigma * sigma))
    {}

    WaveFunction(const WaveFunction& wf) :
        m_param(wf.param()),
        m_inv_sigma_squared(wf.get_inv_sigma_squared())
    {}

    WaveFunction(WaveFunction&& wf) :
        m_param(std::move(wf.m_param)),
        m_inv_sigma_squared(wf.m_inv_sigma_squared)
    {}

    WaveFunction& operator=(const PsiParam& par)
    {
        m_param = par;
        m_inv_sigma_squared = 1. / (par.sigma * par.sigma);
        return *this;
    }

    WaveFunction& operator=(PsiParam&& par)
    {
        m_param = std::move(par);
        m_inv_sigma_squared = 1. / (m_param.sigma * m_param.sigma);
        return *this;
    }

    WaveFunction& operator=(const WaveFunction& wf)
    {
        m_param = wf.param();
        m_inv_sigma_squared = wf.get_inv_sigma_squared();
        return *this;
    }

    WaveFunction& operator=(WaveFunction&& wf)
    {
        m_param = wf.m_param;
        m_inv_sigma_squared = wf.m_inv_sigma_squared;
        return *this;
    }

    double operator()(const double x) const;
    double second_der(const double x) const;
    double PDF(const double x) const;
    double normalized_PDF(const double x) const;
    double hamiltonian(const double x) const;
    double kinetic(const double x) const;
    double potential(const double x) const;

    // Getters
    constexpr double mu() const { return m_param.mu; }
    constexpr double sigma() const { return m_param.sigma; }
    constexpr double get_inv_sigma_squared() const { return m_inv_sigma_squared; }
    constexpr PsiParam param() const { return m_param; }

    // Setters
    void mu(const double mu) { m_param.mu = mu; }
    void sigma(const double sigma)
    {
        assert(sigma != 0. && "Sigma cannot be zero");
        m_param.sigma = sigma;
        m_inv_sigma_squared = 1. / (sigma * sigma);
    }
    void param(const PsiParam& par)
    {
        mu(par.mu);
        sigma(par.sigma);
    }
};

#endif // __WAVE_FUNCTION__
