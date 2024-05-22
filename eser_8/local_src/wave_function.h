#include <cmath>
#include <numeric>

struct PsiParam
{
    double mu, sigma;
};

class WaveFunction
{
private:
    PsiParam m_param{.mu = 0., .sigma = 1.};

public:
    WaveFunction() = default;
    WaveFunction(const PsiParam par) : m_param{par} {}
    WaveFunction(const double mu, const double sigma) : m_param{.mu = mu, .sigma = sigma} {}

    double operator()(const double x) const;
    double second_der(const double x) const;
    double PDF(const double x) const;
    double hamiltonian(const double x) const;

    // Setters
    void mu(const double mu) { m_param.mu = mu; }
    void sigma(const double sigma) { m_param.sigma = sigma; }
    void param(const PsiParam par) { m_param = par; }

    // Getters
    double mu() const { return m_param.mu; }
    double sigma() const { return m_param.sigma; }
    PsiParam param() const { return m_param; }
};

// Old function...
double wave_function(const double x, const double mu, const double sigma);
double wave_function_prob(const double x, const double mu, const double sigma);
double wave_function_second(const double x, const double mu, const double sigma);

double potential(const double x);
double hamiltonian(const double x, const double mu = 0., const double sigma = 1.);

double boltzman_weight(const double beta, const double energy);
