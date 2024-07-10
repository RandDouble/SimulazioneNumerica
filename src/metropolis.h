#ifndef __METROPOLIS__
#define __METROPOLIS__

#include <functional>
#include <numeric>

#include "random.h"

// #define _TEST_

class Metropolis
{
private:
    std::size_t m_n_step{0};
    Random m_rng;
    double m_acceptance{0.};

public:
    Metropolis() = default;

    Random* get_rng_ptr() { return &m_rng; }

    constexpr std::size_t get_n_step() { return m_n_step; }
    void set_n_step(std::size_t n_step) { m_n_step = n_step; }

    double get_acceptance() const { return m_acceptance; }

    template <typename T>
    T generate(T start_pos, const std::size_t n_step, const std::function<double(T)>& PDF, const std::function<T()>& sampler);

    template <typename T>
    T generate(T start_pos, const std::function<double(T)>& PDF, const std::function<T()>& sampler);

    template <typename T>
    double accept(const std::function<double(T)>& PDF, T next, T actual);
};

template <typename T>
inline T Metropolis::generate(T start_pos, const std::size_t n_step, const std::function<double(T)>& PDF, const std::function<T()>& sampler)
{
    T x = start_pos;
    m_acceptance = 0.;
    for (std::size_t i = 0; i < n_step; i++)
    {
        T x_new = x + sampler();
        double a = accept(PDF, x_new, x);

#ifdef _TEST_
        std::cout << "Acceptance : " << a << '\n';
#endif

        double test = m_rng.Rannyu();
        m_acceptance += (test <= a) / static_cast<double>(m_n_step);
        x = (test <= a) ? x_new : x;
    }
    return x;
}

template <typename T>
inline T Metropolis::generate(T start_pos, const std::function<double(T)>& PDF, const std::function<T()>& sampler)
{
    return generate(start_pos, m_n_step, PDF, sampler);
}

template <typename T>
inline double Metropolis::accept(const std::function<double(T)>& PDF, T next, T actual)
{
    return  PDF(next) / PDF(actual);
}

#endif // __METROPOLIS__