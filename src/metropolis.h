#ifndef __METROPOLIS__
#define __METROPOLIS__

#include <functional>
#include <numeric>

#include "random.h"

// #define _TEST_

class Metropolis
{
protected:
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
    T generate(T start_pos,
               const std::size_t n_step,
               const std::function<double(T)>& PDF,
               const std::function<T()>& sampler);

    template <typename T>
    T generate(T start_pos,
               const std::function<double(T)>& PDF,
               const std::function<T()>& sampler);

    template <typename T>
    double accept(const std::function<double(T)>& PDF, T next, T actual);
};

template <typename T>
class MetropolisMemory : public Metropolis
{
private:
    T m_old_position;

public:
    MetropolisMemory(T initial_pos) : Metropolis(), m_old_position(initial_pos) {}

    void set_current_position(T pos) { m_old_position = pos; }
    T get_current_position() const { return m_old_position; }

    T generate(T start_pos,
               const std::function<double(T)>& PDF,
               const std::function<T()>& sampler)
    {
        m_old_position = start_pos;
        m_acceptance = 0.;
        for (std::size_t i = 0; i < m_n_step; i++)
        {
            T x_new = m_old_position + sampler();
            const double a = accept(PDF, x_new, m_old_position);
            const bool test = (m_rng.Rannyu() <= a);
            m_acceptance += static_cast<double>(test) / static_cast<double>(m_n_step);
            m_old_position = (test) ? x_new : m_old_position;
        }
        return m_old_position;
    }

    T generate(const std::function<double(T)>& PDF, const std::function<T()>& sampler)
    {
        return MetropolisMemory::generate(m_old_position, PDF, sampler);
    }
};

template <typename T>
inline T Metropolis::generate(T start_pos,
                              const std::size_t n_step,
                              const std::function<double(T)>& PDF,
                              const std::function<T()>& sampler)
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
inline T Metropolis::generate(T start_pos,
                              const std::function<double(T)>& PDF,
                              const std::function<T()>& sampler)
{
    return generate(start_pos, m_n_step, PDF, sampler);
}

template <typename T>
inline double Metropolis::accept(const std::function<double(T)>& PDF, T next, T actual)
{
    const double p_next = PDF(next);
    const double p_actual = PDF(actual);
    assert(!std::isnan(p_next) && "PDF(next) is NaN");
    assert(!std::isnan(p_actual) && "PDF(actual) is NaN");
    assert(!std::isinf(p_next) && "PDF(next) is Infinite");
    assert(!std::isinf(p_actual) && "PDF(actual) is Infinite");
    return p_next / p_actual;
}

#endif // __METROPOLIS__
