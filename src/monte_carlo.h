#ifndef __MONTE_CARLO__
#define __MONTE_CARLO__

#include <functional>
#include <vector>

#include "utilities.h"
#include "random.h"

class MonteCarlo
{
private:
    unsigned int m_launch_size{100};
    unsigned int m_n_block{100};
    Random *m_rng{nullptr};

public:
    MonteCarlo() = default;
    MonteCarlo(const unsigned int launch_size, const unsigned int n_block) : m_launch_size{launch_size}, m_n_block{n_block} { ; }
    MonteCarlo(Random *rng) : m_rng{rng} { ; }

    inline void set_launch_size(const unsigned int launch) { m_launch_size = launch; }
    inline unsigned int get_launch() const { return m_launch_size; }

    inline void set_n_block(const unsigned int block) { m_n_block = block; }
    inline unsigned int get_n_block() const { return m_n_block; }

    inline void set_rng(Random *rng) { m_rng = rng; }

    values calculate_uniform(const double a, const double b, std::function<double(double)> integrand);
    values calculate_uniform(const double a, const double b, const unsigned int total_throws, std::function<double(double)> integrand);
    std::vector<values> calculate_uniform_step(const double a, const double b, const unsigned int block_size, unsigned int total_throws, std::function<double(double)> integrand);
    std::vector<values> calculate_distribution(const double a, const double b, const unsigned int block_size, unsigned int total_throws, std::function<double(double)> integrand, std::function<double(double)> inverse_cumulative);
    std::vector<values> calculate_integral_general(const double a, const double b, const unsigned int block_size, unsigned total_throws, std::function<double(double)> integrand, std::function<double(void)> ext_rng);
    // double calculate_cumulative(const double a, const double b, std::function<double(double)> integrand, std::function<double(double)> cumulative);
};

#endif // __MONTE_CARLO__
