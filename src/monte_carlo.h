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

    /// @brief Calculate integral using average method, number of throws is default
    /// @param a inferior bound
    /// @param b superior bound
    /// @param integrand function to be integrated
    /// @return integral value
    double calculate_uniform(const double a, const double b, std::function<double(double)> integrand);
    
    /// @brief Calculate integral using average method
    /// @param a inferior buond
    /// @param b superior bound
    /// @param throws nr of extraction from random number generator 
    /// @param integrand function to be integrated
    /// @return integral value
    double calculate_uniform(const double a, const double b, const unsigned int total_throws, std::function<double(double)> integrand);


    /// @brief Return Block Means using average method uniform sampling, number of blocks is a setting in class
    /// @param a inferior bound
    /// @param b superior bound
    /// @param block_size nr of extraction per block
    /// @param integrand function to be integrated
    /// @return collection of value of integral with error
    std::vector<values> calculate_uniform_blocks(const double a, const double b, const unsigned int block_size, std::function<double(double)> integrand);
    
    double calculate_distribution_convert(const double a, const double b, const unsigned int throws, std::function<double(double)> integrand, std::function<double(double)> inverse_cumulative);
    double calculate_distribution_rng(const double a, const double b, const unsigned int throws, std::function<double(double)> integrand, std::function<double(void)> ext_rng);
    
    /// @brief Return Block Means using inverse cumulative sampling, number of blocks is a setting in class
    /// @param a inferior bound
    /// @param b superior bound
    /// @param block_size nr of extraction per block 
    /// @param integrand function to be integrated
    /// @param inverse_cumulative inverse cumulative supplied to random number generator to extract values
    /// @return collection of values of integral with error
    std::vector<values> calculate_dist_conv_block(const double a, const double b, const unsigned int block_size, std::function<double(double)> integrand, std::function<double(double)> inverse_cumulative);
    
    /// @brief Return Block Means using supplied random number generator, number of blocks is a setting in class
    /// @param a inferior bound
    /// @param b superior bpund
    /// @param block_size nr of extraction per block
    /// @param integrand function to be integrated
    /// @param ext_rng external random number generator for sampling
    /// @return collection of values of integral with error
    std::vector<values> calculate_dist_rng_block(const double a, const double b, const unsigned int block_size, std::function<double(double)> integrand, std::function<double(void)> ext_rng);
    // double calculate_cumulative(const double a, const double b, std::function<double(double)> integrand, std::function<double(double)> cumulative);
};

#endif // __MONTE_CARLO__
