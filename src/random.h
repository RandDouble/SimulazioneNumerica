/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>

#ifndef __Random_Sim__
#define __Random_Sim__

// namespace Sim
// {
class Random
{
public:
    using result_type = uint64_t;

private:
    static constexpr result_type m_multiplier = 34522712143931ull;
    result_type l_tot, n_tot;
    static constexpr result_type m_modulus = (1ull << 48);

public:
    // constructors
    Random() = default;
    // destructor
    ~Random() = default;

    // Method to initialize the RNG
    void Initializer(const std::string &seed_file = "seed.in", const std::string &prime_file = "Primes",
                     const std::size_t rows_to_skip = 0ull);
    // Method to set the seed for the RNG
    void SetRandom(int *, int, int);
    // Method to save the seed to a file
    void SaveSeed(const std::string &filename = "seed.out") const;
    // Method to generate a random number in the range [0,1)
    double Rannyu(void);
    // Method to generate a random number in the range [min,max)
    double Rannyu(const double min, const double max);
    // Method to generate random integers in range [0 , 2**48)
    result_type Ranint();
    // Method to generate random integers in range [min, max)]
    result_type Ranint(const uint64_t min, const uint64_t max);
    // Method to generate a random number with a Gaussian distribution
    double Gauss(const double mean, const double sigma);
    // Method to generate a random number with an Exponential distribution
    double Exponential(const double lambda);
    // Method to generate a random number with a Lorentian distribution
    double Lorenztian(const double x_0, const double gamma);
    // Method Accept Reject for extreme case
    double AcceptReject(const double a, const double b, const double max, std::function<double(double)> &PDF);
    double AcceptReject(const double a, const double b, const double max, const std::function<double(double)> &PDF);
    // Method with inverse cumulative
    double ExternalInvCum(std::function<double(double)> &ICDF);

    // uint64_t operator()(const uint64_t max) { return Ranint(0ull, max); }
    // uint64_t operator()(const uint64_t min, const uint64_t max) { return Ranint(min, max); }
    // double operator()(const double min, const double max) { return Rannyu(min, max); }

    result_type operator()()
    {
        return Ranint();
    }
    static constexpr result_type min()
    {
        return 0ull;
    }
    static constexpr result_type max()
    {
        return (m_modulus - 1);
    }
};

// }
#endif // __Random_Sim__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
