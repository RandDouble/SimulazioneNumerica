/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random_Sim__
#define __Random_Sim__

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>

namespace Sim
{
    class Random
    {

    private:
        const uint64_t m_tot{34522712143931ull};
        uint64_t l_tot, n_tot;

    protected:
    public:
        // constructors
        Random() = default;
        // destructor
        ~Random() = default;

        // Method to set the seed for the RNG
        void SetRandom(int *, int, int);
        // Method to save the seed to a file
        void SaveSeed();
        // Method to generate a random number in the range [0,1)
        double Rannyu(void);
        // Method to generate a random number in the range [min,max)
        double Rannyu(const double min, const double max);
        // Method to generate a random number with a Gaussian distribution
        double Gauss(const double mean, const double sigma);
        // Method to generate a random number with an Exponential distribution
        double Exponential(const double lambda);
        // Method to generate a random number with a Lorentian distribution
        double Lorenztian(const double x_0, const double gamma);
        // Method Accept Reject for extreme case
        double AcceptReject(const double a, const double b, const double max, std::function<double(double)> &PDF);
        // Method with inverse cumulative
        double ExternalInvCum(std::function<double(double)> &ICDF);
    };

}
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