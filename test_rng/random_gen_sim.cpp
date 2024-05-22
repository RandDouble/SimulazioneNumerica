/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random_gen_sim.h"

namespace Sim
{

    Random ::Random() {}

    Random ::~Random() {}

    void Random ::SaveSeed()
    {
        std::ofstream WriteSeed;
        WriteSeed.open("../OUTPUT/seed.out");
        if (WriteSeed.is_open())
        {
            WriteSeed << ((l_tot >> 36) & 4095) << " " << ((l_tot >> 24) & 4095) << " "
                      << ((l_tot >> 12) & 4095) << " " << (l_tot & 4095) << std::endl;
        }
        else
            std::cerr << "PROBLEM: Unable to open random.out" << std::endl;
        WriteSeed.close();
        return;
    }

    double Random ::Gauss(const double mean, const double sigma)
    {
        double s = Rannyu();
        double t = Rannyu();
        double x = sqrt(-2. * std::log(1. - s)) * std::cos(2. * M_PI * t);
        return mean + x * sigma;
    }

    double Random ::Rannyu(const double min, const double max)
    {
        return min + (max - min) * Rannyu();
    }

    double Random ::Rannyu(void)
    {
        constexpr double twom48 = 1. / (1ull << 48);
        l_tot = l_tot * m_tot + n_tot;
        l_tot &= ((1ull << 48) - 1);
        double r = twom48 * l_tot;
        return r;
    }

    void Random ::SetRandom(int *s, int p1, int p2)
    {
        l_tot = (static_cast<uint64_t>(s[0]) << (12 * 3)) + (static_cast<uint64_t>(s[1]) << (12 * 2)) + (static_cast<uint64_t>(s[2]) << 12) + static_cast<uint64_t>(s[3]);

        uint16_t n1 = 0;
        uint16_t n2 = 0;
        uint16_t n3 = p1;
        uint16_t n4 = p2;

        n_tot = (static_cast<uint64_t>(n1) << (12 * 3)) + (static_cast<uint64_t>(n2) << (12 * 2)) + (static_cast<uint64_t>(n3) << (12 * 1)) + static_cast<uint64_t>(n4);
    }
    double Random::Exponential(const double lambda)
    {
        assert(lambda > 0 && "Lambda Parameter should be greater than zero\n");
        double y = Rannyu();
        double r = -std::log(1 - y) / lambda;
        return r;
    }

    double Random::Lorenztian(const double x_0, const double gamma)
    {
        assert(gamma > 0 && "Gamma should be greater than zero\n");
        double y = Rannyu();
        double r = gamma * std::tan(M_PI * (y - 0.5)) + x_0;
        return r;
    }

    double Random::AcceptReject(const double a, const double b, const double max, std::function<double(double)> &PDF)
    {
        double x = 0, y = 0;

        do
        {
            x = Rannyu(a, b);
            y = Rannyu(0, max);
        } while (PDF(x) < y);

        return x;
    }

    double Random::ExternalInvCum(std::function<double(double)> &ICDF)
    {
        return ICDF(Rannyu());
    }

} // namespace Sim

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
