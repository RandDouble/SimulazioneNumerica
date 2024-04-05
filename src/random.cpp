/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"

using namespace std;

Random ::Random() {}
// Default constructor, does not perform any action

Random ::~Random() {}
// Default destructor, does not perform any action

void Random ::SaveSeed()
{
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open())
   {
      WriteSeed << "RANDOMSEED	" << l1 << " " << l2 << " " << l3 << " " << l4 << endl;
      ;
   }
   else
      cerr << "PROBLEM: Unable to open random.out" << endl;
   WriteSeed.close();
   return;
}

double Random ::Gauss(const double mean, const double sigma)
{
   // This function generates a random number from a Gaussian distribution with given mean and sigma
   assert(sigma > 0 && "Sigma must be greater than zero\n");
   double s = Rannyu();
   double t = Rannyu();
   double x = sqrt(-2. * log(1. - s)) * cos(2. * M_PI * t);
   return mean + x * sigma;
}

double Random ::Rannyu(const double min, const double max)
{
   // This function generates a random number in the range [min, max)
   assert(max > min && "Max should be greather than min\n");
   return min + (max - min) * Rannyu();
}

double Random ::Rannyu(void)
{
   // This function generates a random number in the range [0,1)
   const double twom12 = 0.000244140625;
   int i1, i2, i3, i4;
   double r;

   i1 = l1 * m4 + l2 * m3 + l3 * m2 + l4 * m1 + n1;
   i2 = l2 * m4 + l3 * m3 + l4 * m2 + n2;
   i3 = l3 * m4 + l4 * m3 + n3;
   i4 = l4 * m4 + n4;
   l4 = i4 % 4096;
   i3 = i3 + i4 / 4096;
   l3 = i3 % 4096;
   i2 = i2 + i3 / 4096;
   l2 = i2 % 4096;
   l1 = (i1 + i2 / 4096) % 4096;
   r = twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * (l4))));

   return r;
}

void Random ::SetRandom(int *s, int p1, int p2)
{
   // This function sets the seed and parameters of the random number generator
   m1 = 502;
   m2 = 1521;
   m3 = 4071;
   m4 = 2107;
   l1 = s[0];
   l2 = s[1];
   l3 = s[2];
   l4 = s[3];
   n1 = 0;
   n2 = 0;
   n3 = p1;
   n4 = p2;

   return;
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

double Random::AcceptReject(const double a, const double b, const double max, std::function<double(double)> PDF)
{
   double x = 0, y = 0;

   do
   {
      x = Rannyu(a, b);
      y = Rannyu(0, max);
   } while (PDF(x) < y);

   return x;
}

double Random::ExternalInvCum(std::function<double(double)> ICDF)
{
   return ICDF(Rannyu());
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
