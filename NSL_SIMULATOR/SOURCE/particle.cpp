/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "particle.h"
#include <iostream>
#include <math.h>

using namespace std;

void Particle ::initialize()
{
    _spin = 1;
    _x.resize(_ndim);
    _xold.resize(_ndim);
    _v.resize(_ndim);
    return;
}

void Particle ::translate(vec delta, vec side)
{
    for (unsigned int i = 0; i < _ndim; i++)
    {
        _x(i) = pbc(_x(i) + delta(i), side(i));
    }
}

void Particle ::setvelocity(int dim, double velocity)
{
    _v(dim) = velocity;
    return;
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
