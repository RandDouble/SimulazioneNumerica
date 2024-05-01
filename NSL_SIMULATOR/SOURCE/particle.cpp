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

void Particle ::acceptmove()
{
    _xold = _x;
}

void Particle ::flip()
{
    _spin = -1 * this->getspin();
}

void Particle ::moveback()
{
    _x = _xold;
}

double Particle ::getposition(const int dim, const bool xnew) const
{
    return (xnew) ? _x(dim) : _xold(dim);
}

void Particle ::translate(const vec &delta, const vec &side)
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

double Particle ::pbc(const double position, const double side)
{
    return position - side * rint(position / side);
}

vec Particle ::getvelocity() const
{
    return _v;
}

double Particle ::getvelocity(const int dim) const
{
    return _v(dim);
}

void Particle ::setspin(int spin)
{
    _spin = spin;
    return;
}

void Particle ::setposition(int dim, double position)
{
    _x(dim) = position;
}

void Particle ::setpositold(int dim, double position)
{
    _xold(dim) = position;
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
