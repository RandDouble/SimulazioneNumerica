/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Particle__
#define __Particle__

#include "random.h"
#include <armadillo>

using namespace std;
using namespace arma;

class Particle
{

private:
    const unsigned int _ndim = 3; // Dimensionality of the system
    int _spin;           // Spin of the particle (+1 or -1)
    vec _x;              // Current position vector
    vec _xold;           // Previous position vector (used in moveback())
    vec _v;              // Velocity vector

public:                                                // Function declarations
    void initialize();                                 // Initialize particle properties
    void translate(vec delta, vec side);               // Translate the particle within the simulation box
    void flip();                                // Flip the spin of the particle
    void moveback();                            // Move particle back to previous position
    void acceptmove();                          // Accept the proposed move and update particle properties
    constexpr int getspin() const;                           // Get the spin of the particle
    void setspin(int spin);                     // Set the spin of the particle
    constexpr double getposition(int dim, bool xnew) const;  // Get the position of the particle along a specific dimension
    void setposition(int dim, double position); // Set the position of the particle along a specific dimension
    void setpositold(int dim, double position); // Set the previous position of the particle along a specific dimension
    constexpr double getvelocity(int dim) const;             // Get the velocity of the particle along a specific dimension
    vec getvelocity() const;                       // Get the velocity vector of the particle
    void setvelocity(const int dim, const double velocity); // Set the velocity of the particle along a specific dimension
    double pbc(double position, double side);   // Apply periodic boundary conditions
};

inline void Particle ::flip()
{
    _spin = -1 * this->getspin();
}

inline void Particle ::moveback()
{
    _x = _xold;
}

inline void Particle ::acceptmove()
{
    _xold = _x;
}

constexpr int Particle ::getspin() const
{
    return _spin;
}

inline void Particle ::setspin(int spin)
{
    _spin = spin;
    return;
}

constexpr double Particle ::getposition(const int dim, const bool xnew) const
{
   return  (xnew)? _x(dim) : _xold(dim);
}

inline void Particle ::setposition(int dim, double position)
{
    _x(dim) = position;
    return;
}

inline void Particle ::setpositold(int dim, double position)
{
    _xold(dim) = position;
    return;
}

constexpr double Particle ::getvelocity(int dim) const
{
    return _v(dim);
}

inline double Particle ::pbc(double position, double side)
{
    return position - side * rint(position / side);
}

inline vec Particle ::getvelocity() const
{
    return _v;
}

#endif // __Particle__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
