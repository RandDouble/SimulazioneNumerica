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

using namespace arma;

class Particle
{

private:
    const unsigned int _ndim = 3; // Dimensionality of the system
    int _spin;                    // Spin of the particle (+1 or -1)
    arma::vec3 _x;                // Current position vector
    arma::vec3 _xold;             // Previous position vector (used in move_back())
    arma::vec3 _v;                // Velocity vector

public:
    Particle() : _spin{1}, _x(arma::fill::zeros), _xold(arma::fill::zeros), _v(arma::fill::zeros){};

    // Function declarations
    void initialize() { _spin = 1; };                                                                                       // Initialize particle properties
    void translate(const arma::vec3& delta, const arma::vec3& side) { _x = pbc(_x + delta, side); };                        // Translate the particle within the simulation box
    void flip() { _spin *= -1; }                                                                                            // Flip the spin of the particle
    void move_back() { _x = _xold; }                                                                                        // Move particle back to previous position
    void accept_move() { _xold = _x; };                                                                                     // Accept the proposed move and update particle properties
    double pbc(const double position, const double side) { return position - side * rint(position / side); }                // Apply periodic boundary conditions
    arma::vec3 pbc(const arma::vec3& position, const arma::vec3& side) { return position - side % floor(position / side); } // Apply periodic boundary conditions

    constexpr int get_spin() const { return _spin; } // Get the spin of the particle
    void set_spin(int spin) { _spin = spin; }        // Set the spin of the particle

    double get_position(int dim, bool xnew) const { return (xnew) ? _x(dim) : _xold(dim); } // Get the position of the particle along a specific dimension
    arma::vec3 get_position(bool xnew) const { return (xnew) ? _x : _xold; }                // Get the position vector of the particle along all dimensions
    void set_position(int dim, double position) { _x(dim) = position; }                     // Set the position of the particle along a specific dimension
    void set_position_old(int dim, double position) { _xold(dim) = position; }              // Set the previous position of the particle along a specific dimension
    void set_position(arma::vec3& position) { _x = position; }                              // Set the position vector of the particle
    void set_position_old(arma::vec3& position) { _xold = position; }                       // Set the previous position vector of the particle

    double get_velocity(const int dim) const { return _v(dim); }                    // Get the velocity of the particle along a specific dimension
    arma::vec3 get_velocity() const { return _v; }                                  // Get the velocity vector of the particle
    void set_velocity(const int dim, const double velocity) { _v(dim) = velocity; } // Set the velocity of the particle along a specific dimension
    void set_velocity(const arma::vec3& velocity) { _v = velocity; }                // Set the velocity vector of the particle

    unsigned int dimension() const { return _ndim; }; // Get the dimensionality of the system
};

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
