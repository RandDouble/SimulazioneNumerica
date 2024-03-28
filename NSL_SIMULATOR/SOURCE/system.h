/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __System__
#define __System__

#include <armadillo>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdlib.h> //exit
#include <string>
#include <vector>
#include <execution>
#include <algorithm>

#include "particle.h"
#include "random.h"

using namespace std;
using namespace arma;


enum class SymType : int
{
    LENNARD_JONES_MD = 0,
    LENNARD_JONES_MC = 1,
    ISING_MRT2 = 2,
    GIBBS = 3,
};

std::istream &operator>>(std::istream &in, SymType &type);

struct measure_flags
{
    bool penergy{false}, kenergy{false}, tenergy{false}; // Flags for measuring different energies
    bool temp{false}, pressure{false}, gofr{false};      // Flags for measuring temperature, pressure, and radial dist. function
    bool magnet{false}, cv{false}, chi{false};           // Flags for measuring magnetization, heat capacity, and susceptibility
    int idx_penergy{0}, idx_kenergy{0}, idx_tenergy{0};  // Indices for accessing energy-related properties in vec _measurement
    int idx_temp{0}, idx_pressure{0}, idx_gofr{0};       // Indices for accessing temperature, pressure, and radial dist. function
    int idx_magnet{0}, idx_cv{0}, idx_chi{0};            // Indices for accessing magnetization, heat capacity, and susceptibility
    vector<stringstream> v_streams;                      // Vector containing stream to output to file
    vector<string> output_names;                         // Vector containing names of output files

    stringstream &stream_penergy() { return v_streams[idx_penergy]; } // Streams for outputting without opening 1000 times a file
    stringstream &stream_kenergy() { return v_streams[idx_kenergy]; }
    stringstream &stream_tenergy() { return v_streams[idx_tenergy]; }
    stringstream &stream_temp() { return v_streams[idx_temp]; }
    stringstream &stream_pressure() { return v_streams[idx_pressure]; }
    stringstream &stream_gofr() { return v_streams[idx_gofr]; }
    stringstream &stream_magnet() { return v_streams[idx_magnet]; }
    stringstream &stream_cv() { return v_streams[idx_cv]; }
    stringstream &stream_chi() { return v_streams[idx_chi]; }

    measure_flags() = default;

};

class System
{

private:
    const int _ndim = 3;       // Dimensionality of the system
    bool _restart;             // Flag indicating if the simulation is restarted
    SymType _sim_type;         // Type of simulation (e.g., Lennard-Jones, Ising)
    int _npart;                // Number of particles
    int _nblocks;              // Number of blocks for block averaging
    int _nsteps;               // Number of simulation steps in each block
    int _nattempts;            // Number of attempted moves
    int _naccepted;            // Number of accepted moves
    double _temp, _beta;       // Temperature and inverse temperature
    double _rho, _volume;      // Density and volume of the system
    double _r_cut;             // Cutoff radius for pair interactions
    double _delta;             // Displacement step for particle moves
    double _J, _H;             // Parameters for the Ising Hamiltonian
    vec _side;                 // Box dimensions
    vec _halfside;             // Half of box dimensions
    Random _rnd;               // Random number generator
    field<Particle> _particle; // Field of particle objects representing the system
    vec _fx, _fy, _fz;         // Forces on particles along x, y, and z directions

    // Properties
    int _nprop;             // Number of properties being measured
    measure_flags _measure; // Container for various flags and relative index to measure properties
    int _n_bins;            // Number of bins for radial distribution function
    double _bin_size;       // Size of bins for radial distribution function
    double _vtail, _ptail;  // Tail corrections for energy and pressure
    vec _block_av;          // Block averages of properties
    vec _global_av;         // Global averages of properties
    vec _global_av2;        // Squared global averages of properties
    vec _average;           // Average values of properties
    vec _measurement;       // Measured values of properties

public:                                                               // Function declarations
    constexpr int get_nbl() const { return _nblocks; }                // Get the number of blocks
    constexpr int get_nsteps() const { return _nsteps; }              // Get the number of steps in each block
    void initialize();                                                // Initialize system properties
    void initialize_properties();                                     // Initialize properties for measurement
    void finalize();                                                  // Finalize system and clean up
    void write_configuration() const;                                 // Write final system configuration to XYZ file
    void write_XYZ(const int nconf) const;                            // Write system configuration in XYZ format on the fly
    void write_velocities() const;                                    // Write final particle velocities to file
    void read_configuration();                                        // Read system configuration from file
    void initialize_velocities();                                     // Initialize particle velocities
    void step();                                                      // Perform a simulation step
    void block_reset(const int blk);                                  // Reset block averages
    void measure();                                                   // Measure properties of the system
    void averages(const int blk);                                     // Compute averages of properties
    double error(const double acc, const double acc2, const int blk); // Compute error
    void move(const int part);                                        // Move a particle
    bool metro(int part);                                             // Perform Metropolis acceptance-rejection step
    double pbc(double position, int i);                               // Apply periodic boundary conditions for coordinates
    int pbc(int i);                                                   // Apply periodic boundary conditions for spins
    void Verlet();                                                    // Perform Verlet integration step
    double Force(const int i, const int dim);                         // Calculate force on a particle along a dimension
    double Boltzmann(const int i, const bool xnew);                   // Calculate Boltzmann factor for Metropolis acceptance
};

#endif // __System__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
