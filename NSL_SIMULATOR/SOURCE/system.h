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

#define ARMA_DONT_PRINT_FAST_MATH_WARNING
#include <algorithm>
#include <armadillo>
#include <execution>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdlib.h> //exit
#include <string>
#include <vector>

#include "particle.h"
#include "random.h"

using namespace arma;

enum class SimType : int
{
    LENNARD_JONES_MD = 0,
    LENNARD_JONES_MC = 1,
    ISING_MRT2 = 2,
    GIBBS = 3,
};

std::istream& operator>>(std::istream& in, SimType& type);

struct measure_flags
{
    bool penergy{false}, kenergy{false}, tenergy{false};         // Flags for measuring different energies
    bool temp{false}, pressure{false}, gofr{false};              // Flags for measuring temperature, pressure, and radial dist. function
    bool magnet{false}, cv{false}, chi{false};                   // Flags for measuring magnetization, heat capacity, and susceptibility
    unsigned int idx_penergy{0}, idx_kenergy{0}, idx_tenergy{0}; // Indices for accessing energy-related properties in vec _measurement
    unsigned int idx_temp{0}, idx_pressure{0}, idx_gofr{0};      // Indices for accessing temperature, pressure, and radial dist. function
    unsigned int idx_magnet{0}, idx_cv{0}, idx_chi{0};           // Indices for accessing magnetization, heat capacity, and susceptibility
    std::vector<std::stringstream> v_streams;                    // Vector containing stream to output to file
    std::vector<std::string> output_names;                       // Vector containing names of output files

    std::stringstream& stream_penergy() { return v_streams[idx_penergy]; } // Streams for outputting without opening 1000 times a file
    std::stringstream& stream_kenergy() { return v_streams[idx_kenergy]; }
    std::stringstream& stream_tenergy() { return v_streams[idx_tenergy]; }
    std::stringstream& stream_temp() { return v_streams[idx_temp]; }
    std::stringstream& stream_pressure() { return v_streams[idx_pressure]; }
    std::stringstream& stream_gofr() { return v_streams[idx_gofr]; }
    std::stringstream& stream_partial_gofr() { return v_streams[idx_gofr + 1]; }
    std::stringstream& stream_magnet() { return v_streams[idx_magnet]; }
    std::stringstream& stream_cv() { return v_streams[idx_cv]; }
    std::stringstream& stream_chi() { return v_streams[idx_chi]; }

    measure_flags() = default;
};

class System
{

private:
    const int _ndim{3};                 // Dimensionality of the system
    bool _restart{false};               // Flag indicating if the simulation is restarted
    SimType _sim_type;                  // Type of simulation (e.g., Lennard-Jones, Ising)
    unsigned int _npart;                // Number of particles
    unsigned int _nblocks;              // Number of blocks for block averaging
    unsigned int _nsteps;               // Number of simulation steps in each block
    unsigned int _nattempts{0};         // Number of attempted moves
    unsigned int _naccepted{0};         // Number of accepted moves
    std::size_t _seed_line{0};          // Line to skip in the seed file
    double _temp, _beta;                // Temperature and inverse temperature
    double _rho, _volume;               // Density and volume of the system
    double _r_cut;                      // Cutoff radius for pair interactions
    double _r_cut_squared;              // Square of cutoff radius for pair interactions
    double _r_gofr_cut_squared;         // Square of cutoff radius for radial distribution function
    double _delta;                      // Displacement step for particle moves
    double _J{0.}, _H{0.};              // Parameters for the Ising Hamiltonian
    arma::vec3 _side{arma::fill::ones}; // Box dimensions
    arma::vec3 _halfside{_side / 2.};   // Half of box dimensions
    Random _rnd;                        // Random number generator
    field<Particle> _particle;          // Field of particle objects representing the system
    vec _fx, _fy, _fz;                  // Forces on particles along x, y, and z directions

    // Properties
    unsigned int _nprop{0};        // Number of properties being measured
    measure_flags _measure;        // Container for various flags and relative index to measure properties
    unsigned int _n_bins{0};       // Number of bins for radial distribution function
    double _bin_size{0.};          // Size of bins for radial distribution function
    double _vtail{0.}, _ptail{0.}; // Tail corrections for energy and pressure
    vec _block_av;                 // Block averages of properties
    vec _global_av;                // Global averages of properties
    vec _global_av2;               // Squared global averages of properties
    vec _average;                  // Average values of properties
    vec _measurement;              // Measured values of properties

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
    bool metro(const unsigned int part);                              // Perform Metropolis acceptance-rejection step

    double pbc(const double position, const unsigned int i); // Apply periodic boundary conditions for coordinates
    int pbc(const int i);                           // Apply periodic boundary conditions for spins
    arma::vec3 pbc(const arma::vec3& position);              // Apply periodic boundary conditions for vectors
    arma::vec3 particle_distance(const arma::vec3& first, const arma::vec3& second);
    double particle_distance_squared(const arma::vec3& first, const arma::vec3& second);

    void Verlet();                                                  // Perform Verlet integration step
    arma::vec3 Force(const unsigned int i);                         // Calculate force on a particle along a dimension
    double Boltzmann(const unsigned int i /* , const bool xnew */); // Calculate Boltzmann factor for Metropolis acceptance

    void general_print(std::ostream& stream, const int blk, const double ave, const double sum_ave, const double sum_ave2);
    void general_print(std::ostream& stream, const double blk, const double ave, const double sum_ave, const double sum_ave2);
    void final_gofr_print(std::ostream& stream, const unsigned int blk, const double bin, const double sum_ave, const double sum_ave2);
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
