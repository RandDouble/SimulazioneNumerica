/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "system.h"
#include <cmath>
#include <cstdlib>
#include <string>

#define NDEBUG_FORCE
#define NDEBUG_TEMPERATURE_PRESSURE
#define NDEBUG_PRESSURE

using namespace arma;

/// @brief Input stream operator overloading for SimType
/// @param in Input stream
/// @param type Output varible
/// @return Input stream
std::istream &operator>>(std::istream &in, SimType &type)
{
    int el; // Temporary variable to store symulation type, that in file is a number
    if (!(in >> el))
    {
        return in;
    }
    if (el > static_cast<int>(SimType::GIBBS))
    {
        in.setstate(in.rdstate() | std::ios::failbit);
        std::cerr << "PROBLEM: unknown simulation type" << std::endl;
        exit(EXIT_FAILURE);
    }

    type = static_cast<SimType>(el);
    return in;
}

/// @brief Perform a simulation step
void System::step()
{
    switch (_sim_type)
    {
    case SimType::LENNARD_JONES_MD:
        Verlet();
        break;

    default:
        for (int i = 0; i < _npart; i++)
        {
            this->move(int(_rnd.Rannyu() * _npart));
        }
        break;
    }

    _nattempts += _npart; // update number of attempts performed on the system
    return;
}

/// @brief Verlet Differential equation integrator
void System::Verlet()
{

#pragma omp parallel for ordered
    for (int i = 0; i < _npart; i++)
    { // Force acting on particle i
        _fx(i) = this->Force(i, 0);
        _fy(i) = this->Force(i, 1);
        _fz(i) = this->Force(i, 2);
    }

    // Change in position and velocity for each particle
    // Do not ask why, but this section when launched with openmp causes  to get only nans...
    // suspect of various entities reading from memory at the same time
#pragma omp barrier

#pragma omp parallel for
    for (int i = 0; i < _npart; i++)
    { // Verlet integration scheme
        double xnew = this->pbc(2.0 * _particle(i).getposition(0, true) - _particle(i).getposition(0, false) + _fx(i) * pow(_delta, 2), 0);
        double ynew = this->pbc(2.0 * _particle(i).getposition(1, true) - _particle(i).getposition(1, false) + _fy(i) * pow(_delta, 2), 1);
        double znew = this->pbc(2.0 * _particle(i).getposition(2, true) - _particle(i).getposition(2, false) + _fz(i) * pow(_delta, 2), 2);
        _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0, false), 0) / (2.0 * _delta));
        _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1, false), 1) / (2.0 * _delta));
        _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2, false), 2) / (2.0 * _delta));
        _particle(i).acceptmove(); // xold = xnew
        _particle(i).setposition(0, xnew);
        _particle(i).setposition(1, ynew);
        _particle(i).setposition(2, znew);
    }

    _naccepted += _npart;
    return;
}

/// @brief Calulate force acting on i-th particle in dim-th direction
/// @param i
/// @param dim
/// @return
double System::Force(const int i, const int dim)
{
    double f = 0.0;
    vec distance;
    distance.resize(_ndim);
#ifndef NDEBUG_FORCE
    std::vector<double> force_history;
    bool proced = false;
#endif
    for (int j = 0; j < _npart; j++)
    {
#ifndef NDEBUG_FORCE
        force_history.push_back(f);
#endif
        if (i != j)
        {
            distance(0) = this->pbc(_particle(i).getposition(0, true) - _particle(j).getposition(0, true), 0);
            distance(1) = this->pbc(_particle(i).getposition(1, true) - _particle(j).getposition(1, true), 1);
            distance(2) = this->pbc(_particle(i).getposition(2, true) - _particle(j).getposition(2, true), 2);
            double dr_square = dot(distance, distance);

            f += (dr_square < _r_cut_squared) * distance(dim) * (pow(dr_square, -7.) - 0.5 * pow(dr_square, -4.)); // Moved * 48. after for loop
                                                                                                                   // f += distance(dim) * (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8));

#ifndef NDEBUG_FORCE
            // std::cout << "Computed Force in " << dim << " direction : " << f << " on particle : " << i << '\n';
            if (std::isnan(f) and !proced)
            {
                std::cout << "current distance (dir : " << dim << " ):\n"
                          << distance << "current dr : " << dr << "\tcurrent particles (i,j) : " << i << '\t' << j << '\n'
                          << "Calculation result : pow(dr, -16.) : " << std::pow(dr, -16.) << "\t0.5 * pow(dr, -8.) : " << 0.5 * pow(dr, -8.) << '\n'
                          << "Increment Computed :" << distance(dim) * (pow(dr, -16.) - 0.5 * pow(dr, -8.)) << '\n'
                          << "Actual Force : " << f << "\tin direction " << dim << '\n';
                std::cin.get();

                std::cout << "Printing force history\n";

                for (auto &&past : force_history)
                {
                    std::cout << past << '\n';
                }

                auto res = std::cin.get();
                proced = (res == ' ');
            }
#endif // NDEBUG_FORCE
        }
    }
    f *= 48.;

    return f;
}

/// @brief Propose a MC move for i-th particle
/// @param i
void System ::move(const int i)
{
    switch (_sim_type)
    {
    case SimType::GIBBS:
    {
        // To be fixed in EXERCISE 6
        // 1. Choosing spin to change at random
        int idx_spin = static_cast<int>(std::floor(_rnd.Rannyu(0., _npart)));

        // 2. Compute energy of nearest spins
        int spin_sum = _particle(pbc(idx_spin - 1)).getspin() + _particle(pbc(idx_spin + 1)).getspin();
        double delta_E = _J * spin_sum + _H;
        // 3. Compute Change
        int new_spin = (_rnd.Rannyu() < (1. / (1. + std::exp(-2. * _beta * delta_E)))) ? 1 : -1;

        _particle(idx_spin).setspin(new_spin);

        _particle(idx_spin).acceptmove();
        _naccepted++;
    }
    break;

    case SimType::LENNARD_JONES_MC:
    {                     // M(RT)^2 LJ system
        vec shift(_ndim); // Will store the proposed translation
        for (int j = 0; j < _ndim; j++)
        {
            shift(j) = _rnd.Rannyu(-1.0, 1.0) * _delta; // uniform distribution in [-_delta;_delta)
        }
        _particle(i).translate(shift, _side); // Call the function Particle::translate
        if (this->metro(i))
        { // Metropolis acceptance evaluation
            _particle(i).acceptmove();
            _naccepted++;
        }
        else
            _particle(i).moveback(); // If translation is rejected, restore the old configuration
    }
    break;
    case SimType::LENNARD_JONES_MD:
        std::cerr << "Not Implemented Yet!!!\n";
        exit(-2);
        break;
    case SimType::ISING_MRT2:
    { // Ising 1D
        if (metro(i))
        {                        // Metropolis acceptance evaluation for a spin flip involving spin i
            _particle(i).flip(); // If accepted, the spin i is flipped
            _naccepted++;
        }
    }

    break;
    default:
        std::cout << "Came to selection with unknown choice\n";
        exit(-1);
        break;
    }
    return;
}

/// @brief Implements Metropolis Algorithm on u-th partcile
/// @param i
/// @return if step was accepted or not.
bool System::metro(const int i)
{
    double delta_E, acceptance;
    switch (_sim_type)
    {
    case SimType::LENNARD_JONES_MC:
        delta_E = this->Boltzmann(i, true) - this->Boltzmann(i, false);
        break;

    case SimType::ISING_MRT2:
        delta_E = 2.0 * _particle(i).getspin() *
                  (_J * (_particle(this->pbc(i - 1)).getspin() + _particle(this->pbc(i + 1)).getspin()) + _H);
        break;
    default:
        std::cerr << "Came in point of metro where you shouldn't be, Aborting\n";
        exit(-2);
    }

    acceptance = std::exp(-_beta * delta_E);
    // Usually acceptace is min(1, p(new) / p(old)), with exponential this is not necessary, with this extraction, min function is achieved by _rnd.Rannyu()
    bool decision = (_rnd.Rannyu() < acceptance); // Metropolis acceptance step

    return decision;
}

double System ::Boltzmann(const int i, const bool xnew)
{
    double energy_i = 0.0;
    // double dx, dy, dz, dr;
// Questo è un buon candidato per la parallizzazione... Inoltre credo che si possa sistemare un filino come codice.
#pragma omp parallel for
    for (int j = 0; j < _npart; j++)
    {
        if (j != i)
        {
            double dx = this->pbc(_particle(i).getposition(0, xnew) - _particle(j).getposition(0, 1), 0);
            double dy = this->pbc(_particle(i).getposition(1, xnew) - _particle(j).getposition(1, 1), 1);
            double dz = this->pbc(_particle(i).getposition(2, xnew) - _particle(j).getposition(2, 1), 2);
            double dr_squared = dx * dx + dy * dy + dz * dz;
            // dr =
            // dr = std::sqrt(dr);

            energy_i += (dr_squared < _r_cut_squared) * (pow(dr_squared, -6.) - pow(dr_squared, -3.)); // Potential energy calculation (pow(dr, -12.) - pow(dr, -6.))
        }
    }
    energy_i *= 4.0;
    return energy_i;
}

// @brief Print block information to stream
void System::general_print(std::ostream &stream, const int blk, const double ave, const double sum_ave, const double sum_ave2)
{
    stream.precision(8);
    stream << std::setw(8) << blk
           << std::setw(14) << ave
           << std::setw(14) << sum_ave / double(blk)
           << std::setw(14) << this->error(sum_ave, sum_ave2, blk) << "\n";
}

void System::general_print(std::ostream &stream, const double blk, const double ave, const double sum_ave, const double sum_ave2)
{
    stream.precision(8);
    stream << std::setw(8) << blk
           << std::setw(14) << ave
           << std::setw(14) << sum_ave / double(blk)
           << std::setw(14) << this->error(sum_ave, sum_ave2, blk) << "\n";
}

/// @brief Enforce periodic boundary conditions
/// @param position
/// @param i
/// @return Returns position after Boudary Condition are applied
double System::pbc(double position, int i)
{
    return position - _side(i) * std::rint(position / _side(i));
}

/// @brief Enforce periodic boundary conditions for spins
/// @param i
/// @return Return Spins after proper Boundary Condition are applied
int System ::pbc(int i)
{
    if (i >= _npart)
        i = i - _npart;
    else if (i < 0)
        i = i + _npart;
    return i;
}

/// @brief Initialize the System object according to the content of the input files in the ../INPUT/ directory
void System ::initialize()
{

    int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
    std::ifstream Primes("../INPUT/Primes");
    Primes >> p1 >> p2;
    Primes.close();
    int seed[4]; // Read the seed of the RNG
    std::ifstream Seed("../INPUT/seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    _rnd.SetRandom(seed, p1, p2);

    std::ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
    couta << "#   N_BLOCK:  ACCEPTANCE:"
          << "\n";
    couta.close();

    std::ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
    std::ofstream coutf;
    coutf.open("../OUTPUT/output.dat");
    std::string property;
    double delta;
    while (!input.eof())
    {
        input >> property;
        if (property == "SIMULATION_TYPE")
        {
            input >> _sim_type;
            if (_sim_type > SimType::LENNARD_JONES_MC) // Ising M(RT)^2 or GIBBS
            {
                input >> _J;
                input >> _H;
            }

            switch (_sim_type)
            {
            case SimType::LENNARD_JONES_MD:
                coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"
                      << "\n";
                break;
            case SimType::LENNARD_JONES_MC:
                coutf << "LJ MONTE CARLO (NVT) SIMULATION"
                      << "\n";
                break;
            case SimType::ISING_MRT2:
                coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION"
                      << '\n'
                      << "SIM_TYPE=" << std::setw(4) << static_cast<int>(SimType::ISING_MRT2)
                      << std::setw(6) << _J
                      << std::setw(6) << _H
                      << '\n';
                break;
            case SimType::GIBBS:
                coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION"
                      << "\n"
                      << "SIM_TYPE=" << std::setw(4) << static_cast<int>(SimType::GIBBS)
                      << std::setw(6) << _J
                      << std::setw(6) << _H
                      << '\n';
                break;
            }
        }
        else if (property == "RESTART")
        {
            input >> _restart;
        }
        else if (property == "TEMP")
        {
            input >> _temp;
            _beta = 1.0 / _temp;
            coutf << "TEMPERATURE= " << _temp << "\n";
        }
        else if (property == "NPART")
        {
            input >> _npart;
            _fx.resize(_npart);
            _fy.resize(_npart);
            _fz.resize(_npart);
            _particle.set_size(_npart);
            for (int i = 0; i < _npart; i++)
            {
                _particle(i).initialize();
                if (_rnd.Rannyu() > 0.5)
                    _particle(i).flip(); // to randomize the spin configuration
            }
            coutf << "NPART= " << _npart << "\n";
        }
        else if (property == "RHO")
        {
            input >> _rho;
            _volume = _npart / _rho;
            _side.resize(_ndim);
            _halfside.resize(_ndim);
            double side = pow(_volume, 1.0 / 3.0);
            for (int i = 0; i < _ndim; i++)
            {
                _side(i) = side;
            }

            _halfside = 0.5 * _side;

            coutf << "SIDE= ";
            for (int i = 0; i < _ndim; i++)
            {
                coutf << std::setw(12) << _side[i];
            }
            coutf << "\n";
        }
        else if (property == "R_CUT")
        {
            input >> _r_cut;
            coutf << "R_CUT= " << _r_cut << "\n";
            _r_cut_squared = _r_cut * _r_cut;
        }
        else if (property == "DELTA")
        {
            input >> delta;
            coutf << "DELTA= " << delta << "\n";
            _delta = delta;
        }
        else if (property == "NBLOCKS")
        {
            input >> _nblocks;
            coutf << "NBLOCKS= " << _nblocks << "\n";
        }
        else if (property == "NSTEPS")
        {
            input >> _nsteps;
            coutf << "NSTEPS= " << _nsteps << "\n";
        }
        else if (property == "ENDINPUT")
        {
            coutf << "Reading input completed!"
                  << "\n";
            break;
        }
        else
            std::cerr << "PROBLEM: unknown input : " << property << std::endl;
    }
    input.close();
    this->read_configuration();
    this->initialize_velocities();
    coutf << "System initialized!"
          << "\n";
    coutf.close();
    return;
}

/// @brief Inizialize velocity using config file `velocities.in`
void System::initialize_velocities()
{
    if (_restart and _sim_type == SimType::LENNARD_JONES_MD)
    {
        std::ifstream cinf;
        cinf.open("../INPUT/CONFIG/velocities.in");
        if (cinf.is_open())
        {
            double vx, vy, vz;
            for (int i = 0; i < _npart; i++)
            {
                cinf >> vx >> vy >> vz;
                _particle(i).setvelocity(0, vx);
                _particle(i).setvelocity(1, vy);
                _particle(i).setvelocity(2, vz);
            }
        }
        else
            cerr << "PROBLEM: Unable to open INPUT file velocities.in" << endl;
        cinf.close();
    }
    else
    {
        vec vx(_npart), vy(_npart), vz(_npart);
        vec sumv(_ndim);
        sumv.zeros();
        for (int i = 0; i < _npart; i++)
        {
            vx(i) = _rnd.Gauss(0., std::sqrt(_temp));
            vy(i) = _rnd.Gauss(0., std::sqrt(_temp));
            vz(i) = _rnd.Gauss(0., std::sqrt(_temp));
            sumv(0) += vx(i);
            sumv(1) += vy(i);
            sumv(2) += vz(i);
        }

        for (int idim = 0; idim < _ndim; idim++)
        {
            sumv(idim) = sumv(idim) / double(_npart);
        }

        double sumv2 = 0.0, scalef;
        for (int i = 0; i < _npart; i++)
        {
            vx(i) = vx(i) - sumv(0);
            vy(i) = vy(i) - sumv(1);
            vz(i) = vz(i) - sumv(2);
            sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
        }

        sumv2 /= double(_npart);
        scalef = std::sqrt(3.0 * _temp / sumv2); // velocity scale factor
        for (int i = 0; i < _npart; i++)
        {
            _particle(i).setvelocity(0, vx(i) * scalef);
            _particle(i).setvelocity(1, vy(i) * scalef);
            _particle(i).setvelocity(2, vz(i) * scalef);
        }
    }

    if (_sim_type == SimType::LENNARD_JONES_MD)
    {
        for (int i = 0; i < _npart; i++)
        {
            double xold = this->pbc(_particle(i).getposition(0, true) - _particle(i).getvelocity(0) * _delta, 0);
            double yold = this->pbc(_particle(i).getposition(1, true) - _particle(i).getvelocity(1) * _delta, 1);
            double zold = this->pbc(_particle(i).getposition(2, true) - _particle(i).getvelocity(2) * _delta, 2);
            _particle(i).setpositold(0, xold);
            _particle(i).setpositold(1, yold);
            _particle(i).setpositold(2, zold);
        }
    }
    return;
}

/// @brief Initialize data members used for measurement of properties using config written in `properties.dat`. Also prepares files for output.
void System ::initialize_properties()
{
    std::string property;
    int index_property = 0;
    _nprop = 0;

    std::ifstream input("../INPUT/properties.dat");
    if (input.is_open())
    {
        while (!input.eof())
        {
            input >> property;
            if (property == "POTENTIAL_ENERGY")
            {
                std::ofstream coutp("../OUTPUT/potential_energy.dat");
                coutp << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:"
                      << "\n";
                coutp.close();
                _nprop++;
                _measure.penergy = true;
                _measure.idx_penergy = index_property;
                _measure.v_streams.emplace_back(std::stringstream(std::ios::out | std::ios::app)); // This will simplify some operations later
                _measure.output_names.emplace_back("../OUTPUT/potential_energy.dat");              // This will simplify some operations later
                index_property++;
                _vtail = 8. * M_PI * _rho * (1. - 3 * std::pow(_r_cut, 6.)) / (9 * std::pow(_r_cut, 9.)); // TO BE FIXED IN EXERCISE 7
            }
            else if (property == "KINETIC_ENERGY")
            {
                std::ofstream coutk("../OUTPUT/kinetic_energy.dat");
                coutk << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:"
                      << "\n";
                coutk.close();
                _nprop++;
                _measure.kenergy = true;
                _measure.idx_kenergy = index_property;
                _measure.v_streams.emplace_back(std::stringstream(std::ios::out | std::ios::app)); // This will simplify some operations later
                _measure.output_names.emplace_back("../OUTPUT/kinetic_energy.dat");                // This will simplify some operations later
                index_property++;
            }
            else if (property == "TOTAL_ENERGY")
            {
                std::ofstream coutt("../OUTPUT/total_energy.dat");
                coutt << "#     BLOCK:   ACTUAL_TE:    TE_AVE:      ERROR:"
                      << "\n";
                coutt.close();
                _nprop++;
                _measure.tenergy = true;
                _measure.idx_tenergy = index_property;
                _measure.v_streams.emplace_back(std::stringstream(std::ios::out | std::ios::app));
                _measure.output_names.emplace_back("../OUTPUT/total_energy.dat");
                index_property++;
            }
            else if (property == "TEMPERATURE")
            {
                std::ofstream coutte("../OUTPUT/temperature.dat");
                coutte << "#     BLOCK:   ACTUAL_T:     T_AVE:       ERROR:"
                       << "\n";
                coutte.close();
                _nprop++;
                _measure.temp = true;
                _measure.idx_temp = index_property;
                _measure.v_streams.emplace_back(std::stringstream(std::ios::out | std::ios::app));
                _measure.output_names.emplace_back("../OUTPUT/temperature.dat");
                index_property++;
            }
            else if (property == "PRESSURE")
            {
                std::ofstream coutpr("../OUTPUT/pressure.dat");
                coutpr << "#     BLOCK:   ACTUAL_P:     P_AVE:       ERROR:"
                       << "\n";
                coutpr.close();
                _nprop++;
                _measure.pressure = true;
                _measure.idx_pressure = index_property;
                _measure.v_streams.emplace_back(std::stringstream(std::ios::out | std::ios::app));
                _measure.output_names.emplace_back("../OUTPUT/pressure.dat");
                index_property++;
                _ptail = 16 * M_PI * _rho * (2. - 3. * std::pow(_r_cut, 6.)) / (9. * std::pow(_r_cut, 9.)); // TO BE FIXED IN EXERCISE 7
            }
            else if (property == "GOFR")
            {
                std::ofstream coutgr("../OUTPUT/gofr.dat");
                coutgr << "# DISTANCE:     AVE_GOFR:        ERROR:"
                       << "\n";
                coutgr.close();
                input >> _n_bins;
                _nprop += _n_bins;
                _bin_size = (_halfside.min()) / (double)_n_bins;
                _measure.gofr = true;
                _measure.idx_gofr = index_property;
                _measure.v_streams.emplace_back(std::stringstream(std::ios::out | std::ios::app));
                _measure.output_names.emplace_back("../OUTPUT/gofr.dat");
                index_property += _n_bins;
            }
            else if (property == "MAGNETIZATION")
            {
                std::ofstream coutpr("../OUTPUT/magnetization.dat");
                coutpr << "#     BLOCK:   ACTUAL_M:     M_AVE:       ERROR:"
                       << "\n";
                coutpr.close();
                _nprop++;
                _measure.magnet = true;
                _measure.idx_magnet = index_property;
                _measure.v_streams.emplace_back(std::stringstream(std::ios::out | std::ios::app));
                _measure.output_names.emplace_back("../OUTPUT/magnetization.dat");
                index_property++;
            }
            else if (property == "SPECIFIC_HEAT")
            {
                std::ofstream coutpr("../OUTPUT/specific_heat.dat");
                coutpr << "#     BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:"
                       << "\n";
                coutpr.close();
                _nprop++;
                _measure.cv = true;
                _measure.idx_cv = index_property;
                _measure.v_streams.emplace_back(std::stringstream(std::ios::out | std::ios::app));
                _measure.output_names.emplace_back("../OUTPUT/specific_heat.dat");
                index_property++;
            }
            else if (property == "SUSCEPTIBILITY")
            {
                std::ofstream coutpr("../OUTPUT/susceptibility.dat");
                coutpr << "#     BLOCK:   ACTUAL_X:     X_AVE:       ERROR:"
                       << "\n";
                coutpr.close();
                _nprop++;
                _measure.chi = true;
                _measure.idx_chi = index_property;
                _measure.v_streams.emplace_back(std::stringstream(std::ios::out | std::ios::app));
                _measure.output_names.emplace_back("../OUTPUT/susceptibility.dat");
                index_property++;
            }
            else if (property == "ENDPROPERTIES")
            {
                std::ofstream coutf;
                coutf.open("../OUTPUT/output.dat", ios::app);
                coutf << "Reading properties completed!"
                      << "\n";
                coutf.close();
                break;
            }
            else
                cerr << "PROBLEM: unknown property" << endl;
        }
        input.close();
    }
    else
        cerr << "PROBLEM: Unable to open properties.dat" << endl;

    // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
    _measurement.resize(_nprop);
    _average.resize(_nprop);
    _block_av.resize(_nprop);
    _global_av.resize(_nprop);
    _global_av2.resize(_nprop);
    _average.zeros();
    _global_av.zeros();
    _global_av2.zeros();
    _nattempts = 0;
    _naccepted = 0;
    return;
}

/// @brief Write final config of the system to file, write all the measure requested to files, writes seed used to file.
void System ::finalize()
{
    this->write_configuration(); // write to file final config of the system
    _rnd.SaveSeed();
    std::ofstream coutf;

    for (size_t i = 0; i < _measure.v_streams.size(); i++) // Write to files each measure result
    {
        coutf.open(_measure.output_names[i], std::ios::app);
        coutf << _measure.v_streams[i].str();
        coutf.close();
    }

    coutf.open("../OUTPUT/output.dat", std::ios::app);
    coutf << "Simulation completed!"
          << "\n";
    coutf.close();
    return;
}

/// @brief Write current configuration as file in directory ../OUTPUT/CONFIG/, if LJ simulation or Montecarlo saves as .xyz file, else if is Ising writes a .spin file
void System ::write_configuration() const
{
    std::ofstream coutf;
    if (_sim_type < SimType::ISING_MRT2) // Select Lennard Jones MD or MONTECARLO
    {
        coutf.open("../OUTPUT/CONFIG/config.xyz");
        if (coutf.is_open())
        {
            coutf << _npart << "\n";
            coutf << "#Comment!"
                  << "\n";
            for (int i = 0; i < _npart; i++)
            {
                coutf << "LJ"
                      << "  "
                      << std::setw(16) << _particle(i).getposition(0, true) / _side(0)          // x
                      << std::setw(16) << _particle(i).getposition(1, true) / _side(1)          // y
                      << std::setw(16) << _particle(i).getposition(2, true) / _side(2) << "\n"; // z
            }
        }
        else
            cerr << "PROBLEM: Unable to open config.xyz" << endl;
        coutf.close();
        this->write_velocities();
    }
    else // Ising Metropolis or Ising with Gibbs
    {
        coutf.open("../OUTPUT/CONFIG/config.spin");
        for (int i = 0; i < _npart; i++)
            coutf << _particle(i).getspin() << " ";
        coutf.close();
    }
    return;
}

/// @brief Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
/// @param nconf
void System ::write_XYZ(const int nconf) const
{
    std::ofstream coutf;
    coutf.open("../OUTPUT/CONFIG/config_" + std::to_string(nconf) + ".xyz");
    if (coutf.is_open())
    {
        coutf << _npart << "\n";
        coutf << "#Comment!\n";
        for (int i = 0; i < _npart; i++)
        {
            coutf << "LJ"
                  << "  "
                  << std::setw(16) << _particle(i).getposition(0, true)          // x
                  << std::setw(16) << _particle(i).getposition(1, true)          // y
                  << std::setw(16) << _particle(i).getposition(2, true) << "\n"; // z
        }
    }
    else
        cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    return;
}

/// @brief Write to file particle velocities.
void System::write_velocities() const
{
    std::ofstream coutf;
    coutf.open("../OUTPUT/CONFIG/velocities.out");
    if (coutf.is_open())
    {
        for (int i = 0; i < _npart; i++)
        {
            coutf << std::setw(16) << _particle(i).getvelocity(0)          // vx
                  << std::setw(16) << _particle(i).getvelocity(1)          // vy
                  << std::setw(16) << _particle(i).getvelocity(2) << "\n"; // vz
        }
    }
    else
    {
        cerr << "PROBLEM: Unable to open velocities.dat" << endl;
    }
    coutf.close();
    return;
}

/// @brief Read configuration from a .xyz or .spin file in directory ../INPUT/CONFIG/
void System ::read_configuration()
{
    std::ifstream cinf;
    cinf.open("../INPUT/CONFIG/config.xyz");
    if (cinf.is_open())
    {
        std::string comment;
        std::string particle;
        double x, y, z;
        int ncoord;
        cinf >> ncoord;
        if (ncoord != _npart)
        {
            cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!"
                 << "\n";
            exit(EXIT_FAILURE);
        }
        cinf >> comment;

        for (int i = 0; i < _npart; i++)
        {
            cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
            _particle(i).setposition(0, this->pbc(_side(0) * x, 0));
            _particle(i).setposition(1, this->pbc(_side(1) * y, 1));
            _particle(i).setposition(2, this->pbc(_side(2) * z, 2));
            _particle(i).acceptmove(); // _x_old = _x_new
        }
    }
    else
        cerr << "PROBLEM: Unable to open INPUT file config.xyz"
             << "\n";
    cinf.close();
    if (_restart and _sim_type > SimType::LENNARD_JONES_MC)
    {
        int spin;
        cinf.open("../INPUT/CONFIG/config.spin");
        for (int i = 0; i < _npart; i++)
        {
            cinf >> spin;
            _particle(i).setspin(spin);
        }
        cinf.close();
    }
    return;
}

/// @brief Reset block accumulators to zero
/// @param blk
void System::block_reset(int blk)
{
    std::ofstream coutf;
    if (blk > 0)
    {
        coutf.open("../OUTPUT/output.dat", ios::app);
        coutf << "Block completed: " << blk << "\n";
        coutf.close();
    }
    _block_av.zeros();
    return;
}

/// @brief Measure properties if relative flag is checked.
void System::measure()
{
    _measurement.zeros();
    // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
    // int bin;
    vec distance;
    distance.resize(_ndim);
    double penergy_temp = 0.0; // temporary accumulator for potential energy
    double kenergy_temp = 0.0; // temporary accumulator for kinetic energy
    double tenergy_temp = 0.0; // temporary accumulator for total energy
    double magnetization = 0.0;
    double virial = 0.0;

    // VIRIAL  ////////////////////////////////////////////////////////////////////
    if (_measure.penergy or _measure.pressure or _measure.gofr)
    {
        std::vector index(_npart - 1, 0); // I want N_part - 1 elements,
        std::iota(index.begin(), index.end(), 0);

        // #pragma omp parallel for ordered reduction(+ : virial)
        //         for (int analyzed = 0; analyzed < _npart - 1; analyzed++)
        //         {
        //             for (int other = analyzed + 1; other < _npart; other++)
        //             {

        //                 distance(0) = this->pbc(_particle(analyzed).getposition(0, true) - _particle(other).getposition(0, true), 0);
        //                 distance(1) = this->pbc(_particle(analyzed).getposition(1, true) - _particle(other).getposition(1, true), 1);
        //                 distance(2) = this->pbc(_particle(analyzed).getposition(2, true) - _particle(other).getposition(2, true), 2);
        //                 dr = std::sqrt(arma::dot(distance, distance));
        //                 // GOFR ... TO BE FIXED IN EXERCISE 7
        //                 if (_measure.gofr)
        //                 {
        //                     // Ragionando sulla condizione ho trovato un modo per trovare l'indice ed evitare quindi il ciclo for
        //                     int index_to_insert_gofr = static_cast<int>(std::floor(dr / _bin_size)); // Voglio essere sicuro che faccia un troncamento verso il basso
        //                     _measurement(_measure.idx_gofr + index_to_insert_gofr) += 2;
        //                     /*
        //                     for (int i = 0; i < _n_bins; i++)
        //                     {
        //                         if (dr > i* _bin_size && dr < (i+1) * _bin_size)
        //                             {
        //                             _measurement(_measure.idx_gofr + i) += 2;
        //                         }
        //                     }
        //                     */
        //                 }

        //                 // POTENTIAL ENERGY
        //                 penergy_temp += (dr < _r_cut and _measure.penergy) * (std::pow(dr, -12.) - std::pow(dr, -6.)); // POTENTIAL ENERGY

        //                 // VIRIAL FOR PRESSURE ... TO BE FIXED IN EXERCISE 4

        //                 virial += (dr < _r_cut and _measure.pressure) * (std::pow(dr, -12.) - 0.5 * std::pow(dr, -6.)); // VIRIAL, multiplication by 48 done after
        //             }
        //         }

        std::for_each(std::execution::par, index.cbegin(), index.cend(), [&](const int &part_analyzed)
                      {
        for (int other_part = part_analyzed + 1; other_part < _npart; other_part++)
        {
            distance(0) = this->pbc(_particle(part_analyzed).getposition(0, true) - _particle(other_part).getposition(0, true), 0);
            distance(1) = this->pbc(_particle(part_analyzed).getposition(1, true) - _particle(other_part).getposition(1, true), 1);
            distance(2) = this->pbc(_particle(part_analyzed).getposition(2, true) - _particle(other_part).getposition(2, true), 2);
            double dr_squared = dot(distance, distance);
            // GOFR ... TO BE FIXED IN EXERCISE 7
            if (_measure.gofr)
            {
                // Ragionando sulla condizione ho trovato un modo per trovare l'indice ed evitare quindi il ciclo for
                int index_to_insert_gofr = static_cast<int>(std::sqrt(dr_squared / (_bin_size * _bin_size))); // Voglio essere sicuro che faccia un troncamento verso il basso
                _measurement(_measure.idx_gofr + index_to_insert_gofr) += 2;
            }

            penergy_temp += (dr_squared < _r_cut_squared and _measure.penergy) *( pow(dr_squared, -6.) - pow(dr_squared, -3.)); // POTENTIAL ENERGY ( 1.0 / pow(dr, 12.) - 1.0 / pow(dr, 6.))

            // PRESSURE ... TO BE FIXED IN EXERCISE 4
            virial +=  (dr_squared < _r_cut_squared and _measure.pressure) * (std::pow(dr_squared, -6.) - 0.5 * std::pow(dr_squared, -3.)); // VIRIAL, multiplication by 48 done after, std::pow(dr, -12.) - 0.5 * std::pow(dr, -6.)


        } });
    }

    // POTENTIAL ENERGY //////////////////////////////////////////////////////////
    if (_measure.penergy)
    {
        penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
        _measurement(_measure.idx_penergy) = penergy_temp;
    }

    // KINETIC ENERGY ////////////////////////////////////////////////////////////
    if (_measure.kenergy)
    {
        for (int i = 0; i < _npart; i++)
        {
            kenergy_temp += 0.5 * dot(_particle(i).getvelocity(), _particle(i).getvelocity());
        }
        kenergy_temp /= double(_npart);
        _measurement(_measure.idx_kenergy) = kenergy_temp;
    }

    // TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
    if (_measure.tenergy)
    {
        switch (_sim_type)
        {
        case SimType::LENNARD_JONES_MC:
        case SimType::LENNARD_JONES_MD:
            _measurement(_measure.idx_tenergy) = kenergy_temp + penergy_temp;
            break;
        case SimType::ISING_MRT2:
        case SimType::GIBBS:
            double s_i, s_j;
            for (int i = 0; i < _npart; i++)
            {
                s_i = double(_particle(i).getspin());
                s_j = double(_particle(this->pbc(i + 1)).getspin());
                tenergy_temp += -_J * s_i * s_j - 0.5 * _H * (s_i + s_j);
            }
            tenergy_temp /= double(_npart);
            _measurement(_measure.idx_tenergy) = tenergy_temp;
            break;
        }
    }
    // TEMPERATURE ///////////////////////////////////////////////////////////////
    if (_measure.temp and _measure.kenergy)
        _measurement(_measure.idx_temp) = (2.0 / 3.0) * kenergy_temp;

    // PRESSURE //////////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 4
    if (_measure.pressure and _measure.temp)
    {
        double temperature = _measurement(_measure.idx_temp);
        _measurement(_measure.idx_pressure) = _ptail + _rho * temperature + 16. * virial / (_volume); // 48. / 3. = 16...
#ifndef NDEBUG_TEMPERATURE_PRESSURE
        std::cout << "current virial value" << std::setw(8) << virial << '\n';
        std::cout << "Actuale Temperature : " << temperature << "\tPressure : " << _measurement(_measure.idx_pressure) << "\n";
#endif // NDEBUG_TEMPERATURE_PRESSURE
    }

    // MAGNETIZATION /////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6

    if (_measure.magnet)
    {

        for (int i = 0; i < _npart; i++)
        {
            double s_i = double(_particle(i).getspin());
            magnetization += s_i;
        }
        magnetization /= static_cast<double>(_npart);
        _measurement(_measure.idx_magnet) = magnetization; // this function as it is
    }

    // SPECIFIC HEAT //////////////////////////////////////////o///////////////////
    // TO BE FIXED IN EXERCISE 6
    if (_measure.cv)
    {
        // Saving total energy squared. The Mean calculation will happen in System::averages.
        double tenergy_squared = (tenergy_temp * _npart) * (tenergy_temp * _npart);
        _measurement(_measure.idx_cv) = tenergy_squared;
    }

    // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6
    if (_measure.chi)
    {
        // double chi_temp = magnetization * magnetization * _beta / static_cast<double>(_npart); // this Somehow functions... lets try calcultion explicitly

        double temp_chi = magnetization * magnetization * _npart;
        // _measurement(_measure.idx_chi) = chi_temp;
        _measurement(_measure.idx_chi) = temp_chi;
    }

    // After all ifs

    _block_av += _measurement; // Update block accumulators

    return;
}

/// @brief Perform averages in block `blk` for each measure obtained
/// @param blk Block index
void System::averages(const int blk)
{

    std::ofstream coutf;
    double average, sum_average, sum_ave2;

    _average = _block_av / double(_nsteps);
    _global_av += _average;
    _global_av2 += _average % _average; // % -> element-wise multiplication

    // POTENTIAL ENERGY //////////////////////////////////////////////////////////
    if (_measure.penergy)
    {
        average = _average(_measure.idx_penergy);
        sum_average = _global_av(_measure.idx_penergy);
        sum_ave2 = _global_av2(_measure.idx_penergy);
        general_print(_measure.stream_penergy(), blk, average, sum_average, sum_ave2);
    }
    // KINETIC ENERGY ////////////////////////////////////////////////////////////
    if (_measure.kenergy)
    {
        average = _average(_measure.idx_kenergy);
        sum_average = _global_av(_measure.idx_kenergy);
        sum_ave2 = _global_av2(_measure.idx_kenergy);
        general_print(_measure.stream_kenergy(), blk, average, sum_average, sum_ave2);
    }
    // TOTAL ENERGY //////////////////////////////////////////////////////////////
    if (_measure.tenergy)
    {
        average = _average(_measure.idx_tenergy);
        sum_average = _global_av(_measure.idx_tenergy);
        sum_ave2 = _global_av2(_measure.idx_tenergy);
        general_print(_measure.stream_tenergy(), blk, average, sum_average, sum_ave2);
    }
    // TEMPERATURE ///////////////////////////////////////////////////////////////
    if (_measure.temp)
    {
        average = _average(_measure.idx_temp);
        sum_average = _global_av(_measure.idx_temp);
        sum_ave2 = _global_av2(_measure.idx_temp);
        general_print(_measure.stream_temp(), blk, average, sum_average, sum_ave2);
    }
    // PRESSURE //////////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 4
    if (_measure.pressure)
    {
#ifndef NDEBUG_PRESSURE
        std::cout << "Pressione Media Blocco : " << _average(_measure.idx_pressure) << "\n";
#endif
        average = _average(_measure.idx_pressure);
        sum_average = _global_av(_measure.idx_pressure);
        sum_ave2 = _global_av2(_measure.idx_pressure);
        general_print(_measure.stream_pressure(), blk, average, sum_average, sum_ave2);
    }

    // GOFR //////////////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 7
    if (_measure.gofr)
    {

        if (blk == get_nbl()) // Last cycle in this moment we write gofr config.
        {
            _measure.stream_gofr().precision(12);
            for (int bin = 0; bin < _n_bins; bin++)
            {
                average = _average(_measure.idx_gofr + bin);
                sum_average = _global_av(_measure.idx_gofr + bin);
                sum_ave2 = _global_av2(_measure.idx_gofr + bin);
                general_print(_measure.stream_gofr(), blk, average, sum_average, sum_ave2);
            }
        }
    }
    // MAGNETIZATION /////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6
    if (_measure.magnet)
    {
        average = _average(_measure.idx_magnet);
        sum_average = _global_av(_measure.idx_magnet);
        sum_ave2 = _global_av2(_measure.idx_magnet);
        general_print(_measure.stream_magnet(), blk, average, sum_average, sum_ave2);
    }

    // SPECIFIC HEAT /////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6

    if (_measure.cv)
    {
        auto beta_squared = _beta * _beta;
        auto tenergy_renorm_squared = _average(_measure.idx_tenergy) * _npart;
        tenergy_renorm_squared *= tenergy_renorm_squared;

        average = beta_squared * (_average(_measure.idx_cv) - tenergy_renorm_squared); //
        average /= static_cast<double>(_npart);

        // Resetting last value and recomputing with actuale average
        _global_av(_measure.idx_cv) -= _average(_measure.idx_cv);
        _global_av2(_measure.idx_cv) -= _average(_measure.idx_cv) * _average(_measure.idx_cv);

        _global_av(_measure.idx_cv) += average;
        _global_av2(_measure.idx_cv) += average * average;

        // Output
        sum_average = _global_av(_measure.idx_cv);
        sum_ave2 = _global_av2(_measure.idx_cv);
        general_print(_measure.stream_cv(), blk, average, sum_average, sum_ave2);
    }

    // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6

    if (_measure.chi)
    {
        double mean_magnetization_squared = _average(_measure.idx_magnet) * _average(_measure.idx_magnet);
        average = _beta * (_average(_measure.idx_chi) - mean_magnetization_squared);

        // Resetting last value and recomputing with actuale average
        _global_av(_measure.idx_chi) -= _average(_measure.idx_chi);
        _global_av2(_measure.idx_chi) -= _average(_measure.idx_chi) * _average(_measure.idx_chi);

        _global_av(_measure.idx_chi) += average;
        _global_av2(_measure.idx_chi) += average * average;

        sum_average = _global_av(_measure.idx_chi);
        sum_ave2 = _global_av2(_measure.idx_chi);
        general_print(_measure.stream_chi(), blk, average, sum_average, sum_ave2);
    }

    // ACCEPTANCE ////////////////////////////////////////////////////////////////
    double fraction;
    coutf.open("../OUTPUT/acceptance.dat", ios::app);
    if (_nattempts > 0)
        fraction = double(_naccepted) / double(_nattempts);
    else
        fraction = 0.0;
    coutf << std::setw(12) << blk << std::setw(12) << fraction << "\n";
    coutf.close();

    return;
}

/// @brief Perform calculation of block error as \f$ \frac{<val^2> - <val>^2}{N} \f$.
/// @param acc First accumulator variable
/// @param acc2 Second Accumulator variable, the accumulator in this case is the sum of the squared values
/// @param blk number of elements in the block
/// @return Computed error
double System::error(const double acc, const double acc2, const int blk)
{
    return (blk <= 1) ? 0.0 : sqrt(fabs(acc2 / double(blk) - pow(acc / double(blk), 2)) / double(blk));
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
