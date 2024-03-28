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

using namespace std;
using namespace arma;

std::istream &operator>>(std::istream &in, SymType &type)
{
    int el;
    if (!(in >> el))
    {
        return in;
    }
    if (el > static_cast<int>(SymType::GIBBS))
    {
        in.setstate(in.rdstate() | std::ios::failbit);
        cerr << "PROBLEM: unknown simulation type" << endl;
        exit(EXIT_FAILURE);
    }

    type = static_cast<SymType>(el);
    return in;
}

void System ::step()
{
    // Perform a simulation step
    if (_sim_type == SymType::LENNARD_JONES_MD) // Perform a MD step
    {
        this->Verlet();
    }
    else
    {
        for (int i = 0; i < _npart; i++)
        {
            this->move(int(_rnd.Rannyu() * _npart));
        }
    }                     // Perform a MC step on a randomly choosen particle
    _nattempts += _npart; // update number of attempts performed on the system
    return;
}

void System ::Verlet()
{
    double xnew, ynew, znew;
    for (int i = 0; i < _npart; i++)
    { // Force acting on particle i
        _fx(i) = this->Force(i, 0);
        _fy(i) = this->Force(i, 1);
        _fz(i) = this->Force(i, 2);
    }
    for (int i = 0; i < _npart; i++)
    { // Verlet integration scheme
        xnew = this->pbc(2.0 * _particle(i).getposition(0, true) - _particle(i).getposition(0, false) + _fx(i) * pow(_delta, 2), 0);
        ynew = this->pbc(2.0 * _particle(i).getposition(1, true) - _particle(i).getposition(1, false) + _fy(i) * pow(_delta, 2), 1);
        znew = this->pbc(2.0 * _particle(i).getposition(2, true) - _particle(i).getposition(2, false) + _fz(i) * pow(_delta, 2), 2);
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

double System::Force(const int i, const int dim)
{
    double f = 0.0, dr;
    vec distance;
    distance.resize(_ndim);
    for (int j = 0; j < _npart; j++)
    {
        if (i != j)
        {
            distance(0) = this->pbc(_particle(i).getposition(0, true) - _particle(j).getposition(0, true), 0);
            distance(1) = this->pbc(_particle(i).getposition(1, true) - _particle(j).getposition(1, true), 1);
            distance(2) = this->pbc(_particle(i).getposition(2, true) - _particle(j).getposition(2, true), 2);
            dr = sqrt(dot(distance, distance));
            if (dr < _r_cut)
            {
                f += distance(dim) * (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8));
            }
        }
    }
    return f;
}

void System ::move(const int i)
{
    // Propose a MC move for particle i
    switch (_sim_type)
    {
    case SymType::GIBBS:
        // To be fixed in esercise 6
        break;
    case SymType::LENNARD_JONES_MC:
        // M(RT)^2
        if (_sim_type == SymType::LENNARD_JONES_MC) // LJ system
        {
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
    case SymType::LENNARD_JONES_MD:
        break;
    case SymType::ISING_MRT2:
    { // Ising 1D
        if (this->metro(i))
        {                        // Metropolis acceptance evaluation for a spin flip involving spin i
            _particle(i).flip(); // If accepted, the spin i is flipped
            _naccepted++;
        }
    }

    break;
    default:
        break;
    }
    return;
}

bool System ::metro(const int i)
{ // Metropolis algorithm
    bool decision = false;
    double delta_E, acceptance;
    if (_sim_type == SymType::LENNARD_JONES_MC)
        delta_E = this->Boltzmann(i, true) - this->Boltzmann(i, false);
    else
        delta_E = 2.0 * _particle(i).getspin() *
                  (_J * (_particle(this->pbc(i - 1)).getspin() + _particle(this->pbc(i + 1)).getspin()) + _H);
    acceptance = exp(-_beta * delta_E);
    if (_rnd.Rannyu() < acceptance)
        decision = true; // Metropolis acceptance step
    return decision;
}

double System ::Boltzmann(const int i, const bool xnew)
{
    double energy_i = 0.0;
    double dx, dy, dz, dr;
    for (int j = 0; j < _npart; j++)
    {
        if (j != i)
        {
            dx = this->pbc(_particle(i).getposition(0, xnew) - _particle(j).getposition(0, 1), 0);
            dy = this->pbc(_particle(i).getposition(1, xnew) - _particle(j).getposition(1, 1), 1);
            dz = this->pbc(_particle(i).getposition(2, xnew) - _particle(j).getposition(2, 1), 2);
            dr = dx * dx + dy * dy + dz * dz;
            dr = sqrt(dr);
            if (dr < _r_cut)
            {
                energy_i += 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6);
            }
        }
    }
    return 4.0 * energy_i;
}

/// @brief Enforce periodic boundary conditions
/// @param position
/// @param i
/// @return Reuturns position after Boudary Condition are applied
double System::pbc(double position, int i)
{
    return position - _side(i) * rint(position / _side(i));
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
    ifstream Primes("../INPUT/Primes");
    Primes >> p1 >> p2;
    Primes.close();
    int seed[4]; // Read the seed of the RNG
    ifstream Seed("../INPUT/seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    _rnd.SetRandom(seed, p1, p2);

    ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
    couta << "#   N_BLOCK:  ACCEPTANCE:"
          << "\n";
    couta.close();

    ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
    ofstream coutf;
    coutf.open("../OUTPUT/output.dat");
    string property;
    double delta;
    while (!input.eof())
    {
        input >> property;
        if (property == "SIMULATION_TYPE")
        {
            input >> _sim_type;
            if (_sim_type > SymType::LENNARD_JONES_MC)
            {
                input >> _J;
                input >> _H;
            }
            switch (_sim_type)
            {
            case SymType::LENNARD_JONES_MD:
                coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"
                      << "\n";
                break;
            case SymType::LENNARD_JONES_MC:
                coutf << "LJ MONTE CARLO (NVT) SIMULATION"
                      << "\n";
                break;
            case SymType::ISING_MRT2:
                coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION"
                      << "\n";
                break;
            case SymType::GIBBS:
                coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION"
                      << "\n";
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
                coutf << setw(12) << _side[i];
            }
            coutf << "\n";
        }
        else if (property == "R_CUT")
        {
            input >> _r_cut;
            coutf << "R_CUT= " << _r_cut << "\n";
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
            cerr << "PROBLEM: unknown input" << endl;
    }
    input.close();
    this->read_configuration();
    this->initialize_velocities();
    coutf << "System initialized!"
          << "\n";
    coutf.close();
    return;
}

void System ::initialize_velocities()
{
    if (_restart and _sim_type == SymType::LENNARD_JONES_MD)
    {
        ifstream cinf;
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
            vx(i) = _rnd.Gauss(0., sqrt(_temp));
            vy(i) = _rnd.Gauss(0., sqrt(_temp));
            vz(i) = _rnd.Gauss(0., sqrt(_temp));
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
        scalef = sqrt(3.0 * _temp / sumv2); // velocity scale factor
        for (int i = 0; i < _npart; i++)
        {
            _particle(i).setvelocity(0, vx(i) * scalef);
            _particle(i).setvelocity(1, vy(i) * scalef);
            _particle(i).setvelocity(2, vz(i) * scalef);
        }
    }

    if (_sim_type == SymType::LENNARD_JONES_MD)
    {
        double xold, yold, zold;
        for (int i = 0; i < _npart; i++)
        {
            xold = this->pbc(_particle(i).getposition(0, true) - _particle(i).getvelocity(0) * _delta, 0);
            yold = this->pbc(_particle(i).getposition(1, true) - _particle(i).getvelocity(1) * _delta, 1);
            zold = this->pbc(_particle(i).getposition(2, true) - _particle(i).getvelocity(2) * _delta, 2);
            _particle(i).setpositold(0, xold);
            _particle(i).setpositold(1, yold);
            _particle(i).setpositold(2, zold);
        }
    }
    return;
}

/// @brief Initialize data members used for measurement of properties
void System ::initialize_properties()
{
    string property;
    int index_property = 0;
    _nprop = 0;

    ifstream input("../INPUT/properties.dat");
    if (input.is_open())
    {
        while (!input.eof())
        {
            input >> property;
            if (property == "POTENTIAL_ENERGY")
            {
                ofstream coutp("../OUTPUT/potential_energy.dat");
                coutp << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:"
                      << "\n";
                coutp.close();
                _nprop++;
                _measure.penergy = true;
                _measure.idx_penergy = index_property;
                _measure.v_streams.emplace_back(stringstream(ios::out | ios::app));   // This will simplify some operations later
                _measure.output_names.emplace_back("../OUTPUT/potential_energy.dat"); // This will simplify some operations later
                index_property++;
                _vtail = 0.0; // TO BE FIXED IN EXERCISE 7
            }
            else if (property == "KINETIC_ENERGY")
            {
                ofstream coutk("../OUTPUT/kinetic_energy.dat");
                coutk << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:"
                      << "\n";
                coutk.close();
                _nprop++;
                _measure.kenergy = true;
                _measure.idx_kenergy = index_property;
                _measure.v_streams.emplace_back(stringstream(ios::out | ios::app)); // This will simplify some operations later
                _measure.output_names.emplace_back("../OUTPUT/kinetic_energy.dat"); // This will simplify some operations later
                index_property++;
            }
            else if (property == "TOTAL_ENERGY")
            {
                ofstream coutt("../OUTPUT/total_energy.dat");
                coutt << "#     BLOCK:   ACTUAL_TE:    TE_AVE:      ERROR:"
                      << "\n";
                coutt.close();
                _nprop++;
                _measure.tenergy = true;
                _measure.idx_tenergy = index_property;
                _measure.v_streams.emplace_back(stringstream(ios::out | ios::app));
                _measure.output_names.emplace_back("../OUTPUT/total_energy.dat");
                index_property++;
            }
            else if (property == "TEMPERATURE")
            {
                ofstream coutte("../OUTPUT/temperature.dat");
                coutte << "#     BLOCK:   ACTUAL_T:     T_AVE:       ERROR:"
                       << "\n";
                coutte.close();
                _nprop++;
                _measure.temp = true;
                _measure.idx_temp = index_property;
                _measure.v_streams.emplace_back(stringstream(ios::out | ios::app));
                _measure.output_names.emplace_back("../OUTPUT/temperature.dat");
                index_property++;
            }
            else if (property == "PRESSURE")
            {
                ofstream coutpr("../OUTPUT/pressure.dat");
                coutpr << "#     BLOCK:   ACTUAL_P:     P_AVE:       ERROR:"
                       << "\n";
                coutpr.close();
                _nprop++;
                _measure.pressure = true;
                _measure.idx_pressure = index_property;
                _measure.v_streams.emplace_back(stringstream(ios::out | ios::app));
                _measure.output_names.emplace_back("../OUTPUT/pressure.dat");
                index_property++;
                _ptail = 0.0; // TO BE FIXED IN EXERCISE 7
            }
            else if (property == "GOFR")
            {
                ofstream coutgr("../OUTPUT/gofr.dat");
                coutgr << "# DISTANCE:     AVE_GOFR:        ERROR:"
                       << "\n";
                coutgr.close();
                input >> _n_bins;
                _nprop += _n_bins;
                _bin_size = (_halfside.min()) / (double)_n_bins;
                _measure.gofr = true;
                _measure.idx_gofr = index_property;
                _measure.v_streams.emplace_back(stringstream(ios::out | ios::app));
                _measure.output_names.emplace_back("../OUTPUT/gofr.dat");
                index_property += _n_bins;
            }
            else if (property == "MAGNETIZATION")
            {
                ofstream coutpr("../OUTPUT/magnetization.dat");
                coutpr << "#     BLOCK:   ACTUAL_M:     M_AVE:       ERROR:"
                       << "\n";
                coutpr.close();
                _nprop++;
                _measure.magnet = true;
                _measure.idx_magnet = index_property;
                _measure.v_streams.emplace_back(stringstream(ios::out | ios::app));
                _measure.output_names.emplace_back("../OUTPUT/magnetization.dat");
                index_property++;
            }
            else if (property == "SPECIFIC_HEAT")
            {
                ofstream coutpr("../OUTPUT/specific_heat.dat");
                coutpr << "#     BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:"
                       << "\n";
                coutpr.close();
                _nprop++;
                _measure.cv = true;
                _measure.idx_cv = index_property;
                _measure.v_streams.emplace_back(stringstream(ios::out | ios::app));
                _measure.output_names.emplace_back("../OUTPUT/specific_heat.dat");
                index_property++;
            }
            else if (property == "SUSCEPTIBILITY")
            {
                ofstream coutpr("../OUTPUT/susceptibility.dat");
                coutpr << "#     BLOCK:   ACTUAL_X:     X_AVE:       ERROR:"
                       << "\n";
                coutpr.close();
                _nprop++;
                _measure.chi = true;
                _measure.idx_chi = index_property;
                _measure.v_streams.emplace_back(stringstream(ios::out | ios::app));
                _measure.output_names.emplace_back("../OUTPUT/susceptibility.dat");
                index_property++;
            }
            else if (property == "ENDPROPERTIES")
            {
                ofstream coutf;
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

void System ::finalize()
{
    this->write_configuration();
    _rnd.SaveSeed();
    ofstream coutf;

    for (size_t i = 0; i < _measure.v_streams.size(); i++)
    {
        coutf.open(_measure.output_names[i], ios::app);
        coutf << _measure.v_streams[i].str();
        coutf.close();
    }
    coutf.open("../OUTPUT/output.dat", ios::app);
    coutf << "Simulation completed!"
          << "\n";
    coutf.close();
    return;
}

/// @brief Write current configuration as a .xyz file in directory ../OUTPUT/CONFIG/
void System ::write_configuration() const
{
    ofstream coutf;
    if (_sim_type < SymType::ISING_MRT2)
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
                      << setw(16) << _particle(i).getposition(0, true) / _side(0)          // x
                      << setw(16) << _particle(i).getposition(1, true) / _side(1)          // y
                      << setw(16) << _particle(i).getposition(2, true) / _side(2) << "\n"; // z
            }
        }
        else
            cerr << "PROBLEM: Unable to open config.xyz" << endl;
        coutf.close();
        this->write_velocities();
    }
    else
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
    ofstream coutf;
    coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");
    if (coutf.is_open())
    {
        coutf << _npart << "\n";
        coutf << "#Comment!\n";
        for (int i = 0; i < _npart; i++)
        {
            coutf << "LJ"
                  << "  "
                  << setw(16) << _particle(i).getposition(0, true)          // x
                  << setw(16) << _particle(i).getposition(1, true)          // y
                  << setw(16) << _particle(i).getposition(2, true) << "\n"; // z
        }
    }
    else
        cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    return;
}

void System ::write_velocities() const
{
    ofstream coutf;
    coutf.open("../OUTPUT/CONFIG/velocities.out");
    if (coutf.is_open())
    {
        for (int i = 0; i < _npart; i++)
        {
            coutf << setw(16) << _particle(i).getvelocity(0)          // vx
                  << setw(16) << _particle(i).getvelocity(1)          // vy
                  << setw(16) << _particle(i).getvelocity(2) << "\n"; // vz
        }
    }
    else
    {
        cerr << "PROBLEM: Unable to open velocities.dat" << endl;
    }
    coutf.close();
    return;
}

/// @brief Read configuration from a .xyz file in directory ../OUTPUT/CONFIG/
void System ::read_configuration()
{
    ifstream cinf;
    cinf.open("../INPUT/CONFIG/config.xyz");
    if (cinf.is_open())
    {
        string comment;
        string particle;
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
    if (_restart and _sim_type > SymType::LENNARD_JONES_MC)
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
    ofstream coutf;
    if (blk > 0)
    {
        coutf.open("../OUTPUT/output.dat", ios::app);
        coutf << "Block completed: " << blk << "\n";
        coutf.close();
    }
    _block_av.zeros();
    return;
}

/// @brief Measure properties
void System::measure()
{
    _measurement.zeros();
    // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
    int bin;
    vec distance;
    distance.resize(_ndim);
    double penergy_temp = 0.0, dr = 0.0; // temporary accumulator for potential energy
    double kenergy_temp = 0.0;           // temporary accumulator for kinetic energy
    double tenergy_temp = 0.0;           // temporary accumulator for total energy
    double magnetization = 0.0;
    double virial = 0.0;

    if (_measure.penergy or _measure.pressure or _measure.gofr)
    {
        for (int i = 0; i < _npart - 1; i++)
        {
            for (int j = i + 1; j < _npart; j++)
            {
                distance(0) = this->pbc(_particle(i).getposition(0, true) - _particle(j).getposition(0, true), 0);
                distance(1) = this->pbc(_particle(i).getposition(1, true) - _particle(j).getposition(1, true), 1);
                distance(2) = this->pbc(_particle(i).getposition(2, true) - _particle(j).getposition(2, true), 2);
                dr = sqrt(dot(distance, distance));
                // GOFR ... TO BE FIXED IN EXERCISE 7
                if (dr < _r_cut and _measure.penergy)
                {
                    penergy_temp += 1.0 / pow(dr, 12) - 1.0 / pow(dr, 6); // POTENTIAL ENERGY
                }
                if (dr < _r_cut and _measure.pressure) // PRESSURE ... TO BE FIXED IN EXERCISE 4
                {
                    virial += 1.0 / pow(dr, 12) - 0.5 / pow(dr, 6);
                }
            }
        }
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
        if (_sim_type < SymType::ISING_MRT2)
        {
            _measurement(_measure.idx_tenergy) = kenergy_temp + penergy_temp;
        }
        else
        {
            double s_i, s_j;
            for (int i = 0; i < _npart; i++)
            {
                s_i = double(_particle(i).getspin());
                s_j = double(_particle(this->pbc(i + 1)).getspin());
                tenergy_temp += -_J * s_i * s_j - 0.5 * _H * (s_i + s_j);
            }
            tenergy_temp /= double(_npart);
            _measurement(_measure.idx_tenergy) = tenergy_temp;
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
        _measurement(_measure.idx_pressure) = _rho * temperature + 16. * virial / (_volume); // 48 / 3 = 16...
    }

    // MAGNETIZATION /////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6
    // SPECIFIC HEAT /////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6
    // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6
    _block_av += _measurement; // Update block accumulators

    return;
}

void System ::averages(const int blk)
{

    ofstream coutf;
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
        _measure.stream_penergy() << setw(12) << blk
                                  << setw(12) << average
                                  << setw(12) << sum_average / double(blk)
                                  << setw(12) << this->error(sum_average, sum_ave2, blk) << "\n";
    }
    // KINETIC ENERGY ////////////////////////////////////////////////////////////
    if (_measure.kenergy)
    {
        average = _average(_measure.idx_kenergy);
        sum_average = _global_av(_measure.idx_kenergy);
        sum_ave2 = _global_av2(_measure.idx_kenergy);
        _measure.stream_kenergy() << setw(12) << blk
                                  << setw(12) << average
                                  << setw(12) << sum_average / double(blk)
                                  << setw(12) << this->error(sum_average, sum_ave2, blk) << "\n";
    }
    // TOTAL ENERGY //////////////////////////////////////////////////////////////
    if (_measure.tenergy)
    {
        average = _average(_measure.idx_tenergy);
        sum_average = _global_av(_measure.idx_tenergy);
        sum_ave2 = _global_av2(_measure.idx_tenergy);
        _measure.stream_tenergy() << setw(12) << blk
                                  << setw(12) << average
                                  << setw(12) << sum_average / double(blk)
                                  << setw(12) << this->error(sum_average, sum_ave2, blk) << "\n";
    }
    // TEMPERATURE ///////////////////////////////////////////////////////////////
    if (_measure.temp)
    {
        average = _average(_measure.idx_temp);
        sum_average = _global_av(_measure.idx_temp);
        sum_ave2 = _global_av2(_measure.idx_temp);
        _measure.stream_temp() << setw(12) << blk
                               << setw(12) << average
                               << setw(12) << sum_average / double(blk)
                               << setw(12) << this->error(sum_average, sum_ave2, blk) << "\n";
    }
    // PRESSURE //////////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 4
    if (_measure.pressure)
    {
        average = _average(_measure.idx_pressure);
        sum_average = _global_av(_measure.idx_pressure);
        sum_ave2 = _global_av2(_measure.idx_pressure);
        _measure.stream_pressure() << setw(12) << blk
                                   << setw(12) << average
                                   << setw(12) << sum_average / double(blk)
                                   << setw(12) << this->error(sum_average, sum_ave2, blk) << "\n";
    }
    // GOFR //////////////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 7
    // MAGNETIZATION /////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6
    // SPECIFIC HEAT /////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6
    // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
    // TO BE FIXED IN EXERCISE 6
    // ACCEPTANCE ////////////////////////////////////////////////////////////////
    double fraction;
    coutf.open("../OUTPUT/acceptance.dat", ios::app);
    if (_nattempts > 0)
        fraction = double(_naccepted) / double(_nattempts);
    else
        fraction = 0.0;
    coutf << setw(12) << blk << setw(12) << fraction << "\n";
    coutf.close();

    return;
}

double System ::error(const double acc, const double acc2, const int blk)
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
