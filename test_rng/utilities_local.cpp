#include "utilities_local.h"

void initializer(Random &rnd, const std::size_t rows_to_skip)
{
    int seed[4];
    int p1, p2;
    std::ifstream Primes("Primes");
    if (Primes.is_open())
    {
        for (std::size_t i = 0; i < rows_to_skip && !Primes.eof(); ++i)
        {
            Primes.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        if (!Primes.eof())
        {
            Primes >> p1 >> p2;
        }
        else
        {
            std::cerr << "ERROR: Failed to go to line" << rows_to_skip << " Aborting\n";
            exit(-1);
        }
    }
    else
        std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
    Primes.close();

    std::ifstream input("seed.in");
    std::string property;
    if (input.is_open())
    {
        while (!input.eof())
        {
            input >> property;
            if (property == "RANDOMSEED")
            {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    }
    else
        std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
}

void initializer(Sim::Random &rnd, const std::size_t rows_to_skip)
{
    int seed[4];
    int p1, p2;
    std::ifstream Primes("Primes");
    if (Primes.is_open())
    {
        for (std::size_t i = 0; i < rows_to_skip && !Primes.eof(); ++i)
        {
            Primes.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        if (!Primes.eof())
        {
            Primes >> p1 >> p2;
        }
        else
        {
            std::cerr << "ERROR: Failed to go to line" << rows_to_skip << " Aborting\n";
            exit(-1);
        }
    }
    else
        std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
    Primes.close();

    std::ifstream input("seed.in");
    std::string property;
    if (input.is_open())
    {
        while (!input.eof())
        {
            input >> property;
            if (property == "RANDOMSEED")
            {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    }
    else
        std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;
}

std::ostream &operator<<(std::ostream &os, const values &val)
{
    os << val.value << "," << val.error;
    return os;
}

double calc_mean(std::vector<double> &v_el)
{
    double mean = std::reduce(PARALLEL v_el.begin(), v_el.end(), 0.) / v_el.size();

    return mean;
}

double calc_mean(std::vector<double>::iterator begin, std::vector<double>::iterator end)
{
    assert((end < begin) && "Error in range given to calc mean\n");
    if (begin == end)
    {
        return *begin;
    }
    double mean = std::reduce(PARALLEL begin, end + 1, 0.) / (end + 1 - begin);

    return mean;
}

double calc_std(std::vector<double> &v_el)
{
    double el_squared = std::transform_reduce(PARALLEL v_el.begin(), v_el.end(), v_el.begin(), 0.) / v_el.size(); // Sum of A_i ^ 2, same as \vec{A} \cdot \vec{A}
    double mean = calc_mean(v_el);
    double mean_squared = mean * mean;
    double variance = el_squared - mean_squared;
    variance /= (v_el.size() - 1);

    return std::sqrt(variance);
}

double calc_std(std::vector<double>::iterator begin, std::vector<double>::iterator end)
{
    assert((end < begin) && "Error in range given to calc std\n");

    if (begin == end)
    {
        return 0.;
    }

    double el_squared = std::transform_reduce(PARALLEL begin, end + 1, begin, 0.) / (end - begin + 1); // Sum of A_i ^ 2, same as \vec{A} \cdot \vec{A}
                                                                                                       // Parallelized version of std::inner_product

    double mean = calc_mean(begin, end);
    double mean_squared = mean * mean;
    double variance = el_squared - mean_squared;
    variance /= (end - begin);

#ifndef NDEBUG
    std::cout << std::setw(10) << mean_squared << std::setw(10) << el_squared << std::setw(10) << (end - begin) << "\t" << std::setw(10) << variance << "\n";
#endif

    return std::sqrt(variance);
}
