#include "random.h"
#include <algorithm>
#include <execution>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

#ifndef __UTILITIES__
#define __UTILITIES__


#ifdef USE_PARALLEL
#define PARALLEL std::execution::par,
#else
#define PARALLEL
#endif

struct values
{
    double value{0.}, error{0.};

    inline std::string format_string() const
    {
        return std::to_string(value) + " +- " + std::to_string(error);
    }

    friend std::ostream &operator<<(std::ostream &os, const values &val);
};

std::ostream &operator<<(std::ostream &os, const values &val);

void initializer(Random &rnd, const std::size_t rows_to_skip = 0);

double calc_mean(std::vector<double> &v_el);
double calc_mean(std::vector<double>::iterator begin, std::vector<double>::iterator end);

double calc_std(std::vector<double> &v_el);
double calc_std(std::vector<double>::iterator begin, std::vector<double>::iterator end);

template <typename T>
void print_file(std::ostream &os, const std::vector<T> &vec)
{
    for (auto &&i : vec)
    {
        os << i << "\n";
    }
}

template <typename T>
void print(const std::vector<T> &vec)
{
    for (auto &&i : vec)
    {
        std::cout << i << "\n";
    }
}

#endif // __UTILITIES__
