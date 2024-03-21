#ifndef __UTILITIES__
#define __UTILITIES__

#define NDEBUG

#include <vector>
#include <iostream>
#include <ostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iomanip>
#include "random.h"

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

void initializer(Random &rnd);

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
