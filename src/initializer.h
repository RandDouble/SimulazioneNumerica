#include <array>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include <armadillo>

#include "random.h"

#ifndef __INITIALIZER__
#define __INITIALIZER__

template <typename T = uint8_t, std::size_t SIZE>
std::array<T, SIZE> create_array()
{
    std::array<T, SIZE> result;
    std::iota(result.begin(), result.end(), 0);
    return result;
}

template <std::size_t SIZE>
std::array<arma::vec2, SIZE> circle_initializer()
{
    std::array<arma::vec2, SIZE> result;

    for (size_t i = 0; i < SIZE; i++)
    {
        double angle = 2. * M_PI * static_cast<double>(i) / static_cast<double>(SIZE);
        result[i] = {std::cos(angle),
                     std::sin(angle)}; // @todo: sistemare, i punti sono da fare a caso
    }
    return result;
}

template <std::size_t SIZE>
std::array<arma::vec2, SIZE> square_initializer(Random& rng)
{
    std::array<arma::vec2, SIZE> elements;
    for (auto& el : elements)
    {
        el = {rng.Rannyu(-1., 1.), rng.Rannyu(-1., 1.)};
    }

    return elements;
}

std::vector<arma::vec2> load_province_position(const std::string& filename)
{
    std::vector<arma::vec2> res(0);
    std::ifstream fin(filename);

    if (!fin)
    {
        std::cerr << "Could Not Open File " << filename << '\n' << "Aborting\n";
        exit(-1);
    }

    for (arma::vec2 appo; fin >> appo(0) >> appo(1);)
    {
        res.push_back(appo);
    }
    return res;
}

std::vector<std::string> load_province_name(const std::string& filename)
{
    std::vector<std::string> res(0);
    std::ifstream fin(filename);

    if (!fin)
    {
        std::cerr << "Could Not Open File " << filename << '\n' << "Aborting\n";
        exit(-1);
    }

    for (std::string appo; std::getline(fin, appo, '\n');)
    {
        res.push_back(appo);
    }
    return res;
}

arma::mat create_matrix(std::vector<arma::vec2>& vec)
{
    arma::mat result_mat;
    result_mat.zeros(vec.size(), vec.size());

    // std::cout << "Matrix size during create_matric is : " << arma::size(result_mat) <<
    // '\n';

    for (size_t i = 0; i < vec.size(); i++)
    {
        for (size_t j = i; j < vec.size(); j++)
        {
            // std::cout << "current i :" << i << "current j : " << j << '\n';
            arma::vec2 distance = vec[i] - vec[j];
            result_mat(i, j) = arma::norm(distance);
            result_mat(j, i) = result_mat(i, j);
        }
    }
    return result_mat;
}

/// @brief Periodic Boundary Conditions
/// @param SIZE Grandezza Effettiva dell'array su cui si può giocare
/// @param idx  Indice dell'oggetto, verrà portato tra [1, SIZE]
/// @return Restituisce l'indice riportato tra [1, SIZE]
std::size_t PBC(const std::size_t SIZE, std::size_t idx)
{
    // Io voglio che idx stia tra 1 e SIZE-1
    while (idx >= SIZE) // Caso migliore non entra neanche in questo ciclo...
    {
        idx -= SIZE;
    }
    while (idx < 1)
    {
        idx += SIZE;
    }

    assert((idx != 0ull)
           && "You stupid man, must not touch the first city!\nHow you dare!!!\n");
    return idx;
}

template <typename T>
struct print_vector
{
    template <typename iterable, int WIDTH = 8>
    void operator()(const iterable& iter)
    {
        for (const T& el : iter)
        {
            std::cout << std::setw(WIDTH) << el << " ";
        }
        std::cout << "\n";
    }

    template <typename iterable>
    void operator()(const iterable& iter, const int width)
    {
        for (const T& el : iter)
        {
            std::cout << std::setw(width) << el << " ";
        }
        std::cout << "\n";
    }
};

template <>
struct print_vector<uint8_t>
{
    template <typename iterable, int WIDTH = 3>
    void operator()(const iterable& iter)
    {
        for (const uint8_t& el : iter)
        {
            std::cout << std::setw(WIDTH) << static_cast<int>(el) << " ";
        }
        std::cout << "\n";
    }

    template <typename iterable>
    void operator()(const iterable& iter, const int width)
    {
        for (const uint8_t& el : iter)
        {
            std::cout << std::setw(width) << static_cast<int>(el) << " ";
        }
        std::cout << "\n";
    }
};

#endif // __INITIALIZER__
