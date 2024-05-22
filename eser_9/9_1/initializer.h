#include <array>
#include <cassert>
#include <cstdlib>
#include <vector>

#include <armadillo>

template <typename T = uint8_t, std::size_t SIZE>
std::array<T, SIZE> create_array()
{
    std::array<T, SIZE> result;
    std::iota(result.begin(), result.end(), 0);
    return result;
}

template <std::size_t SIZE>
std::array<arma::vec2, SIZE> circle_initializer(const int n_points)
{
    std::array<arma::vec2, SIZE> result;

    for (size_t i = 0; i < SIZE; i++)
    {
        double angle = 2. * M_PI * static_cast<double>(i) / static_cast<double>(SIZE);
        result[i] = {arma::cos(angle), arma::sin(angle)}; // @todo: sistemare, i punti sono da fare a caso
    }
    return result;
}

template <std::size_t SIZE>
std::array<arma::vec2, SIZE> square_initializer(const int n_points)
{
    return std::array<arma::vec2, SIZE>;
}

std::size_t PBC(std::size_t SIZE, std::size_t idx)
{
    assert((idx == 0) && "You stupid man, must not touch the first city!\nHow you dare!!!\n");

    while (idx >= SIZE) // Caso migliore non entra neanche in questo ciclo...
    {
        idx -= SIZE + 1;
    }
    return idx;
}