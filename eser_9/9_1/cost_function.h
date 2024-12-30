#include <array>
#include <cstdlib>

#include <armadillo>

#ifndef __COST_FUNCTION__
#define __COST_FUNCTION__

template <typename T = uint8_t, std::size_t SIZE, template <typename> typename iterable>
std::array<double, SIZE> array_cost(const iterable<std::array<T, SIZE>> &elements,
                                    const std::array<arma::vec2, SIZE> &positions)
{
    std::array<double, SIZE> result;
    for (int i = 0; i < SIZE; i++)
    {
        result[i] = cost(elements[i], positions);
    }

    return result;
}

#endif // __COST_FUNCTION__
