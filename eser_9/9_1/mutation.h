#include <algorithm>
#include <array>
#include <cstdlib>

#include <armadillo>

#include "initializer.h"
#include "random.h"

template <std::size_t SIZE>
void pair_permutation(std::array<uint8_t, SIZE> &cities, Random &rng)
{
    // Choosing random indexes

    std::size_t idx_1 = static_cast<size_t>(rng.Rannyu(1, SIZE));
    std::size_t idx_2 = PBC(SIZE, static_cast<size_t>(rng.Rannyu(idx_1 + 1, idx_1 + SIZE))); // Choosing between SIZE - 1 elements
                                                                                             // differerent from idx_1

    // Using Swap
    std::swap(cities[idx_1], cities[idx_2]);
}

template <std::size_t SIZE>
void shift_block(std::array<uint8_t, SIZE> &cities, Random &rng, std::size_t block_size, std::size_t shift_size)
{
    assert((block_size < (SIZE - 1)) && "Block Size to shift greater than Array Size\n");
    assert((shift_block > 0) && "We want to shift only on the left.\n");
    std::size_t starting_idx = idx_1 = static_cast<size_t>(rng.Rannyu(1, SIZE));

    // Moving to the end block to move

    std::rotate(cities.begin() + 1, , cities.end());
}

template <std::size_t SIZE>
void permutate_contiguos(std::array<uint8_t, SIZE> &cities, Random &rng)
{
    std::size_t m_contiguos = static_cast<std::size_t>(rng.Rannyu(1, SIZE / 2));
    std::size_t first_pos_idx = static_cast<std::size_t>(rng.Rannyu(1, SIZE / 2));
    auto first_position = cities.begin() + first_pos_idx;
    auto second_position = cities.begin();

    /// @todo: logica per evitare di intersecare medesima regione...
    std::swap_ranges(first_position, first_position + m_contiguos, second_position);
}

template <std::size_t SIZE>
void inversion(std::array<uint8_t, SIZE> &cities, Random &rng)
{
    auto first_position = cities.begin() + 1;
    std::size_t invertion_lenght = static_cast<double>();

    std::reverse(cities.begin() + 1, cities.end());
}
