#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "initializer.h"
#include "mutation.h"
#include "random.h"
#include "utilities.h"

constexpr std::size_t SIZE = 32;
constexpr std::size_t POP_SIZE = 100;

int main()

{
    Random rng;
    initializer(rng);

    // Positions on a circle
    auto positions = circle_initializer<SIZE>();

    std::for_each(positions.begin(), positions.end(), print_vector<double>());

    // Creating population
    constexpr uint32_t initial_shuffling = 300;
    std::vector <Individual<SIZE>> population(POP_SIZE);

    std::cout << "Pre Swapping Initialization\n";
    // Printing Population
    std::for_each(population.begin(), population.end(), [](const Individual<SIZE>& el)
                  { el.print_DNA(); });

    for (auto &element : population)
    {
        for (uint32_t i = 0; i < initial_shuffling; i++)
        {
            element.pair_permutation(rng);
        }
    }

    std::cout << "After swapping initialization\n";
    // Printing Population
    std::for_each(population.begin(), population.end(), print_vector<uint8_t>());



    return 0;
}
