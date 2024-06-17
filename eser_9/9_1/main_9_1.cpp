#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "initializer.h"
#include "mutation.h"
#include "population.h"
#include "random.h"
#include "utilities.h"

constexpr std::size_t SIZE = 32;
constexpr std::size_t POP_SIZE = 100;
constexpr std::size_t N_GEN = 1000;

int main()

{
    Random rng;
    initializer(rng);

    // Positions on a circle
    auto positions = circle_initializer<SIZE>();

    std::for_each(positions.begin(), positions.end(), print_vector<double>());

    // Creating population
    constexpr uint32_t initial_shuffling = 300;
    Population<SIZE> population(POP_SIZE);

    std::cout << "Pre Swapping Initialization\n";
    // Printing Population
    std::for_each(population.begin(), population.end(), [](const Individual<SIZE>& el)
                  { el.print_DNA(); });

    for (auto& element : population)
    {
        for (uint32_t i = 0; i < initial_shuffling; i++)
        {
            element.pair_permutation(rng);
        }
    }

    std::cout << "After swapping initialization\n";
    // Printing Population
    std::for_each(population.begin(), population.end(), print_vector<uint8_t>());
    population.sort_population(positions);
    for (auto& pop : population)
    {
        pop.print_DNA();
        std::cout << "Cost is : " << pop.cost(positions) << '\n';
    }

    std::cout << "Starting Circle Optimization\n\n";

    // Using Genetic Algorithm

    for (size_t i = 0; i < N_GEN; i++)
    {
        population.new_gen(rng);
        population.sort_population(positions);
        std::cout << "Best is : " << population.begin()->cost(positions) << '\n';
    }

    // Printing best circle configuration

    population.begin()->print_DNA();

    std::ofstream foff("best_config_circle.csv");

    foff << "x,y\n";

    for (auto& x : *(population.begin()))
    {
        foff << positions[x].at(0) << ',' << positions[x].at(1) << '\n';
    }

    foff.close();

    // square test

    positions = square_initializer<SIZE>(rng);

    std::cout << "Starting Square Test\n";

    // Shuffling population for square
    for (auto& element : population)
    {
        for (uint32_t i = 0; i < initial_shuffling; i++)
        {
            element.pair_permutation(rng);
        }
    }

    std::for_each(population.begin(), population.end(), print_vector<uint8_t>());
    population.sort_population(positions);
    for (auto& pop : population)
    {
        pop.print_DNA();
        std::cout << "Cost is : " << pop.cost(positions) << '\n';
    }

    std::cout << "Square Circle Optimization\n\n";

    // Using Genetic Algorithm

    for (size_t i = 0; i < N_GEN; i++)
    {
        population.new_gen(rng);
        population.sort_population(positions);
        std::cout << "Best is : " << population.begin()->cost(positions) << '\n';
    }

    // Printing best square configuration

    population.begin()->print_DNA();

    foff.open("best_config_square.csv");

    foff << "x,y\n";

    for (auto& x : *(population.begin()))
    {
        foff << positions[x].at(0) << ',' << positions[x].at(1) << '\n';
    }

    foff.close();

    return 0;
}
