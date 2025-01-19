#define ARMA_DONT_PRINT_FAST_MATH_WARNING
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

constexpr std::size_t SIZE = 34; // Points to visit
constexpr std::size_t POP_SIZE = 10000;
constexpr std::size_t HALF_POP_SIZE = POP_SIZE / 2;
constexpr std::size_t N_GEN_CIRCLE = 100;
constexpr std::size_t N_GEN_SQUARE = 100;
constexpr double SELECTOR_COEFF = -.4;

constexpr double CROSSOVER_PROB = 0.65;
constexpr double SWAP_PROB = 0.05;
constexpr double SHIFT_PROB = 0.05;
constexpr double PERMUTATE_PROB = 0.05;
constexpr double INVERSE_PROB = 0.05;

int main()
{
#pragma region INITIALIZATION
    Random rng;
    initializer(rng);

    // Positions on a circle
    auto circle_positions = circle_initializer<SIZE>();
    auto square_positions = square_initializer<SIZE>(rng);

    std::for_each(circle_positions.begin(), circle_positions.end(), print_vector<double>());

    // Creating population
    constexpr uint32_t initial_shuffling = 300;
    Population<SIZE> population(POP_SIZE);

    population.crossover_prob(CROSSOVER_PROB);
    population.inverse_prob(INVERSE_PROB);
    population.permutate_prob(PERMUTATE_PROB);
    population.shift_prop(SHIFT_PROB);
    population.swap_prob(SWAP_PROB);
    population.selection_coeff(SELECTOR_COEFF);

    // std::cout << "Pre Swapping Initialization\n";
    // Printing Population
    // std::for_each(population.begin(), population.end(), [](const Individual<SIZE>& el)
    //   { el.print_DNA(); });

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
    population.sort_population(circle_positions);
    for (auto &pop : population)
    {
        // pop.print_DNA();
        std::cout << "Cost is : " << pop.cost(circle_positions) << '\n';
    }

#pragma endregion INITIALIZATION

#pragma region CIRCLE_OPTIMIZATION

    std::cout << "Starting Circle Optimization\n\n";

    // Using Genetic Algorithm
    std::ofstream foff("cost_circle.csv");
    foff << "Gen,Best,Average_on_half\n";

    for (size_t i = 0; i < N_GEN_CIRCLE; i++)
    {
        population.new_gen(rng);
        population.sort_population(circle_positions);
        std::cout << "Best is : " << population.begin()->cost(circle_positions) << '\n';
        double average_on_half = std::accumulate(population.begin(), population.begin() + HALF_POP_SIZE, 0.,
                                             [&circle_positions](double acc, Individual<SIZE> &el) -> double {
                                                 return acc + el.cost(circle_positions);
                                             }) /
                                 HALF_POP_SIZE;
        foff << i << ',' << population.begin()->cost(circle_positions) << ',' << average_on_half << '\n';
    }
    foff.close();

    // Printing best circle configuration

    population.begin()->print_DNA();

    foff.open("best_config_circle.csv");

    foff << "x,y\n";

    for (auto &x : *(population.begin()))
    {
        foff << circle_positions[x].at(0) << ',' << circle_positions[x].at(1) << '\n';
    }

    foff.close();

#pragma endregion CIRCLE_OPTIMIZATION

#pragma region SQUARE_OPTIMIZATION

    // square test
    std::cout << "Starting Square Test\n";

    // Shuffling population for square
    for (auto &element : population)
    {
        for (uint32_t i = 0; i < initial_shuffling; i++)
        {
            element.pair_permutation(rng);
        }
    }

    // std::for_each(population.begin(), population.end(), print_vector<uint8_t>());
    population.sort_population(square_positions);
    for (auto &pop : population)
    {
        // pop.print_DNA();
        std::cout << "Cost is : " << pop.cost(square_positions) << '\n';
    }

    std::cout << "Square Optimization\n\n";

    // Using Genetic Algorithm

    foff.open("cost_square.csv");
    foff << "Gen,Best,Average_on_half\n";
    for (size_t i = 0; i < N_GEN_SQUARE; i++)
    {
        population.new_gen(rng);
        population.sort_population(square_positions);
        std::cout << "Best is : " << population.begin()->cost(square_positions) << '\n';
        
        double average_on_half =
            std::accumulate(population.begin(), population.begin() + HALF_POP_SIZE, 0.,
                        [&square_positions](const double &acc, const Individual<SIZE> &el) -> double {
                            return acc + el.cost(square_positions);
                        }) /
            HALF_POP_SIZE;
        foff << i << ',' << population.begin()->cost(square_positions) << ',' << average_on_half << '\n';
    }
    foff.close();

    // Printing best square configuration
    population.sort_population(square_positions);

    population.begin()->print_DNA();

    foff.open("best_config_square.csv");

    foff << "x,y\n";

    for (auto &x : *(population.begin()))
    {
        std::cout << static_cast<uint16_t>(x) << " ";
        foff << square_positions[x].at(0) << ',' << square_positions[x].at(1) << '\n';
    }

    foff.close();

#pragma endregion SQUARE_OPTIMIZATION

    return 0;
}
