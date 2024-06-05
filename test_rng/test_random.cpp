#include <chrono>
#include <iostream>
#include <vector>

#include "old_random.h"
#include "random_gen_sim.h"
#include "utilities_local.h"

using namespace std::chrono;

int main()
{
    Random rng;
    Sim::Random sim_rng;

    initializer(rng);
    initializer(sim_rng);

    constexpr std::size_t size = 100000;
    std::vector result_unmodified(size, 0.);
    std::vector result_sim(size, 0.);
    std::vector equality_check(size, 0);

    auto start = high_resolution_clock::now();

    for (size_t i = 0; i < size; i++)
    {
        result_unmodified[i] = rng.Rannyu();
    }

    auto end_first = high_resolution_clock::now();

    for (size_t i = 0; i < size; i++)
    {
        result_sim[i] = sim_rng.Rannyu();
    }

    auto end_second = high_resolution_clock::now();

    std::cout << "Print Result original Generator\n";
    print(result_unmodified);
    std::cout << "Print Result modified Generator\n";
    print(result_sim);

    for (size_t i = 0; i < size; i++)
    {
        equality_check[i] = (std::abs(result_unmodified[i] - result_sim[i]) < std::numeric_limits<double>::epsilon());
    }

    int num_errors = size - std::accumulate(equality_check.begin(), equality_check.end(), 0);

    auto duration_first = duration<double>(end_first - start).count();
    auto duration_second = duration<double>(end_second - end_first).count();
    auto gain = duration_first / duration_second;
    std::cout << "Number of errors :" << std::setw(8) << num_errors << '\n';
    std::cout << "Elapsed time extraction unmodified :" << std::setw(15) << duration_first << '\n'
              << "Elapsed time extraction new :" << std::setw(15) << duration_second << '\n'
              << "Time Gain : " << std::setw(10) << gain << '\n';
    // print(equality_check);

    sim_rng.SaveSeed();

    return 0;
}
