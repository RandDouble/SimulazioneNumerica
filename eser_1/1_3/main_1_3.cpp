#include <iostream>
#include <vector>
#include "utilities.h"
#include "random.h"
#include "buffon.h"

#define NDEBUG
/* Simulate Buffon Esperiment:
 * A needle of Lenght L is thrown at random onto a horizontal plane ruled with straight lines a distance d,
 * (must be d > L, but do not use d >> L otherwise P << 1) apart. the probability P that the needle will intersect one of these lines is:
 * P = (2L )/ (pi * d). This could be used to evaluate PI from throws od the needle: if the needle id throwm down N_thr times
 * ans is observed to land on a line N_hit of those times, we can make an estimation fo Pi from:
 *          PI = 2L / (P * d) = lim N_thr -> oo (2L * N_thr) / (N_hit * d )
 */

int main()
{
    const double needle_lenght = 0.1, grid_spacing = 0.2;
    constexpr int n_blocks = 100;
    constexpr int n_throws = 10000; // Number of throws per block
    constexpr int n_section = 20;
    Random rng;
    initializer(rng);
    
    std::vector possible_pi(n_blocks, 0.); // It should inferre that is a vector of doubles
    std::vector pi_blocks(n_blocks, values{0.,0.});

    Buffon buffon(needle_lenght, grid_spacing, n_section, &rng);

    std::cout << "starting calculation\n";
    for (auto &pies : possible_pi)
    {
        buffon.launch(n_throws);
        pies = buffon.compute_pi();
        #ifndef NDEBUG
        std::cout << pies << "\n";
        #endif //NDEBUG
    }
    
    auto vec_begin = possible_pi.begin();
    
    #ifndef NDEBUG
    std::cout << "Starting Mean block\n";
    #endif // NDEBUG


    for (int i = 0; i < n_blocks; i++)
    {
        pi_blocks[i].value = calc_mean(vec_begin, vec_begin + i);
        pi_blocks[i].error = calc_std(vec_begin, vec_begin + i);
        #ifndef NDEBUG 
        std::cout <<std::right<<std::setw(2)<< i <<".\t"<< pi_blocks[i].format_string() << "\n";
        #endif // NDEBUG
    }

    std::ofstream fout("buffon.csv");
    if (!fout)
    {
        std::cerr << "Could not open buffon.csv\nExiting\n";
        exit(1);
    }

    for (auto &el : pi_blocks)
    {
        fout << el << "\n";
    }
    fout.close();
    return 0;
}