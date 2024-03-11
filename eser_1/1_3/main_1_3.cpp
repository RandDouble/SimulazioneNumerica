#include <iostream>
#include <vector>
#include "utilities.h"
#include "random.h"
#include "buffon.h"

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
    const int n_blocks = 100;
    const int n_throws = 1000; // Number of throws per block
    const int n_seciton = 20;
    Random rng;
    std::vector possible_pi(n_blocks, 0.); // It should inferre that is a vector of doubles
    initializer(rng);
    Buffon buffon(needle_lenght, grid_spacing, n_seciton, &rng);

    std::cout << "starting calculation\n";
    for (auto &pies : possible_pi)
    {
        buffon.launch(n_throws);
        pies = buffon.compute_pi();
        std::cout << pies << "\n";
    }
    return 0;
}