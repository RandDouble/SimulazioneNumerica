#include <iostream>
#include <fstream>
#include <vector>

#include "random.h"
#include "utilities.h"
#include "random_walk.h"

/*
    1. on a cubic lattice with lattice constant $a=1$;
    at each discrete time the walker makes a forward or backward step of length equal to $a$
    in one of the 3 principal directions of the lattice: $x$, $y$ or $z$

    2. in the continuum; at each discrete time the walker makes a step of length equal to $a(=1)$
    along a **random direction** obtained by sampling **uniformly**
    the solid angle: $\theta \in [0,\pi]$ and $\phi \in [0,2\pi]$
 */

int main()
{
    Random rng;
    initializer(rng);
    Random_Walk walker(1., {0., 0., 0.}, &rng);

    // constexpr unsigned int n_repetition = 10000;
    constexpr unsigned int n_time_steps = 100;


    // Point 1
        

    std::ofstream f_off("random_walk_lattice.csv");

    f_off << "x,y,z\n";
    for (size_t i = 0; i < n_time_steps; i++)
    {
        walker.discrete_single_walk();
        f_off << walker.get_actual_pos() << "\n";
    }

    f_off.close();
    
    f_off.open("random_walk_continuos.csv");
    walker.reset();
    for (size_t i = 0; i < n_time_steps; i++)
    {
        walker.continuos_single_walk();
        f_off << walker.get_actual_pos() << "\n";
    }
    

    f_off.close();
    return 0;
}
