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

    constexpr unsigned int n_time_steps = 100;
    constexpr unsigned int n_runs = 10000;


    // Point 1
        
    // Vogliamo visualizzazioni fighe per vedere se tutto funziona
    std::ofstream
     f_off("random_walk_lattice.csv");

    f_off << "x,y,z\n";
    for (size_t i = 0; i < n_time_steps; i++)
    {
        walker.discrete_single_walk();
        f_off << walker.get_actual_pos() << "\n";
    }

    f_off.close();

    // Vogliamo visualizzazioni fighe per vedere se tutto funziona
    f_off.open("random_walk_continuos.csv");
    walker.reset();
    for (size_t i = 0; i < n_time_steps; i++)
    {
        walker.continuos_single_walk();
        f_off << walker.get_actual_pos() << "\n";
    }
    f_off.close();

    // Adesso facciamo il vero punto.
    std::vector res(n_runs, 0.);
    for (size_t i = 0; i < n_runs; i++)
    {
        walker.reset();
        walker.discrete_walk(n_time_steps);
        res[i] = walker.get_distance_squared_from_origin();
    }

    f_off.open("result_lattice.csv");
    print_file(f_off, res);
    f_off.close();

    for (size_t i = 0; i < n_runs; i++)
    {
        walker.reset();
        walker.continuos_walk(n_time_steps);
        res[i] = walker.get_distance_squared_from_origin();
    }

    f_off.open("result_continuos.csv");
    print_file(f_off, res);
    f_off.close();


    return 0;
}
