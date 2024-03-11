#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include "random.h"
#include "monte_carlo.h"
#include "function.h"

int main()
{

    std::ifstream primes("../lib/Primes");
    int seed[4];
    int p1, p2;

    if (!primes)
    {
        std::cerr << "Prime file not found, exiting\n";
        exit(-1);
    }

    primes >> p1 >> p2;
    primes.close();

    std::ifstream seed_list("../lib/seed.in");
    std::string property;

    if (!seed_list)
    {
        while (!seed_list.eof())
        {
            seed_list >> property;
            if (property == "RANDOMSEED")
            {
                seed_list >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            }
        }
        seed_list.close();
    }
    Random rng;
    rng.SetRandom(seed, p1, p2);
    rng.SaveSeed();
    double a = 0., b = 1.;

    MonteCarlo integrator(&rng);
    for (int i = 1; i < 20; i++)
    {
        integrator.set_launch_size( 200 * (1 << i));
        std::cout << "Integrale Passo : " << std::setw(8) << integrator.get_launch() << "\tPasso : " << std::setw(8) << integrator.calculate_uniform(a, b, integrand).format_string() << "\n";
    }

    return 0;
}