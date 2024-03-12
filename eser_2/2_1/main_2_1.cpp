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

    
    Random rng;
    initializer(rng);
    double a = 0., b = 1.;

    MonteCarlo integrator(&rng);
    for (int i = 1; i < 20; i++)
    {
        integrator.set_launch_size( 200 * (1 << i));
        std::cout << "Integrale Passo : " << std::setw(8) << integrator.get_launch() << "\tPasso : " << std::setw(8) << integrator.calculate_uniform(a, b, integrand).format_string() << "\n";
    }

    return 0;
}