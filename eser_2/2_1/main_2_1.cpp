#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include "random.h"
#include "monte_carlo.h"
#include "function.h" // To instantiate function to be integrated and distribution

int main()
{

    Random rng;
    initializer(rng);
    double a = 0., b = 1.;

    MonteCarlo integrator(&rng);
    integrator.set_n_block(1000);

    // 2_1_1, calculate using uniform extraction
    std::vector v_uniform = integrator.calculate_uniform_blocks(a, b, 100, integrand);
    
    // File output
    std::ofstream f_out("uniform.csv");
    for (auto &&val : v_uniform)
    {
        f_out << val << "\n";
    }
    f_out.close();

    // 2_1_2, calculate using importance sampling, on top of that I used antithetic integrand, combined with accept reject for extracting numbers
    std::vector v_precise = integrator.calculate_dist_rng_block(a, b, 100, integrand_precise, [&]()
                                                                { return rng.AcceptReject(a, b, integrand_anti_expansion_max(), integrand_anti_expansion_norm); });
    // Lambda function are useful in this context.
    
    // File output
    f_out.open("precise.csv");
    for (auto &&val : v_precise)
    {
        f_out << val << "\n";
    }
    f_out.close();

    return 0;
}