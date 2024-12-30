#include <metropolis.h>
#include <utilities.h>

#include <iostream>
#include <vector>

#include "wave_function.h"

int main()
{
    constexpr PsiParam psi_0{.mu = 1., .sigma = 1.};

    WaveFunction wave(psi_0);

    MetropolisMemory<double> metro(0.);
    metro.get_rng_ptr()->Initializer();
    metro.set_n_step(10);

    std::vector points(100000, 0.);

    auto sampler = [&]() { return metro.get_rng_ptr()->Rannyu(-1., 1.); };

    auto PDF = [&](const double x) { return wave.PDF(x); };

    for (auto &point : points)
    {
        point = metro.generate(PDF, sampler);
    }

    std::ofstream fout("test.dat");
    fout << "x\n";
    for (const auto &point : points)
    {
        fout << point << '\n';
    }
    fout.close();
    return 0;
}
