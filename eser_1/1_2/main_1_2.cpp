/// @author Stefano Pilosio
/// @date 29th March 2024
/// @brief Solution fo exercise 1.2
///
/// 1. Add two probability distributions by using the **method of the inversion
/// of the cumulative distribution** to sample from a **generic** exponential distribution,
/// $p(x) = \lambda \exp(-\lambda x)$, $x\in [0;+\infty]$
/// (see <a href="https://en.wikipedia.org/wiki/Exponential_distribution">
/// this Wikipedia link</a>),
/// and a **generic** Cauchy-Lorentz distribution
/// $p(x)=\frac{1}{\pi}\frac{\Gamma}{(x-\mu)^2+\Gamma^2}$,
/// $x\in [-\infty;+\infty]$
/// (see <a href="https://en.wikipedia.org/wiki/Cauchy_distribution">
/// this Wikipedia link</a>).
/// @paragraph 2. Make 3 pictures with the histograms obtained filling them with $10^4$
/// realizations of $S_N = \frac{1}{N}\sum_{i=1}^N x_i$ (for $N=1, 2, 10, 100$), being $x_i$
/// a random variable sampled throwing a *standard* dice (fig.1), an *exponential* dice
/// (fig.2, use $\lambda=1$) and a *Lorentzian* dice (fig.3, use $\mu=0$ and $\Gamma=1$).
/// **************************************************************************************
/// @todo tha you can try to fit the case $N=100$ with a Gaussian for standard and exponential dices, whereas you should
/// use a Cauchy-Lorentz distribution for the last case.

#include <array>
#include <iomanip>
#include <iostream>
#include <vector>

#include "random.h"
#include "utilities.h"

/// @brief Function to print to file
/// @param file_out Ofstream not yet opened
/// @param filename Output file name
/// @param vec Vector to be printed
void print_file(std::ofstream &file_out, std::string filename, std::vector<double> vec);

int main()
{
    Random rng;
    initializer(rng);
    const std::array iterations{1, 2, 10, 100};
    const int realization = static_cast<int>(1e4);
    const std::string prefix_unif{"uniform_"};
    const std::string prefix_exp{"exponential_"};
    const std::string prefix_cauch{"cauchy_"};
    const double lambda = 1.;
    const double mu = 0.;
    const double gamma = 1.;

    std::ofstream file_off;
    std::vector container(realization, 0.);

    // Uniform

    for (const auto &iter : iterations)
    {
        for (auto &el : container)
        {
            // S_N = sum(x_i) / N, N in {1, 2 , 10 , 100}
            for (int i = 0; i < iter; i++)
            {
                el += std::floor(std::fmod(rng.Rannyu() * std::pow(6, 3), 6)) + 1.;
            }
            el /= iter;
        }
        print_file(file_off, prefix_unif + std::to_string(iter) + ".dat", container);
        container = std::vector(realization, 0.);
    }

    // Exponential

    for (const auto &iter : iterations)
    {
        for (auto &el : container)
        {

            for (int i = 0; i < iter; i++)
            {
                el += rng.Exponential(lambda);
            }
            el /= iter;
        }
        print_file(file_off, prefix_exp + std::to_string(iter) + ".dat", container);
        container = std::vector(realization, 0.);
    }

    // Lorentzian

    for (const auto &iter : iterations)
    {
        for (auto &el : container)
        {

            for (int i = 0; i < iter; i++)
            {
                el += rng.Lorenztian(mu, gamma);
            }
            el /= iter;
        }
        print_file(file_off, prefix_cauch + std::to_string(iter) + ".dat", container);
        container = std::vector(realization, 0.);
    }

    return 0;
}

void print_file(std::ofstream &file_out, std::string filename, std::vector<double> vec)
{
    file_out.open(filename);
    if (!file_out)
    {
        std::cerr << "Could Not open file: " << filename << ", Closing\n";
        exit(-1);
    }

    for (const auto &el : vec)
    {
        file_out << el << "\n";
    }
    file_out.close();
}
