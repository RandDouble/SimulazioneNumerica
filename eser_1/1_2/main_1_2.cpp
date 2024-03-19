#include "random.h"
#include "utilities.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>

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
