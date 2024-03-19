#include <iostream>
#include <fstream>
#include <vector>

#include "utilities.h"

void fill_mean_vector(std::vector<double> &vec, Random &rng, const int n_block_size);
void fill_sigma_vector(std::vector<double> &vec, Random &rng, const int n_block_size);
double calc_chi_square(std::vector<int> &vec_occur, double expected_val);

int main()
{
    Random rng;
    initializer(rng);
    constexpr int n_block = 100;
    constexpr int n_throws = 10000;
    constexpr int n_block_size = n_throws / n_block;

    // Esercizio 1.1.1

    std::vector v_mean(n_block, 0.); // 100 block initialized to zero for mean

    fill_mean_vector(v_mean, rng, n_block_size);

    // Esercizio 1.1.2

    // Resetting Seed
    std::vector v_sigma(n_block, 0.); // 100 block initialized to zero for standard deviation
    initializer(rng);
    fill_sigma_vector(v_sigma, rng, n_block_size);

    // Esercizio 1.1.3
    constexpr int n_subdivision = 100;
    constexpr double exp_value = static_cast<double>(n_throws) / n_subdivision;
    initializer(rng);
    std::vector v_chi_square(n_block, 0.); // 100 block initialized to zero for chi square

    for (auto &chi_square : v_chi_square)
    {
        std::vector v_occurences(n_subdivision, 0);
        for (size_t i = 0; i < n_throws; i++)
        {
            double ran_nmbr = rng.Rannyu();
            int index = static_cast<int>(std::floor(std::fmod(ran_nmbr * n_subdivision, n_subdivision)));
            v_occurences[index]++;
        }
        chi_square = calc_chi_square(v_occurences, exp_value);
    }

    // Medie a Blocchi

    auto beg_mean_vec = v_mean.begin();
    std::vector<values> v_out(n_block);
    auto beg_sigma_vec = v_sigma.begin();
    std::vector<values> v_out_sigma(n_block_size);

    for (int i = 0; i < n_block; i++)
    {
        v_out[i] = {calc_mean(beg_mean_vec, beg_mean_vec + i), calc_std(beg_mean_vec, beg_mean_vec + i)};
        v_out_sigma[i] = {calc_mean(beg_sigma_vec, beg_sigma_vec + i), calc_std(beg_sigma_vec, beg_sigma_vec + i)};
    }

#pragma region OutputToFileMean

    std::ofstream f_out("data.csv");
    if (!f_out)
    {
        std::cerr << "Could not open data.csv\n";
        exit(-1);
    }

    for (auto &el : v_out)
    {
        f_out << el << "\n";
    }

    f_out.close();

#pragma endregion

#pragma region OutputToFileSigma

    f_out.open("sigma.csv");

    if (!f_out)
    {
        std::cerr << "Could not open file sigma.csv\n";
        exit(-1);
    }

    for (auto &el : v_out_sigma)
    {
        f_out << el << "\n";
    }
    f_out.close();

#pragma endregion

#pragma region OutputChiSquare

    f_out.open("chi_square.csv");

    if (!f_out)
    {
        std::cerr << "Could not open file chi_square.csv\n";
        exit(-1);
    }

    for (auto &&chi : v_chi_square)
    {
        f_out << chi << "\n";
    }

    f_out.close();

#pragma endregion

    return 0;
}

void fill_mean_vector(std::vector<double> &vec, Random &rng, const int n_block_size)
{
    for (auto &val : vec)
    {
        for (int i = 0; i < n_block_size; i++)
        {
            val += rng.Rannyu();
        }
        val /= n_block_size;
    }
}

void fill_sigma_vector(std::vector<double> &vec, Random &rng, const int n_block_size)
{
    for (auto &val : vec)
    {
        for (int i = 0; i < n_block_size; i++)
        {
            double value = rng.Rannyu();
            val += (value - 0.5) * (value - 0.5);
        }
        val /= n_block_size;
    }
}

double calc_chi_square(std::vector<int> &vec_occur, double expected_val)
{
    double result = 0.;
    for (auto &occ : vec_occur)
    {
        result += ((occ - expected_val) * (occ - expected_val) / expected_val);
    }
    return result;
}
