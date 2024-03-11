#include <iostream>
#include <fstream>
#include <vector>

#include "utilities.h"

void fill_mean_vector(std::vector<double> &vec, Random &rng, const int n_block_size);
void fill_sigma_vector(std::vector<double> &vec, Random &rng, const int n_block_size);

int main()
{
    Random rng;
    initializer(rng);
    const int n_block = 100;
    const int n_throws = 10000;
    const int n_block_size = n_throws / n_block;
    std::ofstream f_out("data.csv");
    if (!f_out)
    {
        std::cerr << "Could not open data.csv\n";
        exit(-1);
    }

    std::vector v_mean(n_block, 0.);  // 100 blocchi inizializzati a zero
    std::vector v_sigma(n_block, 0.); // 100 blocchi inizializzati a zero

    fill_mean_vector(v_mean, rng, n_block_size);
    
    // Resetting Seed
    initializer(rng);
    fill_sigma_vector(v_sigma, rng, n_block_size);

    auto beg_mean_vec = v_mean.begin();
    std::vector<values> v_out(n_block);
    auto beg_sigma_vec = v_sigma.begin();
    std::vector<values> v_out_sigma(n_block_size);

    for (int i = 0; i < n_block; i++)
    {
        v_out[i] = {calc_mean(beg_mean_vec, beg_mean_vec + i), calc_std(beg_mean_vec, beg_mean_vec + i)};
        v_out_sigma[i] = {calc_mean(beg_sigma_vec, beg_sigma_vec + i), calc_std(beg_sigma_vec, beg_sigma_vec + i)};
    }

    for (auto &el : v_out)
    {
        f_out << el << "\n";
    }
    f_out.close();
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
