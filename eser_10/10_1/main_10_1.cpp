#include "mpi.h"

#include <cassert>
#include <fstream>
#include <iostream>

#ifndef NDEBUG
#define ARMA_PRINT_EXCEPTIONS
#define ARMA_PRINT_ERRORS
#endif

#include <json.hpp>

#include "initializer.h"
#include "population.h"

#include "random.h"
#include "utilities.h"

#define MASTER 0
#define N_PROV 110

// I Hate a little bit to use Macro Constants, they are error-prone, generally unsafe
// and have no type or compile checking, but in this case they can be useful..

struct Evolution_Parameters
{
    std::size_t N_GEN;
    std::size_t POP_SIZE;
    std::size_t N_CONTACTS;
    std::size_t N_TRANSFER_ELEMENTS;
    unsigned int initial_shuffling;
    double selection_coeff;
    void print(const int rank) const;
};

struct Mutation_Probabilities
{
    double crossover_prob;
    double swap_prob;
    double permutate_prob;
    double inverse_prob;
    double shift_prob;

    void print(const int rank) const;
};

arma::mat load_distance_matrix(const int rank,
                               std::vector<arma::vec2>& prov_pos,
                               std::vector<std::string>& prov_name);

void communicate_distance_matrix(const int rank, arma::mat& dist_mat);

int main(int argc, char* argv[])
{
    int size, rank, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(hostname, &len);

    arma::mat distance_matrix;
    std::vector<arma::vec2> prov_pos;
    std::vector<std::string> prov_name;

    // Loading configuration
    std::ifstream fin("input.json");
    nlohmann::json data = nlohmann::json::parse(fin);
    fin.close();

    /**** Good Initial values
    // const std::size_t N_GEN = 100ull;
    // const std::size_t POP_SIZE = 500ull;
    // const std::size_t N_CONTACTS = 15ull;
    // const std::size_t N_TRANSFER_ELEMENTS = 5ull;
    // const unsigned int initial_shuffling = 300;

    // const double selection_coeff = -0.75;

    // const double crossover_prob = 0.6;
    // const double swap_prob = 0.1;
    // const double permutate_prob = 0.1;
    // const double inverse_prob = 0.1;
    // const double shift_prop = 0.1;
    */

    const Evolution_Parameters evol_params({data["number_of_generations"],
                                            data["population_size"],
                                            data["number_of_contacts"],
                                            data["n_transfer_elements"],
                                            data["initial_shuffling"],
                                            data["selection_coeff"]});

    const Mutation_Probabilities probs({data["crossover_probability"],
                                        data["swap_probability"],
                                        data["permutate_probability"],
                                        data["inverse_probability"],
                                        data["shift_probability"]});

    // Using printf to avoid buffering problems, in fact printf writes atomically.
    std::printf("Hello World!!\nI am rank %d of %d\non CPU : %s\n", rank, size, hostname);
    evol_params.print(rank);
    probs.print(rank);

    // Distance matrix calculation
    distance_matrix = load_distance_matrix(rank, prov_pos, prov_name);
    communicate_distance_matrix(rank, distance_matrix);

    // Generating Random Number Generator
    Random rng;
    initializer(rng,
                rank); // selecting different rows based on rank to have indipendent RNG

    // Initializing Population... Due to implementation this is not a generic algorithm,
    // we need to know size... But we know Italian provinces are 110
    Population<N_PROV> population(evol_params.POP_SIZE);

    // Setting Algorithm Parameters
    population.selection_coeff(evol_params.selection_coeff); // Selection Coefficient
    population.crossover_prob(probs.crossover_prob);         // Crossover Probability

    population.swap_prob(probs.swap_prob);
    population.permutate_prob(probs.permutate_prob);
    population.inverse_prob(probs.inverse_prob);
    population.shift_prop(probs.shift_prob);

    // Checking for health
    for (auto& pop : population)
    {
        assert(pop.check_health() && "Found ill vector\n");
    }

    // Initial shuffling
    for (auto& element : population)
    {
        for (unsigned int i = 0; i < evol_params.initial_shuffling; i++)
        {
            element.pair_permutation(rng);
        }
    }

    std::printf(
        "Rank %d => Population initialized, RNG initialized\nStarting Evolution\n", rank);

    std::ofstream gen_report("gen_report_thread_" + std::to_string(rank) + ".csv");
    gen_report << "Generation,Cost,Sequence\n";

    population.sort_population(distance_matrix);
    // Evolution
    for (size_t n = 0; n < evol_params.N_CONTACTS; n++)
    {
        for (size_t i = 0; i < evol_params.N_GEN; i++)
        {

            // Reproducing and generating new generation
            population.new_gen(rng);

            // Reporting Best for this generation
            population.sort_population(distance_matrix);
            gen_report << i + evol_params.N_GEN * n << ','
                       << population.begin()->cost(distance_matrix) << ",\"";
            for (auto& gen : *population.begin())
            {
                gen_report << static_cast<uint16_t>(gen) << ' ';
            }
            gen_report << "\"\n";
        }

        // Transfering Best of each continent to all other continents
        auto initial_pos = population.begin();
        std::vector<uint8_t> best_el(initial_pos->size()
                                     * evol_params.N_TRANSFER_ELEMENTS);

        for (size_t i = 0; i < evol_params.N_TRANSFER_ELEMENTS; i++)
        {
            for (size_t j = 0; j < initial_pos->size(); j++)
            {
                best_el[i * N_PROV + j] = (*(initial_pos + i))[j];
            }
        }

        std::vector<uint8_t> receiver(initial_pos->size()
                                      * evol_params.N_TRANSFER_ELEMENTS * size);

        // Contact between continents
        MPI_Allgather(best_el.data(),
                      best_el.size(),
                      MPI_UINT8_T,
                      receiver.data(),
                      best_el.size(),
                      MPI_UINT8_T,
                      MPI_COMM_WORLD);

        auto copying_pos
            = population.end() - (evol_params.N_TRANSFER_ELEMENTS * (size - 1));
        size_t counter = 0;

        for (size_t provenience = 0; provenience < static_cast<size_t>(size);
             provenience++)
        {
            if (provenience == static_cast<size_t>(rank)) // Do not copy your datas
            {
                continue;
            }

            // copy elements from Individual i-th come from provenience-th rank
            for (size_t i = 0; i < evol_params.N_TRANSFER_ELEMENTS; i++)
            {
                for (size_t j = 0; j < copying_pos->size(); j++)
                {
                    (*(copying_pos + counter * evol_params.N_TRANSFER_ELEMENTS + i))[j]
                        = receiver[counter * N_PROV * evol_params.N_TRANSFER_ELEMENTS
                                   + i * N_PROV + j];
                }
            }
            counter++;
        }

        population.sort_population(distance_matrix);
        // std::cout << '\n'
        //           << "Transferring Complete\n";
    }

    // Last independent Evolution before final results
    for (size_t i = 0; i < evol_params.N_GEN; i++)
    {

        // Reproducing and generating new generation
        population.new_gen(rng);

        // Reporting Best for this generation
        population.sort_population(distance_matrix);
        gen_report << i + evol_params.N_GEN * evol_params.N_CONTACTS << ','
                   << population.begin()->cost(distance_matrix) << ",\"";
        for (auto& gen : *population.begin())
        {
            gen_report << static_cast<uint16_t>(gen) << ' ';
        }
        gen_report << "\"\n";
    }
    gen_report.close();

    // population.sort_population(distance_matrix);

    std::printf("Rank %d => Evolution Ended\nBest is : %f\n",
                rank,
                population.begin()->cost(distance_matrix));

    // Getting Best of each continent comparing and then outputting best of all

    Population<N_PROV> best_pop((rank == MASTER) ? size : 0);
    std::vector<uint8_t> recv_best_vec((rank == MASTER) ? size * N_PROV : 0, 0);

    if (rank == MASTER)
    {
        std::printf("Rank %d => Asking for other best\n", rank);
        // ss << " : Asking for other best\n";
    }
    else
    {
        std::printf("Rank %d => Sending my best\n", rank);
        // ss << " : Sending my best\n";
    }

    MPI_Gather(population.begin()->data(),
               N_PROV,
               MPI_UINT8_T,
               recv_best_vec.data(),
               population.begin()->size(),
               MPI_UINT8_T,
               MASTER,
               MPI_COMM_WORLD);

    if (rank == MASTER)
    {
        std::printf("Rank %d => Communication Succeed\nCalculating best\n", rank);
        // std::cout << "Rank :" << std::setw(3) << rank << " : Communication Succeed\n"
        //           << "Calculating best\n";

        size_t counter = 0;
        for (auto& pop : best_pop)
        {
            for (size_t i = 0; i < N_PROV; i++)
            {
                pop[i] = recv_best_vec[counter * N_PROV + i];
            }
            counter++;
        }

        best_pop.sort_population(distance_matrix);
        auto best = *(best_pop.begin());

        std::printf("Rank %d => Best Calculated\nStart Writing\n", rank);
        // std::cout << "Rank :" << std::setw(3) << rank << ": Best Calculated\n"
        //           << "Start Writing\n";

        std::ofstream foff("best_route.csv");

        if (!foff)
        {
            std::cerr << "Rank :" << std::setw(3) << rank
                      << "Could not open file, Aborting\n";
            exit(-3);
        }

        foff << "INDEX,PROVINCE,LONGITUDE,LATITUDE\n";
        for (size_t i = 0; i < N_PROV; i++)
        {
            const uint16_t idx = static_cast<uint16_t>(best[i]);
            const std::string name = prov_name[idx];
            const double longitude = prov_pos[idx](0);
            const double latitude = prov_pos[idx](1);
            foff << idx << ',' << name << ',' << longitude << ',' << latitude << '\n';
        }
        foff.close();
        std::printf(
            "Rank %d => Absolute best is : %f\n", rank, best.cost(distance_matrix));
        // std::cout << "Rank : " << std::setw(3) << rank
        //           << " : Absolute best is : " << best.cost(distance_matrix) << '\n';
    }

    MPI_Finalize();
    return 0;
}

void Evolution_Parameters::print(const int rank) const
{
    std::printf("Rank %d => Population Size is %lu\n", rank, POP_SIZE);
    std::printf("Rank %d => Number of generation to evolve : %lu\n", rank, N_GEN);
    std::printf(
        "Rank %d => Number of contacts between continents : %lu\n", rank, N_CONTACTS);
    std::printf("Rank %d => Number of elements to transfer for continent : %lu\n",
                rank,
                N_TRANSFER_ELEMENTS);
    std::printf("Rank %d => Number of shuffle to generate population : %u\n",
                rank,
                initial_shuffling);
    std::printf("Rank %d => Selection Coefficent : %.2f\n", rank, selection_coeff);
}

void Mutation_Probabilities::print(const int rank) const
{
    std::printf(
        "Rank %d => Probability of crossover : %.2f\n", rank, crossover_prob * 100.);
    std::printf("Rank %d => Probability of swap : %.2f\n", rank, swap_prob * 100.);
    std::printf(
        "Rank %d => Probability of permutation : %.2f\n", rank, permutate_prob * 100.);
    std::printf(
        "Rank %d => Probability of Inversion : %.2f\n", rank, inverse_prob * 100.);
    std::printf("Rank %d => Probability of shift : %.2f\n", rank, shift_prob * 100.);
}

arma::mat load_distance_matrix(const int rank,
                               std::vector<arma::vec2>& prov_pos,
                               std::vector<std::string>& prov_name)
{
    arma::mat distance_matrix;
    if (rank == MASTER)
    {
        prov_pos = load_province_position("cap_prov_ita.dat");

        prov_name = load_province_name("prov_ita.txt");

        if (prov_name.size() != N_PROV)
        {
            std::cerr << "Error in number of provinces, Aborting\n"
                      << "Detected Number : " << prov_name.size() << '\n'
                      << "Expected : " << N_PROV << '\n';

            exit(-1);
        }

        distance_matrix = create_matrix(prov_pos);
        std::printf("Rank %d => Matrix Shape is : %llu x %llu\n",
                    rank,
                    distance_matrix.n_rows,
                    distance_matrix.n_cols);
    }
    else
    {
        distance_matrix.zeros(N_PROV, N_PROV);
    }

    return distance_matrix;
}

void communicate_distance_matrix(const int rank, arma::mat& distance_matrix)
{ // Creating a scope to have a temporary vector to be killed after going out of scope
    std::printf("Rank %d => Before Transfer %f\n", rank, distance_matrix(1, 0));

    std::vector<double> data_transfer(N_PROV * N_PROV, 0.);
    if (rank == MASTER)
    {
        for (unsigned i = 0; i < N_PROV; i++)
        {
            for (unsigned j = 0; j < N_PROV; j++)
            {
                data_transfer[i * N_PROV + j] = distance_matrix(i, j);
            }
        }
    }

    MPI_Bcast(
        data_transfer.data(), data_transfer.size(), MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    if (rank != MASTER)
    {
        for (unsigned i = 0; i < N_PROV; i++)
        {
            for (unsigned j = 0; j < N_PROV; j++)
            {
                distance_matrix(i, j) = data_transfer[i * N_PROV + j];
            }
        }
    }
}
