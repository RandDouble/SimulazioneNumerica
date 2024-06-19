#include "mpi.h"

#include <cassert>
#include <fstream>
#include <iostream>

#define ARMA_PRINT_EXCEPTIONS
#define ARMA_PRINT_ERRORS

#include "initializer.h"
#include "population.h"

#include "random.h"
#include "utilities.h"

#define MASTER 0

// I Hate a little bit to use Macro Constants, they are error-prone, generally unsafe
// and have no type or compile checking, but in this case they can be useful..

#define POP_SIZE 200ull
#define N_GEN 300ull
#define N_PROV 110
#define N_CONTACTS 15ull
#define N_TRANSFER_ELEMENTS 5ull

constexpr unsigned int initial_shuffling = 300;

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

    std::cout << "Hello World!!\n"
              << "I am rank : " << rank << " of size : " << size << '\n'
              << "on CPU : " << hostname << '\n';

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
        std::cout << "Rank : " << std::setw(3) << " : Matrix Shape is : " << arma::size(distance_matrix) << '\n';
    } // Quando capirai come fare la comunicazione tra i vari thread possiamo spostare l'inizializzazione al primo thread e tutti gli altri aspettano
    else
    {
        distance_matrix.zeros(N_PROV, N_PROV);
    }

    { // Creating a scope to have a temporary vector to be killed after going out of scope
        std::cout << "Rank :" << std::setw(3) << rank << " : Before Transfer : " << distance_matrix(1, 0) << '\n';

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

        MPI_Bcast(data_transfer.data(), data_transfer.size(), MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

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
        std::cout << "Rank :" << std::setw(3) << rank << " : After Transfer : " << distance_matrix(1, 0) << '\n';

        std::cout << "Rank :" << std::setw(3) << rank << " : Successfull data transfer\n";
    }

    // Generating Random Number Generator
    Random rng;
    initializer(rng, rank); // selecting different rows based on rank to have indipendent RNG

    // Initializing Population... Due to implementation this is not a generic algorithm, we need to know size... But we know Italian provinces are 110
    Population<N_PROV> population(POP_SIZE);

    // Checking for health
    for (auto& pop : population)
    {
        assert(pop.check_health() && "Found ill vector\n");
    }

    // Initial shuffling
    for (auto& element : population)
    {
        for (unsigned int i = 0; i < initial_shuffling; i++)
        {
            element.pair_permutation(rng);
        }
    }

    std::cout << "Rank : " << rank << ", Population initialized, RNG initialized\n"
              << "Starting Evolution\n";

    // Evolution
    for (size_t n = 0; n < N_CONTACTS; n++)
    {
        for (size_t i = 0; i < N_GEN; i++)
        {
            // std::cout << "Rank : " << rank << "\tGeneration nr :" << std::setw(4) << i << '\n';

            population.sort_population(distance_matrix);
            population.new_gen(rng);
        }

        population.sort_population(distance_matrix);

        auto initial_pos = population.begin();
        std::vector<uint8_t> best_el(initial_pos->size() * N_TRANSFER_ELEMENTS);

        for (size_t i = 0; i < N_TRANSFER_ELEMENTS; i++)
        {
            for (size_t j = 0; j < initial_pos->size(); j++)
            {
                best_el[i * N_PROV + j] = (*(initial_pos + i))[j];
            }
        }

        std::vector<uint8_t> receiver(initial_pos->size() * N_TRANSFER_ELEMENTS * size);

        // Contact between continents
        MPI_Allgather(best_el.data(), best_el.size(), MPI_UINT8_T, receiver.data(), best_el.size(), MPI_UINT8_T, MPI_COMM_WORLD);

        // Checking if copy was happening
        // std::for_each(receiver.begin(), receiver.end(), [](const uint8_t& i)
        //                              { std::cout << std::setw(4) << static_cast<uint16_t>(i); });

        auto copying_pos = population.end() - (N_TRANSFER_ELEMENTS * (size - 1));
        size_t counter = 0;

        for (size_t provenience = 0; provenience < static_cast<size_t>(size); provenience++)
        {
            if (provenience == static_cast<size_t>(rank)) // Do not copy your datas
            {
                continue;
            }

            // copy elements from Individual i-th come from provenience-th rank
            for (size_t i = 0; i < N_TRANSFER_ELEMENTS; i++)
            {
                for (size_t j = 0; j < copying_pos->size(); j++)
                {
                    (*(copying_pos + counter * N_TRANSFER_ELEMENTS + i))[j] = receiver[counter * N_PROV * N_TRANSFER_ELEMENTS + i * N_PROV + j];
                }
            }
            counter++;
        }

        population.sort_population(distance_matrix);
        // std::cout << '\n'
        //           << "Transferring Complete\n";
    }
    std::cout << "Rank : " << rank << "\tEvolution Ended\n"
              << "\tBest is : " << population.begin()->cost(distance_matrix) << '\n';

    // Getting Best of each continent comparing and then outputting best of all

    Population<N_PROV> best_pop((rank == MASTER) ? size : 0);
    std::vector<uint8_t> recv_best_vec((rank == MASTER) ? size * N_PROV : 0, 0);

    std::cout << "Rank :" << std::setw(3) << rank << (rank == MASTER) ? " : Asking for other best\n" : "Sending my best\n";

    MPI_Gather(population.begin()->data(), N_PROV, MPI_UINT8_T, recv_best_vec.data(), population.begin()->size(), MPI_UINT8_T, MASTER, MPI_COMM_WORLD);

    if (rank == MASTER)
    {
        std::cout << "Rank :" << std::setw(3) << rank << " : Communication Succeed\n"
                  << "Calculating best\n";

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

        std::cout << "Rank :" << std::setw(3) << rank << ": Best Calculated\n"
                  << "Start Writing\n";

        std::ofstream foff("best_route.csv");

        if (!foff)
        {
            std::cerr << "Rank :" << std::setw(3) << rank << "Could not open file, Aborting\n";
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
        std::cout << "Rank : " << std::setw(3) << rank << " : Absolute best is : " << best.cost(distance_matrix) << '\n';
    }

    MPI_Finalize();
    return 0;
}
