#include <algorithm>
#include <cstdlib>
#include <initializer_list>
#include <vector>

#include <armadillo>

#include "initializer.h"
#include "random.h"
#include "swap_functions.h"

#ifndef __MUTATIONS__
#define __MUTATIONS__

template <size_t SIZE, typename GenomeType = uint8_t> class Individual : public std::vector<GenomeType>
{
    static_assert(std::is_integral<GenomeType>::value, "GenomeType must be an integral type");

public:
    using individual = Individual<SIZE, GenomeType>;
    using vector = std::vector<GenomeType>;

    Individual() : vector(SIZE, GenomeType{})
    {
        std::iota(individual::begin(), individual::end(), GenomeType{});
    }

    Individual(individual &other) : vector(SIZE, GenomeType{})
    {
        std::copy(other.begin(), other.end(), individual::begin());
    }

    Individual(individual &&other) : vector(other)
    {
        assert(individual::size() == SIZE && "Array is malformed\n");
    }

    Individual(std::initializer_list<GenomeType> l) : vector( l )
    {
        assert(check_health() && "You have duplicates in the initializer list");
    }

    // template <typename T/* , typename = std::enable_if_t<std::is_integral_v<T>> */> Individual(T &&l...) :
    // vector(l...)
    // {
    //     ;
    // }


    individual &operator=(const individual &other)
    {
        if (this == &other)
            return *this;
        if (individual::size() != other.size())
            individual::resize(SIZE);
        assert(individual::size() == SIZE && "Array has size different from individual\n");
        assert(other.size() == individual::size() && "You are trying to assign a different size array\n");
        std::copy(other.begin(), other.end(), individual::begin());
        return *this;
    }


    void pair_permutation(Random &rng)
    {
        // Choosing random indexes

        std::size_t idx_1 = rng.Ranint(1, SIZE);
        std::size_t idx_2 = PBC(SIZE - 1, rng.Ranint(idx_1 + 1, idx_1 + SIZE)); // Choosing between SIZE - 1 elements
                                                                                // differerent from idx_1

        // Using Swap
        std::swap(individual::operator[](idx_1), individual::operator[](idx_2));
    }

    void shift_block(Random &rng)
    {
        std::size_t start_idx = rng.Ranint(1, SIZE);
        std::size_t middle_idx = start_idx + rng.Ranint(1, SIZE - 1);
        std::size_t end_idx = middle_idx + rng.Ranint(1, SIZE - (middle_idx - start_idx));

        assert((end_idx - middle_idx > 0) && "We want to shift only on the left.\n");
        assert((end_idx - middle_idx < SIZE) && "Shift size should be less than array size");

        // Moving to the end block to move
        // const auto begin = m_DNA.begin() + starting_idx;
        // const auto middle = begin + block_size;
        // const auto end = begin + block_size + shift_size;

        // std::rotate(begin, middle, end);

        PBC_swap::rotate<SIZE>(individual::begin(), start_idx, middle_idx, end_idx, 1);
    }

    void permutate_contiguos(Random &rng)
    {
        std::size_t m_contiguos = rng.Ranint(1, SIZE / 2);
        std::size_t first_pos_idx = rng.Ranint(1, SIZE); // Selected in the first half - m_contiguos
        std::size_t second_pos_idx =
            first_pos_idx + rng.Ranint(m_contiguos, SIZE - m_contiguos); // Select in second half - m_contiguos
        // auto first_position = m_DNA.begin() + first_pos_idx;
        // auto second_position = m_DNA.begin() + second_pos_idx;

        // std::swap_ranges(first_position, first_position + m_contiguos, second_position);
        PBC_swap::swap_ranges<SIZE>(individual::begin(), first_pos_idx, second_pos_idx, m_contiguos, 1);
    }

    void inversion(Random &rng)
    {
        auto inversion_lenght = rng.Ranint(2, SIZE);
        auto idx_start_inversion = rng.Ranint(1, SIZE);

        assert((inversion_lenght < SIZE) && "You choose an invertion lenght higher than array lenght");
        assert((inversion_lenght != 0ull) && "What sense has to have an invertion lenght equals to zero");
        // assert(((idx_start_inversion + inversion_lenght) <= SIZE) && "You cannot go out of array bounds");

        // auto start = m_DNA.begin() + idx_start_inversion;
        // auto end = start + inversion_lenght; // Last element is not included in rotation
        // std::reverse(start, end);
        PBC_swap::reverse<SIZE>(individual::begin(), idx_start_inversion, idx_start_inversion + inversion_lenght, 1);
    }

    void crossover(Individual &mother, Individual &daughter, Individual &son, Random &rng)
    {
#ifdef TEST_ENV
        auto cut_position = 3;
#else
        auto cut_position = rng.Ranint(1, SIZE);
#endif
        std::copy(individual::begin(), individual::begin() + cut_position, son.begin());
        std::copy(mother.begin(), mother.begin() + cut_position, daughter.begin());

        std::size_t counter_son = cut_position;
        std::size_t counter_daughter = cut_position;
        // Filling son with mother missing parts
        for (std::size_t j = 1; j < SIZE; j++)
        {
            for (std::size_t k = cut_position; k < SIZE; k++)
            {
                if (individual::operator[](k) ==
                    mother[j]) // Father DNA in k position is equal to mother's DNA at j position
                {
                    son[counter_son++] = mother[j]; // Mother's DNA is copied in son.
                }
                if (mother[k] ==
                    individual::operator[](j)) // Father DNA in j position is equal to mother's DNA at k position
                {
                    daughter[counter_daughter++] = individual::operator[](j); // Father DNA is copied in daughter's DNA
                }
            }
        }
    }

    void print_DNA() const
    {
        print_vector<GenomeType> print;
        print(*this);
    }

    bool check_health()
    {
        bool result = true;

        result = result && (individual::operator[](0) == 0);

        vector copy_vector(SIZE, 0);

        for (size_t i = 0; i < SIZE; i++)
        {
            if(copy_vector[individual::operator[](i)] == 0 )
                copy_vector[individual::operator[](i)] += 1;
            else
                result = false;
        }
        if (!result)
        {
            std::cout << "There are duplicates in the DNA\n";    
        }

        // std::copy(individual::begin(), individual::end(), copy_vector.begin());
        // std::sort(copy_vector.begin(), copy_vector.end());

        // auto end_pos = std::unique(copy_vector.begin(), copy_vector.end());
        // result = result && (SIZE == std::distance(copy_vector.begin(), end_pos));

        return result;
    }

    double cost(const std::array<arma::vec2, SIZE> &positions) const
    {
        double acc = 0.;
        for (size_t i = 0; i < SIZE; i++)
        {
            arma::vec2 difference = positions[individual::operator[](PBC_swap::PBC<SIZE>(i + 1))] -
                                    positions[individual::operator[](PBC_swap::PBC<SIZE>(i))];
            // std::cout << '(' << PBC_swap::PBC<SIZE>(i)<< ',' << PBC_swap::PBC<SIZE>(i+1) << ")  ";
            acc += arma::norm(difference);
        }
        // std::cout << '\n';
        return acc;
    }
};

#endif // __MUTATIONS__