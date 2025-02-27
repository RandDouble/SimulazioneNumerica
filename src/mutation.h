#include <algorithm>
#include <array>
#include <cstdlib>
#include <initializer_list>

#include <armadillo>

#include "initializer.h"
#include "random.h"
#include "swap_functions.h"

#ifndef __MUTATIONS__
#define __MUTATIONS__

template <std::size_t SIZE, typename POP_TYPE = uint8_t> class Individual
{
private:
    std::array<POP_TYPE, SIZE> m_DNA;

public:
    Individual()
    {
        std::iota(m_DNA.begin(), m_DNA.end(), 0);
    }

    Individual(Individual<SIZE, POP_TYPE> &other)
    {
        std::copy(other.m_DNA.begin(), other.m_DNA.end(), m_DNA.begin());
    }
    Individual(Individual<SIZE, POP_TYPE> &&other) : m_DNA(std::move(other.m_DNA))
    {
        ;
    }

#ifdef ENABLE_INDIVIDUAL_INITIALIZER_LIST
    template <typename... T> Individual(T &&...l) : m_DNA{{static_cast<uint8_t>(std::forward<T>(l))...}}
    {
        ;
    }
#endif

    Individual<SIZE, POP_TYPE> &operator=(Individual<SIZE, POP_TYPE> &other)
    {
        std::copy(other.m_DNA.begin(), other.m_DNA.end(), m_DNA.begin());
        return *this;
    }

    Individual<SIZE, POP_TYPE> &operator=(Individual<SIZE, POP_TYPE> &&other)
    {
        m_DNA = std::move(other.m_DNA);
        return *this;
    }

    constexpr decltype(m_DNA.begin()) begin()
    {
        return m_DNA.begin();
    }
    constexpr decltype(m_DNA.end()) end()
    {
        return m_DNA.end();
    }

    constexpr decltype(m_DNA.cbegin()) begin() const
    {
        return m_DNA.cbegin();
    }
    constexpr decltype(m_DNA.cend()) end() const
    {
        return m_DNA.cend();
    }

    constexpr decltype(m_DNA.cbegin()) cbegin() const
    {
        return m_DNA.cbegin();
    }
    constexpr decltype(m_DNA.cend()) cend() const
    {
        return m_DNA.cend();
    }

    POP_TYPE operator[](const std::size_t idx) const
    {
        return m_DNA[idx];
    }
    POP_TYPE &operator[](const std::size_t idx)
    {
        return m_DNA[idx];
    }

    constexpr std::size_t size() const
    {
        return m_DNA.size();
    }
    decltype(m_DNA.data()) data()
    {
        return m_DNA.data();
    }

    void pair_permutation(Random &rng)
    {
        // Choosing random indexes

        std::size_t idx_1 = rng.Ranint(1, SIZE);
        std::size_t idx_2 = PBC(SIZE - 1, rng.Ranint(idx_1 + 1, idx_1 + SIZE)); // Choosing between SIZE - 1 elements
                                                                                // differerent from idx_1

        // Using Swap
        std::swap(m_DNA[idx_1], m_DNA[idx_2]);
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

        PBC_swap::rotate<SIZE>(m_DNA.begin(), start_idx, middle_idx, end_idx, 1);
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
        PBC_swap::swap_ranges<SIZE>(m_DNA.begin(), first_pos_idx, second_pos_idx, m_contiguos, 1);
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
        PBC_swap::reverse<SIZE>(m_DNA.begin(), idx_start_inversion, idx_start_inversion + inversion_lenght, 1);
    }

    void crossover(Individual &mother, Individual &daughter, Individual &son, Random &rng)
    {
#ifdef TEST_ENV
        auto cut_position = 3;
#else
        auto cut_position = rng.Ranint(1, SIZE);
#endif
        assert(check_health() && "Father ha genetic defects");
        assert(mother.check_health() && "Mother has genetic defects");

        std::copy(m_DNA.begin(), m_DNA.begin() + cut_position, son.begin());
        std::copy(mother.begin(), mother.begin() + cut_position, daughter.begin());

        std::size_t counter_son = cut_position;
        std::size_t counter_daughter = cut_position;
        // Filling son with mother missing parts
        for (std::size_t j = 1; j < SIZE; j++)
        {
            for (std::size_t k = cut_position; k < SIZE; k++)
            {
                if (m_DNA[k] == mother[j]) // Father DNA in k position is equal to mother's DNA at j position
                {
                    son[counter_son++] = mother[j]; // Mother's DNA is copied in son.
                }
                if (mother[k] == m_DNA[j]) // Father DNA in j position is equal to mother's DNA at k position
                {
                    daughter[counter_daughter++] = m_DNA[j]; // Father DNA is copied in daughter's DNA
                }
            }
        }
        assert(daughter.check_health() && "Daughter dead after birth");
        assert(son.check_health() && "Son dead after birth");
    }

    void print_DNA() const
    {
        print_vector<POP_TYPE> print;
        print(this->m_DNA);
    }


    /// @brief Check health of the Individual
    /// @return False if the individual is ill, True otherwise
    bool check_health() const
    {
        bool result = true;

        result = result && (m_DNA[0] == 0);

        std::array<POP_TYPE, SIZE> copy_vector;
        std::copy(begin(), end(), copy_vector.begin());
        std::sort(copy_vector.begin(), copy_vector.end());

        auto end_pos = std::unique(copy_vector.begin(), copy_vector.end());

        result = result && (SIZE == std::distance(copy_vector.begin(), end_pos));

        return result;
    }

    double cost(const std::array<arma::vec2, SIZE> &positions) const
    {
        double acc = 0.;
        assert(check_health() && "This vector is ill\n");
        for (size_t i = 0; i < SIZE; i++)
        {
            arma::vec2 difference =
                positions[m_DNA[PBC_swap::PBC<SIZE>(i + 1)]] - positions[m_DNA[PBC_swap::PBC<SIZE>(i)]];
            // std::cout << '(' << PBC_swap::PBC<SIZE>(i)<< ',' << PBC_swap::PBC<SIZE>(i+1) << ")  ";
            acc += arma::norm(difference);
        }
        // std::cout << '\n';
        return acc;
    }

    double cost(const std::vector<arma::vec2> &positions) const
    {
        assert(check_health() && "This vector is ill\n");

        double acc = 0.;
        for (size_t i = 0; i < SIZE; i++)
        {
            arma::vec2 difference =
                positions[m_DNA[PBC_swap::PBC<SIZE>(i + 1)]] - positions[m_DNA[PBC_swap::PBC<SIZE>(i)]];
            // std::cout << '(' << PBC_swap::PBC<SIZE>(i)<< ',' << PBC_swap::PBC<SIZE>(i+1) << ")  ";
            acc += arma::norm(difference);
        }
        // std::cout << '\n';
        return acc;
    }

    double cost(const arma::mat &matrix) const
    {
        assert(check_health() && "This vector is ill\n");

        double acc = 0.;
        for (size_t i = 0; i < SIZE; i++)
        {
            // std::cout << "Actual i : " << PBC_swap::PBC<SIZE>(i) << "\tActual i+1 : " << PBC_swap::PBC<SIZE>(i + 1)
            // << '\n'
            //           << "DNA i-th : " << m_DNA[PBC_swap::PBC<SIZE>(i)] << "\tDNA i+1-th" <<
            //           m_DNA[PBC_swap::PBC<SIZE>(i + 1)] << '\n';
            try
            {
                const double dist = matrix(m_DNA[PBC_swap::PBC<SIZE>(i + 1)], m_DNA[PBC_swap::PBC<SIZE>(i)]);
                acc += dist;
            }
            catch (const std::out_of_range &range)
            {
                if constexpr (std::is_same<POP_TYPE, uint8_t>::value)
                {
                    std::for_each(m_DNA.begin(), m_DNA.end(),
                                  [](const uint8_t &i) { std::cerr << std::setw(4) << static_cast<uint16_t>(i); });
                    std::cerr << "\nActual i, i+1 : " << std::setw(4)
                              << static_cast<uint16_t>(m_DNA[PBC_swap::PBC<SIZE>(i)]) << std::setw(4)
                              << static_cast<uint16_t>(m_DNA[PBC_swap::PBC<SIZE>(i + 1)]) << '\n'
                              << "Matrix size : " << arma::size(matrix) << '\n'
                              << range.what() << '\n';
                }
                else
                {
                    std::for_each(m_DNA.begin(), m_DNA.end(),
                                  [](const uint8_t &i) { std::cerr << std::setw(4) << static_cast<POP_TYPE>(i); });
                    std::cerr << "\nActual i, i+1 : " << std::setw(4)
                              << m_DNA[PBC_swap::PBC<SIZE>(i)] << std::setw(4)
                              << m_DNA[PBC_swap::PBC<SIZE>(i + 1)] << '\n'
                              << "Matrix size : " << arma::size(matrix) << '\n'
                              << range.what() << '\n';
                }
                exit(-2);
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << '\n';
            }

            // std::cout << '(' << PBC_swap::PBC<SIZE>(i)<< ',' << PBC_swap::PBC<SIZE>(i+1) << ")  ";
        }
        // std::cout << '\n';
        return acc;
    }
};

#endif // __MUTATIONS__