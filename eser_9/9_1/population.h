#include <algorithm>
#include <vector>

#include "mutation.h"
#include "random.h"

#ifndef __POPULATION__
#define __POPULATION__

template <std::size_t SIZE>
class Population
{
private:
    std::vector<Individual<SIZE>> m_pop;
    std::vector<Individual<SIZE>> m_old_gen;
    double m_swap_prob{0.2}, m_shift_prob{0.1},
        m_permutate_prob{0.1}, m_invers_prob{0.1}, m_crossover_prob{0.6};
    double m_selection_coeff{-.5};

private:
    static double ICDF(double in, const double sel_coeff)
    {
        const double power = 1. / (1. + sel_coeff);
        return std::pow(in, power);
    }

    std::size_t selection_operator(Random& rng)
    {
        // std::function<double(double)> F = [&](double in) -> double
        // { return ICDF(in, m_selection_coeff); };

        //     return static_cast<size_t>(m_old_gen.size() *
        //            rng.ExternalInvCum(F));
        // }
        return rng.Ranint(0, m_pop.size() / 2);
    }

public:
    Population(const std::size_t size) : m_pop(size), m_old_gen(size) { ; }
    Population(std::vector<Individual<SIZE>>& vec) : m_pop{vec}, m_old_gen{vec} { ; }
    Population(std::vector<Individual<SIZE>>&& vec) : m_pop{vec} { m_old_gen = std::copy(m_pop); }

    void crossover_prob(const double& prob) { m_crossover_prob = prob; }
    void swap_prob(const double& prob) { m_swap_prob = prob; }
    void permutate_prob(const double& prob) { m_permutate_prob = prob; }
    void inverse_prob(const double& prob) { m_invers_prob = prob; }
    void shift_prop(const double& prob) { m_shift_prob = prob; }
    void selection_coeff(const double& coeff) { m_selection_coeff = coeff; }

    double crossover_prob() const { return m_crossover_prob; }
    double swap_prob() const { return m_swap_prob; }
    double permutate_prob() const { return m_permutate_prob; }
    double inverse_prob() const { return m_invers_prob; }
    double shift_prop() const { return m_shift_prob; }
    double selection_coeff() const { return m_selection_coeff; }

    constexpr decltype(m_pop.begin()) begin() { return m_pop.begin(); }
    constexpr decltype(m_pop.end()) end() { return m_pop.end(); }
    constexpr decltype(m_pop.cbegin()) cbegin() const { return m_pop.cbegin(); }
    constexpr decltype(m_pop.cend()) cend() const { return m_pop.cend(); }
    constexpr decltype(m_pop.cbegin()) begin() const { return cbegin(); }
    constexpr decltype(m_pop.cend()) end() const { return cend(); }

    void sort_population(std::array<arma::vec2, SIZE>& positions)
    {
        std::stable_sort(m_pop.begin(), m_pop.end(), [&positions](const Individual<SIZE>& a, const Individual<SIZE>& b)
                         { return (a.cost(positions) < b.cost(positions)); });
    }

    void new_gen(Random& rng)
    {
        m_old_gen = m_pop;
        for (size_t i = 0; i < m_pop.size(); i += 2)
        {

            size_t idx_mother = selection_operator(rng);
            size_t idx_father = selection_operator(rng);
            if (idx_mother == idx_father)
            {
                idx_father++;
            }
            Individual<SIZE> mother = m_old_gen[idx_mother];
            Individual<SIZE> father = m_old_gen[idx_father];
            Individual<SIZE>* son = &(m_pop[i]);
            Individual<SIZE>* daughter = &(m_pop[i + 1]);

            // Pair Permutation
            if (rng.Rannyu() < m_swap_prob)
                mother.pair_permutation(rng);
            if (rng.Rannyu() < m_swap_prob)
                father.pair_permutation(rng);

            // Shift Block
            if (rng.Rannyu() < m_shift_prob)
                mother.shift_block(rng);
            if (rng.Rannyu() < m_shift_prob)
                father.shift_block(rng);

            // Permutation
            if (rng.Rannyu() < m_permutate_prob)
                mother.permutate_contiguos(rng);
            if (rng.Rannyu() < m_permutate_prob)
                father.permutate_contiguos(rng);

            // Inversion
            if (rng.Rannyu() < m_invers_prob)
                mother.inversion(rng);
            if (rng.Rannyu() < m_invers_prob)
                father.inversion(rng);

            // Crossover
            if (rng.Rannyu() < m_crossover_prob)
            {
                mother.crossover(father, *son, *daughter, rng);
            }
            else
            {
                std::copy(mother.begin(), mother.end(), daughter->begin());
                std::copy(father.begin(), father.end(), son->begin());
            }
        }
    }
};

#endif