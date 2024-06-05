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
    double mutation_prob, crossover_prob;

public:
    Population(const std::size_t size) : m_pop(size), m_old_gen(size) { ; }
    Population(std::vector<Individual<SIZE>>& vec) : m_pop{vec}, m_old_gen{vec} { ; }
    Population(std::vector<Individual<SIZE>>&& vec) : m_pop{vec} { m_old_gen = std::copy(m_pop); }

    constexpr decltype(m_pop.begin()) begin() { return m_pop.begin(); }
    constexpr decltype(m_pop.end()) end() { return m_pop.end(); }
    constexpr decltype(m_pop.cbegin()) cbegin() const { return m_pop.cbegin(); }
    constexpr decltype(m_pop.cend()) cend() const { return m_pop.cend(); }
    constexpr decltype(m_pop.cbegin()) begin() const { return cbegin(); }
    constexpr decltype(m_pop.cend()) end() const { return cend(); }

    void sort_population(std::array<arma::vec2, SIZE>& positions)
    {
        std::sort(m_pop.begin(), m_pop.end(), [](const Individual<SIZE>& a, const Individual<SIZE>& b)
                  { return a.cost(positions) < b.cost(positions); });
    }

    void new_gen()
    {
    }
};

#endif