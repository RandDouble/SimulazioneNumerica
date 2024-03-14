#include <numeric>
#include <vector>
#include <array>

#include "random.h"

class random_walk
{
private:
    double m_a = 1.;
    std::array<double, 3> m_start_pos{0.0, 0.0, 0.0};
    std::array<double, 3> m_actual_pos = m_start_pos;
    Random *m_rng{nullptr};

public:
    random_walk(Random *rng) : m_rng{rng} { ; };
    ~random_walk() = default;

    void discrete_walk(const unsigned int steps);
    void continuos_walk(const unsigned int step, )
};
