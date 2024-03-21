#include <numeric>
#include <vector>
#include <array>

#include "random.h"
#include "vector3D.h"

class Random_Walk
{
private:
    double m_a{1.};
    Vector3D m_start_pos{0.0, 0.0, 0.0};
    Vector3D m_actual_pos = m_start_pos;
    Random *m_rng{nullptr};

public:
    Random_Walk(Random *rng) : m_rng{rng} { ; };
    Random_Walk(Vector3D start_pos, Random *rng) : m_start_pos{start_pos}, m_actual_pos{start_pos}, m_rng{rng} { ; };
    Random_Walk(double a, Vector3D start_pos, Random *rng) : m_a{a}, m_start_pos{start_pos}, m_actual_pos{start_pos}, m_rng{rng} { ; };
    ~Random_Walk() = default;

    Vector3D get_actual_pos() const { return m_actual_pos; }
    void reset() { m_actual_pos = m_start_pos; }
    void set_start_pos(Vector3D pos) { m_start_pos = pos; }

    double get_distance_from_origin() const { return m_actual_pos.distance(m_start_pos); }
    double get_distance_squared_from_origin() const { return m_actual_pos.distance_squared(m_start_pos); }

    void discrete_single_walk();
    void discrete_walk(const unsigned int steps);

    void continuos_single_walk();
    void continuos_walk(const unsigned int steps);
};
