#include "random_walk.h"


void Random_Walk::discrete_single_walk()
{
    int dir = static_cast<int>(std::floor(m_rng->Rannyu() * 3)) % 3;
    switch (dir)
    {
    case 0: // x direction
        m_actual_pos.x += (m_rng->Rannyu() > 0.5) ? m_a : -m_a;
        break;
    case 1: // y direction
        m_actual_pos.y += (m_rng->Rannyu() > 0.5) ? m_a : -m_a;
        break;
    case 2: // z direction
        m_actual_pos.z += (m_rng->Rannyu() > 0.5) ? m_a : -m_a;
        break;
    default:
        std::cerr << "Cannot be another direction, error, exiting\n";
        exit(-1);
        break;
    }
}

void Random_Walk::discrete_walk(const unsigned int steps)
{
    for (unsigned int i = 0; i < steps; i++)
    {
        discrete_single_walk();
    }
}

void Random_Walk::continuos_single_walk()
{
    m_actual_pos += Vector3D::generate_vector(m_rng) * m_a;
}

void Random_Walk::continuos_walk(const unsigned int steps)
{
    for (unsigned int i = 0; i < steps; i++)
    {
        continuos_single_walk();
    }
}

std::ostream &operator<<(std::ostream &oss, Vector3D vec)
{
    oss << vec.x << ',' << vec.y << ',' << vec.z;
    return oss;
}   
