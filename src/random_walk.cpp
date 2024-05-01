#include "random_walk.h"

void Random_Walk::discrete_single_walk()
{
    int dir = static_cast<int>(m_rng->Rannyu() * 6.);
    switch (dir)
    {
    case 0: // x direction
        m_actual_pos.x += m_a;
        break;
    case 1: // y direction
        m_actual_pos.y += m_a;
        break;
    case 2: // z direction
        m_actual_pos.z += m_a;
        break;
    case 3:
        m_actual_pos.x -= m_a;
        break;
    case 4:
        m_actual_pos.y -= m_a;
        break;
    case 5:
        m_actual_pos.z -= m_a;
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
