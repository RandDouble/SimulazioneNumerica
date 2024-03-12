#include "buffon.h"

double Buffon::starting_position()
{
    double x = m_rng->Rannyu(0, m_n_section * m_d); // Number extraction between 0 and max_lenght
    return x;
}

double Buffon::ending_position(double start_pos)
{
    double x = 0., y = 0., norm = 0.;

    // Search if unit cirle a random vector, with accept reject
    do
    {
        x = m_rng->Rannyu(-1, 1);
        y = m_rng->Rannyu(-1, 1);
        norm = std::hypot(x, y);
    } while (norm > 1.);
    x /= norm;

    // double theta = m_rng->Rannyu(0, std::numeric_limits<double>::max());
    // x = std::cos(theta);

    double end_x = start_pos + m_L * x;
    return end_x;
}

bool Buffon::is_intersecting(double start, double end) const
{
    int first_sec = static_cast<int>(std::floor(start / m_d));
    int second_sec = static_cast<int>(std::floor(end / m_d));
    return (first_sec != second_sec) ? true : false;
}

void Buffon::launch(int launch)
{
    m_intersect = 0;
    m_throws = launch;
    for (int i = 0; i < m_throws; i++)
    {
        double start = starting_position();
        double end = ending_position(start);
        m_intersect += (is_intersecting(start, end)) ? 1 : 0;
    }
}

double Buffon::compute_pi()
{
    if (m_throws != 0)
    {
        return 2. * m_L * m_throws / (m_intersect * m_d);
    }
    else
    {
        std::cerr << "Launch needles before trying to calculate pi!!!\n";
        exit(-1);
    }
}
