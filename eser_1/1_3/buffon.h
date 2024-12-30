#ifndef __BUFFON__
#define __BUFFON__

#include "random.h"
#include <cmath>
#include <iostream>
#include <limits>

class Buffon
{
private:
    double m_L{1.}, m_d{10.};
    int m_n_section{10}; // Number of grids excluded the first in position
    int m_throws{0}, m_intersect{0};
    Random *m_rng{NULL};

    /*
     *  |      |      |      |      |      |      |
     *  |      |      |      |      |      |      |
     *  |      |      |      |      |      |      |
     *  |  +---|-->   |      |      |      |      |
     *  |      |      |      |      |      |      |
     *  |      |      |      |      |      |      |
     *  |      |      |      |      |      |      |
     *  |      |      |      |      |      |      |
        ^      ^      ^      ^      ^      ^      ^
     x  0      1      2      3    ...     ...    n_section
     */

public:
    Buffon() = default;
    Buffon(double L, double d, int n_section, Random *rng) : m_L{L}, m_d{d}, m_n_section{n_section}, m_rng(rng)
    {
        ;
    }

    // ~Buffon() { delete m_rng; }

    void set_rng(Random *rng)
    {
        m_rng = rng;
    }
    int get_intersection() const
    {
        return m_intersect;
    }

    /// @brief Use Random Number generator to extract starting position form 0 to n_section * d
    /// @return starting position form 0 to n_section * d
    double starting_position();

    /// @brief Extract and angle from 0 to 2 PI to calculate needle ending position as start_pos + L * cos(theta)
    /// @param start_pos starting position
    /// @return ending position of the needle
    double ending_position(double start_pos);

    /// @brief check if the needle is intersecting a grid
    /// @param start needle starting position
    /// @param end needle ending position
    /// @return if needle is intersecting
    bool is_intersecting(double start, double end) const;

    /// @brief launch the needle n_throws times and update internal state to number of intersection obtained
    /// @param n_throws number of throws to execute
    void launch(int n_throws);

    /// @brief Compute pi as (2 * L * n_throws) / (n_intersection * d)
    /// @return An approximation of pi
    double compute_pi();
};

#endif