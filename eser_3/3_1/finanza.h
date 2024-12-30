#ifndef __FINANZA__
#define __FINANZA__

#define NDEBUG

#include <numeric>

#include "random.h"
#include "utilities.h"

inline double Wiener(Random *rng, double time)
{
    return rng->Gauss(0., time);
}

inline double Wiener_step(Random *rng)
{
    return rng->Gauss(0., 1.);
}

class FinancialOption
{

    const double S_0{100.};
    const double m_r{0.1}, m_sigma{.25};
    const double m_K{100.};
    double m_T{1.};
    Random *m_rng{nullptr};
    unsigned int m_extraction{100};

public:
    FinancialOption(Random *rng) : m_rng{rng}
    {
        ;
    }

    void set_extraction(unsigned int n_extraction)
    {
        m_extraction = n_extraction;
    }
    unsigned int get_extraction()
    {
        return m_extraction;
    }

    double call_option_price_conv(double ending_price) const;
    double put_option_price_conv(double ending_price) const;

    double call_option_final_mean();
    double put_option_final_mean();
    double call_option_step_mean(const unsigned int n_step);
    double put_option_step_mean(const unsigned int n_step);

    double ending_price_final();
    double ending_price_step(const unsigned int n_step);
};

#endif // __FINANZA__