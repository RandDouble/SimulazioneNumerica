#include <cmath>
#include <vector>

#ifndef __FUNCTION__
#define __FUNCTION__

inline double integrand(double x)
{
    return M_PI_2 * std::cos(M_PI_2 * x);
}

inline double integrand_anti(double x)
{
    double y = M_PI_2 / std::sqrt(2) * std::cos(M_PI_2 * (x - 0.5));
    return y;
}

double integrand_anti_expansion_norm(double x)
{
    constexpr double denom = 1. - (M_PI_4 * M_PI_4 / 6.);
    constexpr double other_piece = M_PI_2 * M_PI_4; // pi^2 / 8
    double y = (1 - other_piece * (x - 0.5) * (x - 0.5)) / denom;
    return y;
}

constexpr inline double integrand_anti_expansion_max()
{
    return (1 / (1. - (M_PI_4 * M_PI_4 / 6.)));
}

double integrand_precise(double x)
{
    return integrand_anti(x) / integrand_anti_expansion_norm(x);
}

#endif //__FUNCTION__