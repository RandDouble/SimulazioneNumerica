#include <cmath>
#include <vector>

#ifndef __FUNCTION__
#define __FUNCTION__

double integrand(double x)
{
    return M_PI_2 * std::cos(M_PI_2 * x);
}

#endif //__FUNCTION__