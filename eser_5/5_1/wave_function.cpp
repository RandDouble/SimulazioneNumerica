#include "wave_function.h"

double ground_state(const Vector3D &pos)
{
    const double probability = std::exp(-2. * pos.distance());
    return probability;
}

double first_excited(const Vector3D &pos)
{
    const double probability = pos.z * pos.z * std::exp(-pos.distance());
    return probability;
}

double second_excited(const Vector3D &pos)
{
    constexpr double factor = 2. / 3.;
    const double distance = pos.distance();
    const double probability = std::pow(3 * pos.z * pos.z - distance, 2) * std::exp(-factor * distance);
    return probability;
}
