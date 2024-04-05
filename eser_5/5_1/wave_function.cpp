#include "wave_function.h"

double ground_state(Vector3D pos)
{
    double probability = std::exp(-pos.distance() * 2.);
    return probability;
}

double first_excited(Vector3D pos)
{
    double probability = pos.z * pos.z * std::exp(-pos.distance());
    return probability;
}
