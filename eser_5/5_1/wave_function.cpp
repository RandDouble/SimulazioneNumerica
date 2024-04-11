#include "wave_function.h"

double ground_state(Vector3D pos)
{
    double probability = std::exp(-2. * pos.distance());
    return probability;
}

double first_excited(Vector3D pos)
{
    double probability = pos.z * pos.z * std::exp(-pos.distance());
    return probability;
}

double second_excited(Vector3D pos)
{
    double probability = std::pow(2 * pos.z * pos.z - pos.x * pos.x - pos.y * pos.y, 2) * std::exp(-pos.distance()) ;
    return probability;
}
