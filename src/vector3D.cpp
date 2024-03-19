#include "vector3D.h"

Vector3D &Vector3D::operator+=(const Vector3D &rhs)
{
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
}

Vector3D Vector3D::generate_vector(Random *rng)
{
    double phi = rng->Rannyu() * M_PI * 2.;
    double theta = rng->ExternalInvCum([](double r)
                                       { return std::acos(1 - 2 * r); });
    return Vector3D{std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};
}