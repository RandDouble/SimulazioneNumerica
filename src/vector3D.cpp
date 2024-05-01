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
    double phi = rng->Rannyu() * 2 * M_PI;
    double theta = std::acos(1. - 2. * rng->Rannyu());

    double x = std::sin(theta) * std::cos(phi);
    double y = std::sin(theta) * std::sin(phi);
    double z = std::cos(theta);

    return Vector3D{x, y, z};
}

Vector3D Vector3D::generate_unif(Random *rng, const double delta)
{
    double x{rng->Rannyu(-delta, delta)}, y{rng->Rannyu(-delta, delta)}, z{rng->Rannyu(-delta, delta)};
    return Vector3D{x, y, z};
}

Vector3D Vector3D::generate_gauss(Random *rng, const double delta)
{
    double x{rng->Gauss(0, delta)}, y{rng->Gauss(0, delta)}, z{rng->Gauss(0, delta)};
    return Vector3D{x, y, z};
}

std::ostream &operator<<(std::ostream &oss, const Vector3D &vec)
{
    oss << vec.x << ',' << vec.y << ',' << vec.z;
    return oss;
}
