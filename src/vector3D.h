#ifndef __VECTOR_3D__
#define __VECTOR_3D__

#include <numeric>
#include <iostream>

#include "random.h"

struct Vector3D
{
    double x{0.}, y{0.}, z{0.};

    double distance(const Vector3D &ref = {0.0, 0.0, 0.0}) const { return std::hypot((x - ref.x), (y - ref.y), (z - ref.z)); }
    double distance_squared(const Vector3D &ref = {0.0, 0.0, 0.0}) const { return distance(ref) * distance(ref); }
    double operator*(const Vector3D &rhs) { return x * rhs.x + y * rhs.y + z * rhs.z; }
    Vector3D operator*(const double &rhs) { return Vector3D{x * rhs, y * rhs, z * rhs}; }
    Vector3D operator+(const Vector3D &rhs) const { return Vector3D{x + rhs.x, y + rhs.y, z + rhs.z}; }
    Vector3D &operator+=(const Vector3D &rhs);
    static Vector3D generate_vector(Random *rng);
    friend std::ostream &operator<<(std::ostream &oss, Vector3D vec);
};

std::ostream &operator<<(std::ostream &oss, Vector3D vec);

#endif // __VECTOR_3D__