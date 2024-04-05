#ifndef __VECTOR_3D__
#define __VECTOR_3D__

#include <iostream>
#include <numeric>

#include "random.h"

struct Vector3D
{
    double x{0.}, y{0.}, z{0.};

    double distance(const Vector3D &ref = {0.0, 0.0, 0.0}) const { return std::hypot((x - ref.x), (y - ref.y), (z - ref.z)); }
    double distance_squared(const Vector3D &ref = {0.0, 0.0, 0.0}) const { return distance(ref) * distance(ref); }

    /// @brief Dot Product
    /// @param rhs Other Vector
    /// @return scalar product
    double operator*(const Vector3D &rhs) { return x * rhs.x + y * rhs.y + z * rhs.z; }

    /// @brief Scalar Product
    /// @param rhs A scalar
    /// @return The original vector scaled by rhs
    Vector3D operator*(const double &rhs) { return Vector3D{x * rhs, y * rhs, z * rhs}; }
    
    /// @brief Sum of Vectors
    /// @param rhs Another Vector
    /// @return 
    Vector3D operator+(const Vector3D &rhs) const { return Vector3D{x + rhs.x, y + rhs.y, z + rhs.z}; }
    Vector3D &operator+=(const Vector3D &rhs);

    /// @brief Generate a vector with radius one
    /// @param rng an exteral random generator
    /// @return 
    static Vector3D generate_vector(Random *rng);
    
    
    /// @brief Generate a vector using a uniform distribution between delta and - delta.
    /// @param rng An external Random number generator 
    /// @param delta Max number generable
    /// @return A new vector with x, y, z between delta and -delta
    static Vector3D generate_unif(Random *rng, const double delta);

    /// @brief Generate a vector using a normal distribution with mean zero and variance delta
    /// @param rng An external random number generator
    /// @param delta Variance to be used in the normal distribution
    /// @return A new vector with x, y and  sampled from a normal distribution with mean 0 and variance delta
    static Vector3D generate_gauss(Random *rng, const double delta);

    friend std::ostream &operator<<(std::ostream &oss, const Vector3D& vec);
};

// std::ostream &operator<<(std::ostream &oss, Vector3D& vec);

#endif // __VECTOR_3D__