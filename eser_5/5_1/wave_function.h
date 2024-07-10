#ifndef __WAVE_FUNCTION__
#define __WAVE_FUNCTION__

#include <numeric>

#include "vector3D.h"

double ground_state(const Vector3D& pos);

double first_excited(const Vector3D& pos);

double second_excited(const Vector3D& pos);

#endif // __WAVE_FUNCTION__
