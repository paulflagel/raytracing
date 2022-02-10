#ifndef RAY
#define RAY

#include "../Vector/Vector.h"

class Ray
{
public:
    Ray(const Vector &origin, const Vector &direction);
    ~Ray();

    Vector C, u;
};

#endif