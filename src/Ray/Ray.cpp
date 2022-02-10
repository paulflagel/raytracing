#include "../Vector/Vector.h"
#include "Ray.h"

Ray::Ray(const Vector &origin, const Vector &direction)
{
    C = origin;
    u = direction;
}

Ray::~Ray() {}
