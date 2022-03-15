#ifndef OBJECT
#define OBJECT

#include "Ray.h"
#include "Vector.h"

class Object
{
public:
    Vector albedo;
    bool isMirror, isTransp, isLight;
    double refraction_index;

    virtual bool intersect(const Ray &r, Vector &P, Vector &N, double &t, Vector &color) const = 0;
};

#endif