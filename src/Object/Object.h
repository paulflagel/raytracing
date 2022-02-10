#ifndef OBJECT
#define OBJECT

#include "../Ray/Ray.h"
#include "../Vector/Vector.h"

class Object
{
public:
    Vector O, albedo;
    bool isMirror, isTransp;
    double refraction_index;

    virtual bool intersect(const Ray &r, Vector &P, Vector &N, double &t) const = 0;
};

#endif