#ifndef SPHERE
#define SPHERE

#include "../Vector/Vector.h"
#include "../Ray/Ray.h"
#include "../Object/Object.h"

class Sphere : public Object
{
public:
    Sphere(const Vector &origin, double radius, const Vector &albedo, bool isMirror = false, bool isTransp = false, bool isLight = false, double refraction_index = 1.4);
    ~Sphere();

    virtual bool intersect(const Ray &r, Vector &P, Vector &N, double &t) const;

    Vector O;
    double R;

    // Already defined in Object !

    // Vector O;
    // Vector albedo;
    // bool isMirror, isTransp, isLight;
    // double refraction_index;
};

#endif