#ifndef SPHERE
#define SPHERE

#include "../Vector/Vector.h"
#include "../Ray/Ray.h"

class Sphere
{
public:
    Sphere(const Vector &origin, double radius, const Vector &albedo, bool isMirror = false, bool isTransp = false, bool isLight = false, double refraction_index = 1.4);
    ~Sphere();
    bool intersect(const Ray &r, Vector &P, Vector &N, double &t);

    Vector O;
    double R;
    Vector albedo;
    bool isMirror, isTransp, isLight;
    double refraction_index;
};

#endif