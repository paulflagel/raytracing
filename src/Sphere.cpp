#include "Sphere.h"
#include <cmath>

Sphere::Sphere(const Vector &origin, double radius, const Vector &albedo, bool isMirror, bool isTransp, bool isLight, double refraction_index)
{
    this->O = origin;
    this->R = radius;
    this->albedo = albedo;
    this->isMirror = isMirror;
    this->isTransp = isTransp;
    this->refraction_index = refraction_index;
    this->isLight = isLight;
}

Sphere::~Sphere() {}

bool Sphere::intersect(const Ray &r, Vector &P, Vector &N, double &t, Vector &color) const // Tester si la sphère est sur la trajectoire
{
    color = this->albedo;
    // On résoud l'équation du second degré a*T^2 + b*t + c = 0
    double a = 1;
    double b = 2 * dot(r.u, r.C - O);
    double c = (r.C - O).norm2() - R * R;

    double delta = b * b - 4 * a * c;

    if (delta >= 0)
    {
        double t1 = (-b - sqrt(delta)) / (2 * a);
        double t2 = (-b + sqrt(delta)) / (2 * a);
        if (t2 < 0) // Les deux solutions sont négatives
            return false;
        if (t1 < 0)
            t = t2;
        else
            t = t1;

        P = r.C + t * r.u;
        N = P - O;
        N.normalize();
        return true;
    }
    return false;
}