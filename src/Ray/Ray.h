#ifndef RAY
#define RAY

class Ray
{
public:
    Ray(const Vector &origin, const Vector &direction);
    ~Ray();

    Vector C, u;
};

#endif