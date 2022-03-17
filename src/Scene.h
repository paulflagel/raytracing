#ifndef SCENE
#define SCENE

#include "Vector.h"
#include "Ray.h"
#include "Sphere.h"
#include "RandomHelper.h"
#include "Object.h"
#include <vector>

class Scene
{
public:
    Scene(RandomHelper &r);
    ~Scene();

    std::vector<Object *> objects_list;
    double I; // Intensité de la lumière
    Sphere *Light;

    void add(Object *obj);                                                                                    // Ajoute la sphere à la fin de la liste
    bool intersect(const Ray &r, Vector &P, Vector &N, int &sphere_id, double &distance_min, Vector &albedo); // Check si une sphère est intersectée, si oui renvoie la plus proche
    Vector getColor(Ray &r, int rebond, bool showLight = false);                                              // Renvoie l'intensité du pixel

    bool indirect_light = true;
    bool soft_shadows = true;
    bool fresnel = true;

    RandomHelper randh;
};

#endif