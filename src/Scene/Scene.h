#ifndef SCENE
#define SCENE

#include "../Vector/Vector.h"
#include "../Ray/Ray.h"
#include "../Sphere/Sphere.h"
#include "../RandomHelper/RandomHelper.h"
#include "../Object/Object.h"
#include <vector>

class Scene
{
public:
    Scene(double I);
    ~Scene();

    std::vector<Object *> objects_list;
    double I; // Intensité de la lumière
    Sphere *Light;

    void add(Object *obj);                                                                    // Ajoute la sphere à la fin de la liste
    bool intersect(const Ray &r, Vector &P, Vector &N, int &sphere_id, double &distance_min); // Check si une sphère est intersectée, si oui renvoie la plus proche
    Vector getColor(Ray &r, int rebond);                                                      // Renvoie l'intensité du pixel

    bool indirect_light = true;
    bool soft_shadows = true;
    int max_rebonds;
};

#endif