#include "Scene.h"
#include <cmath>
#include <random>
#include "omp.h"
#include <iostream>

#define EPSILON 0.01

Scene::Scene(RandomHelper &r)
{
    this->randh = r;
}

Scene::~Scene() {}

void Scene::add(Object *obj) { objects_list.push_back(obj); } // Ajoute la sphere à la fin de la liste

bool Scene::intersect(const Ray &r, Vector &P, Vector &N, int &sphere_id, double &distance_min, Vector &albedo) // Check si une sphère est intersectée, si oui renvoie la plus proche
{
    distance_min = 1E99;
    bool scene_intersection = false; // Est-ce que la scène possède au moins une intersection
    for (int k = 0; k < objects_list.size(); k++)
    {
        Vector Pobject, Nobject;
        Vector albedoObject;
        double local_distance;

        bool local_intersection = objects_list[k]->intersect(r, Pobject, Nobject, local_distance, albedoObject);
        if (local_intersection)
        {
            scene_intersection = true;
            if (local_distance < distance_min)
            {
                distance_min = local_distance;
                P = Pobject;
                N = Nobject;
                sphere_id = k;
                albedo = albedoObject;
            }
        }
    }
    return scene_intersection;
}

Vector Scene::getColor(Ray &r, int rebond, bool showLight) // Renvoie l'intensité du pixel
{
    if (rebond <= 0) // Condition de fin de récursion
    {
        return Vector(0., 0., 0.);
    }

    Vector P, n, albedo; // Point d'intersection rayon/sphère et vecteur normal à la surface en P
    int id_sphere;       // Id de la sphère intersectée
    double t;            // distance entre la caméra et l'intersection
    bool intersection_bool = this->intersect(r, P, n, id_sphere, t, albedo);

    Vector color; // Intensité de la lumière

    if (intersection_bool)
    {
        Object *oCurrent = objects_list[id_sphere]; // Sphere intersectée

        // Si on arrive directement sur la lumière :
        if ((oCurrent->isLight) && showLight)
        {
            return I * Vector(1., 1., 1.);
        }

        //  SURFACE MIROIR :
        if (oCurrent->isMirror)
        {
            Vector uReflected = r.u - 2 * dot(r.u, n) * n;
            Ray rReflected(P + (EPSILON * n), uReflected);
            return this->getColor(rReflected, rebond - 1, showLight);
        }

        // SURFACE TRANSPARENTE :
        else if (oCurrent->isTransp)
        {
            double dot_incident_normal = dot(r.u, n);
            double n1, n2;
            if (dot_incident_normal < 0)
            {
                n1 = 1.0;
                n2 = oCurrent->refraction_index;
            }
            else
            {
                n1 = oCurrent->refraction_index;
                n2 = 1.0;
                n = -n;
                dot_incident_normal = -dot_incident_normal;
            }

            double normal_sqrt_term = 1 - (n1 / n2) * (n1 / n2) * (1 - dot_incident_normal * dot_incident_normal);

            double k0 = ((n1 - n2) / (n1 + n2)) * ((n1 - n2) / (n1 + n2));
            double R = k0 + (1 - k0) * std::pow(1 + dot_incident_normal, 5);
            bool is_fresnel = fresnel && (randh.random_double() < R);

            if (normal_sqrt_term < 0 || is_fresnel) // Réflexion totale
            {
                Vector uReflected = r.u - 2 * dot(r.u, n) * n;
                Ray rReflected(P + (EPSILON * n), uReflected);
                return this->getColor(rReflected, rebond - 1, showLight);
            }
            else
            {
                Vector uRefractedN = -sqrt(normal_sqrt_term) * n;
                Vector uRefractedT = (n1 / n2) * (r.u - dot_incident_normal * n);
                Vector uRefracted = uRefractedN + uRefractedT;
                Ray rRefracted(P - (EPSILON * n), uRefracted);
                return this->getColor(rRefracted, rebond - 1);
            }
        }

        // SURFACE DIFFUSE :
        else
        {
            // ECLAIRAGE DIRECT

            // Grandeurs calculées pour l'ombre portée
            Vector Plum, nlum, albedoLum;
            double tlum;
            int idlum;

            // Ombres douces
            if (this->soft_shadows)
            {
                Vector v = P - this->Light->O;
                v.normalize();
                Vector omega_random = randh.random_cos(v);                        // On met -> parce que c'est un pointeur
                omega_random.normalize();                                         // Direction aléatoire sur l'hémisphère de sLum
                Vector x_random = omega_random * this->Light->R + this->Light->O; // Point à la surface de sLum issu de la direction aléatoire
                Vector omega_i = (x_random - P);
                omega_i.normalize();
                double distlum2 = (x_random - P).norm2();

                // Calcul d'ombre portée : on regarde s'il y a un obstacle entre le point et la lumière
                if (this->intersect(Ray(P + (EPSILON * n), omega_i), Plum, nlum, idlum, tlum, albedoLum) && tlum * tlum < distlum2 * 0.99)
                {
                    // Ombre portée
                    color = Vector(0., 0., 0.);
                }
                else
                {
                    double pdf = std::max(0., dot(v, omega_random)) / M_PI / (this->Light->R * this->Light->R);
                    Vector BRDF = albedo / M_PI;
                    double jacobien = std::max(0., dot(-omega_i, omega_random)) / distlum2;
                    color = this->I / (4 * M_PI * this->Light->R * this->Light->R * M_PI) * std::max(0., dot(n, omega_i)) * BRDF * jacobien / pdf;
                }
            }
            // Pas d'ombres douces
            else
            {
                Vector l = this->Light->O - P;
                double distlum2 = l.norm2();
                l.normalize();
                if (this->intersect(Ray(P + (EPSILON * n), l), Plum, nlum, idlum, tlum, albedoLum) && (tlum + this->Light->R) * (tlum + this->Light->R) < distlum2 * 0.99)
                {
                    color = Vector(0., 0., 0.);
                }
                else
                {
                    Vector BRDF = albedo / M_PI;
                    color = this->I / (4 * M_PI * distlum2) * std::max(0., dot(l, n)) * BRDF;
                }
            }

            // ECLAIRAGE INDIRECT
            if (this->indirect_light)
            {
                Vector random_u = randh.random_cos(n);
                Ray rRandom(P + (EPSILON * n), random_u);
                color = color + albedo * this->getColor(rRandom, rebond - 1);
            }
            return color;
        }
    }

    // PAS D'INTERSECTION : pixel noir
    return Vector(0., 0., 0.);
};
