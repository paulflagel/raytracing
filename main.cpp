#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <algorithm>
#include <iostream>
#include <cmath>

// & : au lieu de passer en paramètre le vecteur soit 3 fois 64 bits, on passe une adresse qui est 1 fois 64 bits

class Vector
{
public:
    explicit Vector(double x = 0, double y = 0, double z = 0)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const
    {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const
    {
        return sqrt(norm2());
    }
    void normalize()
    {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; }; //renvoie un double constant (qu'on ne peut pas modifier) genre a[0]
    double &operator[](int i) { return data[i]; };      //renvoie la référence vers un double pas constant (pas de const pour pouvoir le modifier) genre a[0] = 2
    double data[3];                                     //tableau statique de 3 éléments (en précision double), qu'on remplit avec les coordonnées dans le costructeur
};

// Opérateurs définis pour les vecteurs
Vector operator+(const Vector &a, const Vector &b)
{
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b)
{
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector &b)
{
    return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector &a, const double b)
{
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator/(const Vector &a, const double b)
{
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector &a, const Vector &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector &a, const Vector &b)
{
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray
{
public:
    explicit Ray(const Vector &origin, const Vector &direction)
    {
        C = origin;
        u = direction;
    }
    Vector C, u;
};

class Sphere
{
public:
    explicit Sphere(const Vector &origin, double radius, const Vector &albedo, bool mirror = false)
    {
        this->O = origin;
        this->R = radius;
        this->albedo = albedo;
        this->mirror = mirror;
    }

    bool intersect(const Ray &r, Vector &P, Vector &N, double &t) // Tester si la sphère est sur la trajectoire
    {
        //On résoud l'équation du second degré a*T^2 + b*t + c = 0
        double a = 1;
        double b = 2 * dot(r.u, r.C - O);
        double c = (r.C - O).norm2() - R * R;

        double delta = b * b - 4 * a * c;

        if (delta >= 0)
        {
            double t1 = (-b - sqrt(delta)) / (2 * a);
            double t2 = (-b + sqrt(delta)) / (2 * a);
            if (t2 < 0) // Les deux solutions sont négatives
            {
                return false;
            }
            if (t1 < 0)
            {
                t = t2;
            }
            else
            {
                t = t1;
            }

            P = r.C + t * r.u;
            N = P - O;
            N.normalize();
            return true;
        }
        return false;
    }

    Vector O;
    double R;
    Vector albedo;
    bool mirror;
};

class Scene
{
public:
    Scene(double I, Vector &L)
    {
        this->I = I;
        this->L = L;
    }

    std::vector<Sphere> liste_spheres;
    double I; // Intensité de la lumière
    Vector L; //Position de la lumière

    void add(Sphere &s) { liste_spheres.push_back(s); } //Ajoute la sphere à la fin de la liste

    bool intersect(const Ray &r, Vector &P, Vector &N, int &sphere_id, double &distance_min)
    {
        distance_min = 1e99;
        bool scene_intersection = false; // est-ce que la scène possède au moins une intersection
        for (int k = 0; k < liste_spheres.size(); k++)
        {
            Vector Psphere, Nsphere;
            double local_distance;

            bool local_intersection = liste_spheres[k].intersect(r, Psphere, Nsphere, local_distance);
            if (local_intersection)
            {
                scene_intersection = true;
                if (local_distance < distance_min)
                {
                    distance_min = local_distance;
                    P = Psphere;
                    N = Nsphere;
                    sphere_id = k;
                }
            }
        }
        return scene_intersection;
    }

    Vector getColor(Ray &r, int rebond)
    {
        if (rebond <= 0)
        {
            return Vector(0, 0, 0);
        }
        Vector P, n; // Point d'intersection rayon/sphère et vecteur normal à la surface en P
        int id_sphere;
        double t; // distance entre la caméra et l'intersection

        if (this->intersect(r, P, n, id_sphere, t))
        {
            if (!this->liste_spheres[id_sphere].mirror)
            {
                Vector l = this->L - P;        // vecteur unitaire en P dirigé vers la source de lumière
                double lnorm2 = l.norm2();     // norme 2 de l
                double distlum = sqrt(lnorm2); // distance entre la lumière L et le point d'intersection de la boule P
                l.normalize();                 // l est un vecteur unitaire
                Vector rho = this->liste_spheres[id_sphere].albedo;

                Vector Plum, nlum;
                double tlum;
                int idlum;

                if (this->intersect(Ray(P + (0.01 * n), l), Plum, nlum, idlum, tlum) && tlum < distlum) // Calcul d'ombre portée
                {
                    return Vector(0, 0, 0);
                }
                else
                {
                    return rho * this->I * std::max(0., dot(l, n)) / (lnorm2 * 4 * M_PI); //on clamp à 0 pour pas avoir une intensité négative
                }
            }
            else
            {
                Vector uReflected = r.u - 2 * dot(r.u, n) * n;
                Ray rReflected(P + (0.01 * n), uReflected);
                return this->getColor(rReflected, rebond - 1);
            }
        }
        return Vector(0, 0, 0);
    }
};

int main()
{
    int W = 512;
    int H = 512;

    Vector C(0, 0, 55);           // Camera
    Vector L(-10, 20, 40);        // Source de lumière
    double fov = 60 * M_PI / 180; // Field of view 60°
    double tanfov2 = tan(fov / 2);
    double I = 5000000000; // Intensité de la lumière
    bool mirror = true;

    Scene scene(I, L);

    // Sphere 1 à l'origine, de rayon 10
    Vector albedo1(0.4, 0.1, 0.);
    Sphere s1(Vector(0, 0, 0), 10, albedo1, mirror);
    // Sphere 2
    Vector albedo2(0.3, 0., 0.3);
    Sphere s2(Vector(-20, 15, -20), 5, albedo2, mirror);
    // Sphere 3
    Vector albedo3(0.3, 0.3, 0.1);
    Sphere s3(Vector(20, -2, 5), 3, albedo3, mirror);
    // Sphère derrière la caméra (magenta)
    Vector albedoBack(0.4, 0.3, 0.5);
    Sphere sBack(Vector(0, 0, -3040), 3000, albedoBack);
    // Sphère en haut (rouge)
    Vector albedoUp(0.4, 0.2, 0.6);
    Sphere sUp(Vector(0, 3040, 0), 3000, albedoUp);
    // Sphère en bas (bleue)
    Vector albedoDown(0.4, 0.2, 0.6);
    Sphere sDown(Vector(0, -3015, 0), 3000, albedoDown);

    scene.add(s1);
    scene.add(s2);
    scene.add(s3);
    // scene.add(sFront);
    scene.add(sBack);
    scene.add(sUp);
    scene.add(sDown);

    std::vector<unsigned char> image(W * H * 3, 0); // Crée un tableau 1D de W*H*3 éléments initialisés à 0 (l'image)

    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            Vector u(j - W / 2 + 0.5, (H - i) - (H / 2) + 0.5, -W / (2 * tanfov2)); // On ajoute 0.5 pour être au centre des pixels plutôt que dans les coins
            u.normalize();
            Ray r(C, u); // Rayon issu de C dans la direction u
            int rebonds = 5;
            Vector intensity = scene.getColor(r, rebonds);

            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(intensity[0], 1 / 2.2)); // R - La puissance 1/2.2 correspond à la correction gamma
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(intensity[1], 1 / 2.2)); // G
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(intensity[2], 1 / 2.2)); // B
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}
