#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include "omp.h"

#define NB_THREADS 8

std::default_random_engine rng[NB_THREADS]; // On en crée 8 pour le calcul parallèle sur 8 threads
std::uniform_real_distribution<double> uniform(0.0, 1.0);

const float epsilon = 0.01;

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
Vector operator*(const Vector &a, const Vector &b) // Produit terme à terme
{
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
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
    explicit Sphere(const Vector &origin, double radius, const Vector &albedo, bool isMirror = false, bool isTransp = false, bool isLight = false, double refraction_index = 1.4)
    {
        this->O = origin;
        this->R = radius;
        this->albedo = albedo;
        this->isMirror = isMirror;
        this->isTransp = isTransp;
        this->isLight = isLight;
        this->refraction_index = refraction_index;
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
    bool isMirror, isTransp, isLight;
    double refraction_index;
};

Vector random_cos(Vector &N) // Retourne omega_i
{
    int thread_id = omp_get_thread_num();

    double r1 = uniform(rng[thread_id]);
    double r2 = uniform(rng[thread_id]);
    double s = sqrt(1 - r2); // Pour éviter de calculer 2 fois la racine

    double x = cos(2 * M_PI * r1) * s;
    double y = sin(2 * M_PI * r1) * s;
    double z = sqrt(r2);

    Vector T1;

    if ((N[2] < N[0]) && (N[2] < N[1])) // Nz plus petit
    {
        T1 = Vector(-N[1], N[0], 0.);
    }
    else if ((N[1] < N[0]) && (N[1] < N[2])) // Ny plus petit
    {
        T1 = Vector(N[2], 0., -N[0]);
    }
    else // Nx plus petit
    {
        T1 = Vector(0., N[2], -N[1]);
    }
    T1.normalize();
    Vector T2 = cross(N, T1);
    return x * T1 + y * T2 + z * N;
}

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
    Vector L; // Position de la lumière

    void add(Sphere &s) { liste_spheres.push_back(s); } // Ajoute la sphere à la fin de la liste

    bool intersect(const Ray &r, Vector &P, Vector &N, int &sphere_id, double &distance_min) // Check si une sphère est intersectée, si oui renvoie la plus proche
    {
        distance_min = 1e99;
        bool scene_intersection = false; // Est-ce que la scène possède au moins une intersection
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

    Vector getColor(Ray &r, int rebond) // Renvoie l'intensité du pixel
    {
        Vector color; // Intensité de la lumière

        if (rebond <= 0) // Condition de fin de récursion
        {
            return Vector(0., 0., 0.);
        }

        Vector P, n;   // Point d'intersection rayon/sphère et vecteur normal à la surface en P
        int id_sphere; // Id de la sphère intersectée
        double t;      // distance entre la caméra et l'intersection

        if (this->intersect(r, P, n, id_sphere, t))
        {
            // Sphère de lumière
            if (this->liste_spheres[id_sphere].isLight)
            {
                return this->liste_spheres[id_sphere].albedo * this->I / (4 * M_PI * this->liste_spheres[id_sphere].R * this->liste_spheres[id_sphere].R);
            }

            //  SURFACE MIROIR :
            if (this->liste_spheres[id_sphere].isMirror)
            {
                Vector uReflected = r.u - 2 * dot(r.u, n) * n;
                Ray rReflected(P + (epsilon * n), uReflected);
                return this->getColor(rReflected, rebond - 1);
            }

            // SURFACE TRANSPARENTE :
            else if (this->liste_spheres[id_sphere].isTransp)
            {
                double dot_incident_normal = dot(r.u, n);
                double n1, n2;
                if (dot_incident_normal < 0)
                {
                    n1 = 1.0;
                    n2 = this->liste_spheres[id_sphere].refraction_index;
                }
                else
                {
                    n1 = this->liste_spheres[id_sphere].refraction_index;
                    n2 = 1.0;
                    n = -1 * n;
                    dot_incident_normal = -dot_incident_normal;
                }

                double normal_sqrt_term = 1 - (n1 / n2) * (n1 / n2) * (1 - dot_incident_normal * dot_incident_normal);

                if (normal_sqrt_term < 0) // Réflexion totale
                {
                    Vector uReflected = r.u - 2 * dot(r.u, n) * n;
                    Ray rReflected(P + (epsilon * n), uReflected);
                    return this->getColor(rReflected, rebond - 1);
                }
                else
                {
                    Vector uRefractedN = -sqrt(normal_sqrt_term) * n;
                    Vector uRefractedT = (n1 / n2) * (r.u - dot_incident_normal * n);
                    Vector uRefracted = uRefractedN + uRefractedT;
                    Ray rRefracted(P - (epsilon * n), uRefracted);
                    return this->getColor(rRefracted, rebond - 1);
                }
            }

            // SURFACE DIFFUSE :
            else
            {
                Vector l = this->L - P;        // vecteur unitaire en P dirigé vers la source de lumière
                double lnorm2 = l.norm2();     // norme 2 de l
                double distlum = sqrt(lnorm2); // distance entre la lumière L et le point d'intersection de la boule P
                l.normalize();                 // on transforme l en vecteur unitaire

                // Grandeurs calculées pour l'ombre portée
                Vector Plum, nlum;
                double tlum;
                int idlum;
                int visibility = 1; // Pas d'ombre portée
                                    /*
                // ECLAIRAGE DIRECT

                // Calcul d'ombre portée
                if (this->intersect(Ray(P + (epsilon * n), l), Plum, nlum, idlum, tlum) && tlum < distlum)
                {
                    visibility = 0;
                }

                // // Soft shadows
                // Vector V = -1 * l;
                // Vector omega = random_cos(V);
                // Vector xprime = omega * this->liste_spheres[0].R + L;
                // Vector nprime = (xprime - L) / (xprime - L).norm();
                // double pdf = dot(l, omega) / (M_PI * this->liste_spheres[0].R * this->liste_spheres[0].R);
                // double jacobien = M_PI * std::max(0., dot(omega, n)) * std::max(0., dot(-1 * omega, nprime)) / ((xprime - P).norm2() * std::max(0., dot(V, nprime)));

                color = 1 * visibility * this->liste_spheres[id_sphere].albedo * this->I * std::max(0., dot(l, n)) / (lnorm2 * 4 * M_PI * M_PI); // clamp à 0 pour pas avoir une intensité négative
*/
                // ECLAIRAGE INDIRECT
                Vector omega_i = random_cos(n);
                Ray rRandom(P + (epsilon * n), omega_i);
                color = color + this->liste_spheres[id_sphere].albedo * this->getColor(rRandom, rebond - 1);
                return color;
            }
        }

        // PAS D'INTERSECTION : pixel noir
        return Vector(0., 0., 0.);
    }
};

// Génère deux nombres aléatoires suivant une loi normale centrée réduite à partir de deux nombres suivant une loi uniforme
const float sigma_bm = 0.3; // Std deviation of the normal distribution (Box Muller)
void box_muller(double &dx, double &dy)
{
    int thread_id = omp_get_thread_num();
    double r1 = uniform(rng[thread_id]);
    double r2 = uniform(rng[thread_id]);
    double r = sqrt(-2 * log(r2));
    double t = 2 * M_PI * r1;
    dx = cos(t) * r * sigma_bm;
    dy = sin(t) * r * sigma_bm;
};

int main()
{
    int W = 512;
    int H = 512;

    Vector C(0, 0, 55);           // Camera
    Vector L(-10, 20, 40);        // Source de lumière
    double fov = 60 * M_PI / 180; // Field of view 60°
    double tanfov2 = tan(fov / 2);
    double I(10000000000);        // Intensité de la lumière
    int nb_rays_monte_carlo = 64; // Nombre de rayons envoyés pour Monte Carlo
    int rebonds = 5;

    Scene scene(I, L);

    // Sphere 1 à l'origine, de rayon 10
    Sphere s1(Vector(0, 0, 0), 10, Vector(0.4, 0.1, 0.), false, false);
    // Sphere 2
    Sphere s2(Vector(-25, 20, -35), 10, Vector(0.4, 0.4, 0.4), false, false);
    // Sphère derrière la boule
    Sphere sBack(Vector(0, 0, -1000), 940, Vector(0., 0.5, 0.));
    // Sphère derrière la caméra
    Sphere sFront(Vector(0, 0, 1000), 940, Vector(0.5, 0., 0.5));
    // Plafond
    Sphere sUp(Vector(0, 1000, 0), 940, Vector(0.5, 0., 0.));
    // Sol
    Sphere sDown(Vector(0, -1000, 0), 990, Vector(0., 0., 0.5));
    // Mur droit
    Sphere sRight(Vector(1000, 0, 0), 940, Vector(0.5, 0.5, 0.5));
    // Mur gauche
    Sphere sLeft(Vector(-1000, 0, 0), 940, Vector(0.5, 0.5, 0.5));

    // Sphere de lumiere
    Sphere sLum(L, 10, Vector(1., 1., 1.), false, false, true);
    scene.add(sLum);

    scene.add(s1);
    scene.add(s2);
    scene.add(sFront);
    scene.add(sBack);
    scene.add(sUp);
    scene.add(sDown);
    scene.add(sRight);
    scene.add(sLeft);

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<unsigned char> image(W * H * 3, 0); // Crée un tableau 1D de W*H*3 éléments initialisés à 0 (l'image)
#pragma omp parallel for                            // Parallélisation du calcul, mettre un flag openmp à gcc. Ici on fait sur toutes les lignes de l'image
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            Vector intensity;
            for (int k = 0; k < nb_rays_monte_carlo; k++)
            {
                //Antialiasing : Plutôt que d'envoyer le rayon pile au centre des pixels, on l'envoie aléatoirement dans le pixel
                double dx = 0.0;
                double dy = 0.0;
                box_muller(dx, dy);

                Vector u(j - W / 2 + 0.5 + dx, (H - i) - (H / 2) + 0.5 + dy, -W / (2 * tanfov2)); // On ajoute 0.5 pour être au centre des pixels plutôt que dans les coins
                u.normalize();
                Ray r(C, u); // Rayon issu de C dans la direction u
                intensity = intensity + scene.getColor(r, rebonds);
            }
            intensity = intensity / nb_rays_monte_carlo;

            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(intensity[0], 1 / 2.2)); // R - La puissance 1/2.2 correspond à la correction gamma
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(intensity[1], 1 / 2.2)); // G
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(intensity[2], 1 / 2.2)); // B
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = end - start;
    auto diff_sec = std::chrono::duration_cast<std::chrono::seconds>(diff);
    std::cout << "Run time : " << diff_sec.count() << "s" << std::endl;

    return 0;
}
