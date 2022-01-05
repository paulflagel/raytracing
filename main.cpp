#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <algorithm>
#include <iostream>

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
    explicit Sphere(const Vector &origin, double radius)
    {
        O = origin;
        R = radius;
    }

    bool intersect(const Ray &r, Vector &P, Vector &N) // Tester si la sphère est sur la
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
            double t;
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
};

int main()
{
    int W = 512;
    int H = 512;

    Vector C(0, 0, 55);            // Camera
    Sphere s(Vector(0, 0, 0), 10); // Sphere à l'origine, de rayon 10
    Vector L(-10, 20, 40);         // Source de lumière
    double fov = 60 * M_PI / 180;  // Field of view 60°
    double tanfov2 = tan(fov / 2);
    Vector rhoSphere(0.1, 0.4, 0.1); // Albedo de la sphère
    double I = 100000;               // Intensité de la lumière

    std::vector<unsigned char>
        image(W * H * 3, 0); // Crée un tableau 1D de W*H*3 éléments initialisés à 0 (l'image)

    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            Vector u(j - W / 2 + 0.5, (H - i) - (H / 2) + 0.5, -W / (2 * tanfov2)); // On ajoute 0.5 pour être au centre des pixels plutôt que dans les coins
            u.normalize();
            Ray r(C, u); // Rayon issu de C dans la direction u
            Vector P, n; // Point d'intersection et vecteur normal à la surface en P

            if (s.intersect(r, P, n))
            {
                Vector l = L - P;
                double lnorm2 = l.norm2();
                l.normalize();

                Vector col = rhoSphere * I * std::max(0., dot(l, n) / lnorm2 * 4 * M_PI);
                image[(i * W + j) * 3 + 0] = std::min(255., col[0]); // R
                image[(i * W + j) * 3 + 1] = std::min(255., col[1]); // G
                image[(i * W + j) * 3 + 2] = std::min(255., col[2]); // B
            }
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}