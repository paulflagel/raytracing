#define _CRT_SECURE_NO_WARNINGS 1
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <chrono>
#include <string>
#include "omp.h"

#include "Vector.h"
#include "Ray.h"
#include "Sphere.h"
#include "Scene.h"
#include "RandomHelper.h"
#include "Object.h"
#include "Triangle.h"

#include "progressbar.hpp"

#define INDIRECT_LIGHT true
#define SOFT_SHADOWS true
#define NUM_RAYS_MC 128
#define ANTIALIASING true
#define DEPTH_OF_FIELD false
#define DDOF 55
#define USE_BVH true
#define NORMALS_INTERPOLATION true
#define FRESNEL false
#define CAMERA_MOVE false

void print(Vector a)
{
    std::cout << a[0] << " " << a[1] << " " << a[2] << std::endl;
}

int main(int argc, char *argv[])
{
    int W = 512;
    int H = 512;

    Vector C(0, 0, 55);           // Camera
    Vector L(-10, 20, 40);        // Light source
    double fov = 80 * M_PI / 180; // Field of view 60°
    double tanfov2 = tan(fov / 2);
    double I = 1E10;         // Light intensity (Watts)
    int nb_rays_monte_carlo; // Number of rays for Monte Carlo (error decreases in 1/sqrt(N))
    if (INDIRECT_LIGHT || SOFT_SHADOWS || ANTIALIASING)
        nb_rays_monte_carlo = NUM_RAYS_MC;
    else
        nb_rays_monte_carlo = 1;
    int rebonds = 5;
    double ddof = DDOF;

    RandomHelper randh; // Class instance containing thread-safe random generators

    // ========== Scene init ==========
    Scene scene(randh);
    scene.I = I;
    scene.indirect_light = INDIRECT_LIGHT;
    scene.soft_shadows = SOFT_SHADOWS;
    scene.fresnel = FRESNEL;

    // ========== Load mesh ==========

    TriangleMesh tri(Vector(0.4, 0.1, 0.1), false, false);
    TriangleMesh tri2(Vector(0.4, 0.1, 0.1), false, false);
    TriangleMesh tri3(Vector(0.4, 0.1, 0.1), false, false);

    // Final scene
    tri.readOBJ("mesh/saturn/13906_Saturn_v1_l3.obj");
    tri.load_texture("mesh/saturn/Saturn_diff.jpg");
    tri.rotate(0, 90);
    tri.rotate(2, -45);
    tri.rotate(0, -20);
    tri.scale(0.04, Vector(15, 10, -10));

    tri2.readOBJ("mesh/satellite/Satellite.obj");
    tri2.load_texture("mesh/satellite/Textures/satellite_Antenna_BaseColor.jpg");
    tri2.load_texture("mesh/satellite/Textures/satellite_Satélite_BaseColor.jpg");
    tri2.load_texture("mesh/satellite/Textures/satellite_Pinos_BaseColor.jpg");
    tri2.load_texture("mesh/satellite/Textures/satellite_Placas_BaseColor.jpg");
    tri2.load_texture("mesh/satellite/Textures/satellite_Couro_BaseColor.jpg");
    tri2.rotate(0, -25);
    tri2.rotate(1, 45);
    tri2.rotate(0, -20);
    tri2.scale(3, Vector(-10, 20, 0));

    tri3.readOBJ("mesh/rocket/10475_Rocket_Ship_v1_L3.obj");
    tri3.load_texture("mesh/rocket/10475_Rocket_Ship_v1_Diffuse.jpg");
    tri3.rotate(1, -90);
    tri3.rotate(0, 180);
    tri3.rotate(2, -45);
    tri3.rotate(1, -30);
    tri3.scale(0.04, Vector(-20, -5, 20));

    tri.normals_interpolation = NORMALS_INTERPOLATION;
    tri2.normals_interpolation = NORMALS_INTERPOLATION;
    tri3.normals_interpolation = NORMALS_INTERPOLATION;
    if (USE_BVH)
    {
        tri.use_bvh = USE_BVH;
        tri2.use_bvh = USE_BVH;
        tri3.use_bvh = USE_BVH;
        std::cout << "Init BVH..." << std::endl;
        tri.init_BVH();
        tri2.init_BVH();
        tri3.init_BVH();
        std::cout << "Init BVH done" << std::endl;
    }
    else
    {
        tri.use_bvh = USE_BVH;
        tri.get_bbox();
    }

    // ========== Spheres init ==========
    // Sphere s1(Vector(0, 0, 0), 10, Vector(0.4, 0.1, 0.), false, true);
    // Sphere s2(Vector(-20, 0, -5), 10, Vector(0.4, 0.4, 0.4), false, false);
    // Sphere s3(Vector(20, 0, -5), 10, Vector(0.1, 0., 0.5), false, false);

    // Sphere sBack(Vector(0, 0, -1000), 940, Vector(0., 0.5, 0.));   // Sphère derrière la boule
    // Sphere sFront(Vector(0, 0, 1000), 940, Vector(0.5, 0., 0.5));  // Sphère derrière la caméra
    // Sphere sUp(Vector(0, 1000, 0), 940, Vector(0.5, 0., 0.));      // Plafond
    // Sphere sDown(Vector(0, -1000, 0), 990, Vector(0., 0., 0.5));   // Sol
    // Sphere sRight(Vector(1000, 0, 0), 940, Vector(0.5, 0.5, 0.5)); // Mur droit
    // Sphere sLeft(Vector(-1000, 0, 0), 940, Vector(0.5, 0.5, 0.5)); // Mur gauche

    Vector albedo_spheres = Vector(0.5, 0.5, 0.6);

    Sphere sBack(Vector(0, 0, -1000), 940, albedo_spheres); // Sphère derrière la boule
    Sphere sFront(Vector(0, 0, 1000), 940, albedo_spheres); // Sphère derrière la caméra
    Sphere sUp(Vector(0, 1000, 0), 940, albedo_spheres);    // Plafond
    Sphere sDown(Vector(0, -1000, 0), 990, albedo_spheres); // Sol
    Sphere sRight(Vector(1000, 0, 0), 940, albedo_spheres); // Mur droit
    Sphere sLeft(Vector(-1000, 0, 0), 940, albedo_spheres); // Mur gauche

    Sphere sLum(L, 5, Vector(1., 1., 1.), false, false, true); // Sphere de lumiere
    scene.add(&sLum);
    scene.Light = &sLum;

    // scene.add(&s1);
    // scene.add(&s2);
    // scene.add(&s3);

    scene.add(&tri);
    scene.add(&tri2);
    scene.add(&tri3);

    scene.add(&sFront);
    scene.add(&sBack);
    scene.add(&sUp);
    scene.add(&sDown);
    // scene.add(&sRight);
    // scene.add(&sLeft);

    const int dTheta = 3; // Degrees for shifting camera
    int nb_images = 1;
    if (CAMERA_MOVE)
        int nb_images = 360 / dTheta;

    progressbar bar(H * nb_images); // Progress bar from https://github.com/gipert/progressbar

    auto start = std::chrono::high_resolution_clock::now(); // Clock counting render time

    for (int l = 0; l < nb_images; l++)
    {
        double angle = dTheta * l * M_PI / 180; // 0° to 360°
        double cosinus = cos(angle);
        double sinus = sin(angle);
        C[0] = 55 * sinus;
        C[1] = 10 * sinus;
        C[2] = 55 * cosinus;

        std::vector<unsigned char> image(W * H * 3, 0); // Crée un tableau 1D de W*H*3 éléments initialisés à 0 (l'image)
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < H; i++)
        {
#pragma omp critical
            bar.update();
            for (int j = 0; j < W; j++)
            {
                Vector intensity;
                for (int k = 0; k < nb_rays_monte_carlo; k++)
                {
                    Vector dPixel;
                    if (ANTIALIASING)
                        dPixel = randh.box_muller(0.3);
                    Vector u(j - W / 2 + 0.5 + dPixel[0], (H - i) - (H / 2) + 0.5 + dPixel[1], -W / (2 * tanfov2));
                    u.normalize();

                    // ========== Rotation ==========
                    double u0 = u[0];
                    double u2 = u[2];
                    u[0] = u0 * cosinus + u2 * sinus;
                    u[2] = -u0 * sinus + u2 * cosinus;
                    u.normalize();
                    Ray r(C, u);

                    if (DEPTH_OF_FIELD)
                    {
                        Vector dAperture = randh.box_muller(1.);
                        Vector Cprime = C + dAperture;
                        Vector uprime = C + ddof / abs(u[2]) * u - Cprime;
                        uprime.normalize();
                        r = Ray(Cprime, uprime); // Rayon issu de C dans la direction u
                    }

                    intensity = intensity + scene.getColor(r, rebonds, true);
                }
                intensity = intensity / nb_rays_monte_carlo;

                image[(i * W + j) * 3 + 0] = std::min(255., std::pow(intensity[0], 1 / 2.2)); // R - Pow 1/2.2 is the gamma correction
                image[(i * W + j) * 3 + 1] = std::min(255., std::pow(intensity[1], 1 / 2.2)); // G
                image[(i * W + j) * 3 + 2] = std::min(255., std::pow(intensity[2], 1 / 2.2)); // B
            }
        }
        if (CAMERA_MOVE)
        {
            std::string s = "output/img" + std::to_string(l) + ".png";
            const char *c = s.c_str();
            stbi_write_png(c, W, H, 3, &image[0], 0);
        }
        else
        {
            stbi_write_png("image.png", W, H, 3, &image[0], 0);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto diff = end - start;
    auto diff_sec = std::chrono::duration_cast<std::chrono::milliseconds>(diff);
    std::cout << std::endl
              << "Run time : " << diff_sec.count() << "ms (" << diff_sec.count() / 1000. << "s)" << std::endl;

    return 0;
}
