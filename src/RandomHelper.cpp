#include <random>
#include <iostream>
#include <omp.h>
#include <cmath>
#include "Vector.h"
#include "RandomHelper.h"

RandomHelper::RandomHelper()
{
    std::random_device r;
    for (int i = 0, N = omp_get_max_threads(); i < N; ++i)
    {
        this->generators.push_back(std::default_random_engine(r()));
        this->uniforms.push_back(std::uniform_real_distribution<double>(0.0, 1.0));
    }
}

Vector RandomHelper::box_muller(float sigma)
{
    std::default_random_engine &engine = generators[omp_get_thread_num()];
    std::uniform_real_distribution<double> &uniform = uniforms[omp_get_thread_num()];
    // Génère deux nombres aléatoires suivant une loi normale centrée réduite à partir de deux nombres suivant une loi uniforme
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double r = sqrt(-2 * log(r2));
    double t = 2 * M_PI * r1;

    double dx = cos(t) * r * sigma;
    double dy = sin(t) * r * sigma;
    return Vector(dx, dy, 0.);
};

Vector RandomHelper::random_cos(Vector &N) // Retourne omega_i
{
    std::default_random_engine &engine = generators[omp_get_thread_num()];
    std::uniform_real_distribution<double> &uniform = uniforms[omp_get_thread_num()];

    double r1 = uniform(engine);
    double r2 = uniform(engine);
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

double RandomHelper::random_double()
{
    std::default_random_engine &engine = generators[omp_get_thread_num()];
    std::uniform_real_distribution<double> &uniform = uniforms[omp_get_thread_num()];

    return uniform(engine);
}
