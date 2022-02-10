#include <random>
#include <iostream>
#include <omp.h>
#include <cmath>
#include "../Vector/Vector.h"
#include "RandomHelper.h"

#define NUM_THREADS 8

std::default_random_engine engine[NUM_THREADS];
std::uniform_real_distribution<double> uniform(0.0, 1.0);

Vector randh::box_muller(float sigma)
{
    // Génère deux nombres aléatoires suivant une loi normale centrée réduite à partir de deux nombres suivant une loi uniforme
    int thread_id = omp_get_thread_num();
    double r1 = uniform(engine[thread_id]);
    double r2 = uniform(engine[thread_id]);
    double r = sqrt(-2 * log(r2));
    double t = 2 * M_PI * r1;

    double dx = cos(t) * r * sigma;
    double dy = sin(t) * r * sigma;
    return Vector(dx, dy, 0.);
};

Vector randh::random_cos(Vector &N) // Retourne omega_i
{
    int thread_id = omp_get_thread_num();

    double r1 = uniform(engine[thread_id]);
    double r2 = uniform(engine[thread_id]);
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
