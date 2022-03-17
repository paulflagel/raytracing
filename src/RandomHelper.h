#ifndef RANDOMHELPER
#define RANDOMHELPER

#include <random>
#include "Vector.h"

class RandomHelper
{
public:
    std::vector<std::default_random_engine> generators;
    std::vector<std::uniform_real_distribution<double>> uniforms;

    RandomHelper();
    void init_helper();
    Vector box_muller(float sigma = 1);
    Vector random_cos(Vector &N);
    double random_double();
};

#endif