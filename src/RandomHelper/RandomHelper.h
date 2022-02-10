#ifndef RANDOMHELPER
#define RANDOMHELPER

#include <random>
#include "../Vector/Vector.h"

namespace randh
{
    Vector box_muller(float sigma = 1);
    Vector random_cos(Vector &N);
}

#endif