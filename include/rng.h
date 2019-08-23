#pragma once

#include <stdint.h>
#include <stdlib.h>

inline void seed_rng(uint64_t seed)
{
    srand(seed);
}

inline double gen_rand()
{
    return (double)rand()/(double)RAND_MAX;
}

inline double rand_range(double min, double max)
{
    double r = (double)gen_rand();
    return (max-min)*r + min;
}

