#include <stdio.h>
#include <stdlib.h>
#include "rng.h"
#include "params.h"
#include "traversal.h"
#include "type.h"
#include "kernels.h"
#include "timer.h"
#include "verify.h"
#include "initialise.h"
#include "finalise.h"
#include "util.h"

double *minifmm(double *weights, double *points, long num_points, int num_terms, int ncri, float theta)
{
    t_fmm_params params;
    init_with_known_data(&params, weights, points, num_points, num_terms, ncri, theta);
    printf("starting computation\n");
    printf(SEPERATOR);
    perform_fmm(&params);
    printf(SEPERATOR);
    return params.pot;
}

int minifmm_no_args()
{
    long num_points = 100000;
    TYPE* points = (TYPE*)malloc(sizeof(TYPE)*num_points*3);
    TYPE* weights = (TYPE*)malloc(sizeof(TYPE)*num_points);
    for (size_t i = 0; i < num_points*3; ++i)
    {
       points[i] = rand_range(-1.0, 1.0);
        if(i%3==0)
            weights[i] = rand_range(-1.0, 1.0);
    }
    minifmm(weights, points, num_points, 7, 200, 0.5);
    return 0;
}
