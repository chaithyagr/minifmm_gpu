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
#include "parse_args.h"
#include "util.h"
const static size_t DEFAULT_NUM_POINTS      = 10000;
const static size_t DEFAULT_NCRIT           = 200;
const static int DEFAULT_NUM_TERMS          = 10;
const static TYPE DEFAULT_THETA             = 0.5;
const static size_t DEFAULT_NUM_SAMPLES     = 100;

void init_defaults(t_fmm_params* params)
{
    params->num_points     = DEFAULT_NUM_POINTS;
    params->ncrit          = DEFAULT_NCRIT;
    params->num_terms      = DEFAULT_NUM_TERMS;
    params->theta          = DEFAULT_THETA;
    params->num_samples    = DEFAULT_NUM_SAMPLES;
}

double * minifmm()
{
    struct t_fmm_params params = t_fmm_params();
    init_defaults(&params);
    params.num_multipoles = params.num_terms*params.num_terms;
    params.num_spharm_terms = params.num_terms*params.num_terms;
    init_data(&params);
    precompute(&params);
    build_tree(&params);
    print_args(&params);
    printf("starting computation\n");
    printf(SEPERATOR);
    perform_fmm(&params);
    printf(SEPERATOR);
    double *ptrOut =(double *) malloc(sizeof(double)*params.num_points);
    for(int l = 0 ; l < params.num_points ; l++)
    {
        ptrOut[l] = params.pot[l];
    }
    finalise(&params);
    return ptrOut;
}

int  main()
{
    double *ptrOut = minifmm();
    ptrOut = minifmm();
    return 0;
}
