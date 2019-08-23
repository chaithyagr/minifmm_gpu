#pragma once

#include <vector>
#include "type.h"

struct t_node
{
    TYPE center[3];
    size_t num_children;
    size_t child[8];
    size_t num_points;
    TYPE rad;
    size_t arr_idx;
    int level;
    int num_desc_nodes;

    size_t point_idx;
    size_t mult_idx;
    size_t offset;
} ;


struct t_fmm_params
{
    size_t root;
    int *reorder, *reorder_ordered;
    bool dodebug;
    TYPE* points;
    TYPE* weights;
    TYPE* points_ordered;
    TYPE* weights_ordered;
    size_t num_points;
    size_t ncrit;
    int num_terms;
    TYPE theta;
    TYPE theta2;
    size_t num_samples;

    size_t num_multipoles;
    size_t num_nodes;
    size_t num_spharm_terms;

    TYPE* acc;
    TYPE* pot;
    TYPE* inner_factors_real;
    TYPE* inner_factors_imag;
    TYPE* outer_factors_real;
    TYPE* outer_factors_imag;

    std::vector<TYPE> x;
    std::vector<TYPE> y;
    std::vector<TYPE> z;
    std::vector<TYPE> w;
    std::vector<TYPE> ax;
    std::vector<TYPE> ay;
    std::vector<TYPE> az;
    std::vector<TYPE> p;

    std::vector<TYPE> M_array_real;
    std::vector<TYPE> M_array_imag;
    std::vector<TYPE> L_array_real;
    std::vector<TYPE> L_array_imag;
    std::vector<t_node> node_array;
};

