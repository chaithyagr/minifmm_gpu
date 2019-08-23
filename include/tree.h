#pragma once

#include <stdlib.h>

#include "type.h"
#include "params.h"

/* typedef struct t_node
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
} t_node; */

inline int is_leaf(t_node* node) { return node->num_children == 0; }

void build_tree(t_fmm_params* params);

void free_tree(t_fmm_params* params);

void tree_to_result(t_fmm_params* params);

inline
t_node* get_node(t_fmm_params* params, size_t node) 
{ 
    return &params->node_array[node]; 
}
