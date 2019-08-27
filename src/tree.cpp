#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <immintrin.h>

#include "type.h"
#include "util.h"
#include "tree.h"

void bound_box(TYPE* mins, TYPE* maxs, TYPE* points, size_t num_points)
{
    for (int i = 0; i < 3; ++i) mins[i] = TYPE_MAX;
    for (int i = 0; i < 3; ++i) maxs[i] = -TYPE_MAX;

    for (size_t i = 0; i < num_points; ++i)
    {
        for (int d = 0; d < 3; ++d)
        {
            mins[d] = MIN(mins[d], points[d*num_points+i]);
            maxs[d] = MAX(maxs[d], points[d*num_points+i]);
        }
    }
}

static size_t point_head = 0;

size_t construct_tree(TYPE cx, TYPE cy, TYPE cz, TYPE r, size_t start, size_t end, size_t ncrit, TYPE* points, TYPE* weights,
        TYPE* points_ordered, TYPE* weights_ordered, int *reorder, int *reorder_ordered, TYPE* acc, TYPE* pot,
        int points_switched, size_t num_points, int level, t_fmm_params* params)
{
    for (int i = 0; i < params->num_multipoles; ++i)
    {
        params->M_array_real.emplace_back(0.0);
        params->M_array_imag.emplace_back(0.0);
        params->L_array_real.emplace_back(0.0);
        params->L_array_imag.emplace_back(0.0);
    }

    params->node_array.push_back(t_node());    
    size_t parent_idx = params->node_array.size()-1;
    t_node* parent = &(params->node_array.back());

    parent->mult_idx = params->num_nodes*params->num_multipoles;
    parent->offset = params->node_array.size()-1;

    if (end - start <= ncrit)
    {
        parent->point_idx = point_head;
        point_head += (end-start);
    } else parent->point_idx = 0;

    params->num_nodes++;

    parent->center[0] = cx;
    parent->center[1] = cy;
    parent->center[2] = cz;
    parent->rad = r;
    parent->level = level;
    parent->num_points = end - start;
    parent->num_children = 0;

    // if leaf node
    if (end - start <= ncrit)
    {
        if (!points_switched)
        {
            for (size_t i = start; i < end; ++i)
            {
                for (int d = 0; d < 3; ++d) points_ordered[d*num_points+i] = points[d*num_points+i];
                weights_ordered[i] = weights[i];
                reorder_ordered[i] = reorder[i];
            }
        }
    }
    else
    {
        size_t num_points_per_oct[8] = {0};
        for (size_t i = start; i < end; ++i)
        {
            // courtesy of exaFMM
            int oct = (points[0*num_points+i] > parent->center[0]) + 
                ((points[1*num_points+i] > parent->center[1]) << 1) + 
                ((points[2*num_points+i] > parent->center[2]) << 2);
            num_points_per_oct[oct]++;           
        }

        size_t oct_pointers[8] = {start};
        size_t oct_pointers_copy[8] = {start};
        for (int i = 1; i < 8; ++i) oct_pointers[i] = oct_pointers[i-1] + num_points_per_oct[i-1];
        for (int i = 0; i < 8; ++i) oct_pointers_copy[i] = oct_pointers[i];

        for (size_t j = start; j < end; ++j)
        {
            int oct = (points[0*num_points+j] > parent->center[0]) + 
                ((points[1*num_points+j] > parent->center[1]) << 1) + 
                ((points[2*num_points+j] > parent->center[2]) << 2);

            size_t i = oct_pointers_copy[oct];
            for (int d = 0; d < 3; ++d) 
            {
                points_ordered[d*num_points+i] = points[d*num_points+j];
            }
            weights_ordered[i] = weights[j];
            reorder_ordered[i] = reorder[j];
            oct_pointers_copy[oct]++;
        }

        TYPE new_r = parent->rad/(TYPE_TWO);
        for (int i = 0; i < 8; ++i)
        {
            if (num_points_per_oct[i])
            {
                TYPE ncx, ncy, ncz;
                ncx = ((i >> 0) & 1) ? (parent->center[0] + new_r) : (parent->center[0] - new_r);
                ncy = ((i >> 1) & 1) ? (parent->center[1] + new_r) : (parent->center[1] - new_r);
                ncz = ((i >> 2) & 1) ? (parent->center[2] + new_r) : (parent->center[2] - new_r);
    
                size_t child = construct_tree(ncx, ncy, ncz, new_r, oct_pointers[i], oct_pointers[i]+num_points_per_oct[i], 
                    ncrit, points_ordered, weights_ordered, points, weights, reorder_ordered, reorder, acc, pot, !points_switched, num_points, level+1, params);
                parent = &params->node_array[parent_idx];
                parent->child[parent->num_children++] = child;
            }
        }
    }
    return parent->offset;
}

//int count_desc_nodes(t_node* root)
//{
//    int count = root->num_children;
//    for (int i = 0; i < root->num_children; ++i)
//    {
//        count += count_desc_nodes(root->child[i]);
//    }
//    root->num_desc_nodes = count;
//    return count;
//}

void build_tree(t_fmm_params* params)
{
    //init_node_array(params);
    point_head=0;
    TYPE mins[3], maxs[3];
    bound_box(mins, maxs, params->points, params->num_points);

    //printf("bound box --- \n");
    //for (int d = 0; d < 3; ++d) printf("%f %f\n", mins[d], maxs[d]);

    TYPE max_rad = TYPE_ZERO;
    TYPE center[3];
    for (int d = 0; d < 3; ++d) 
    {
        center[d] = (maxs[d] + mins[d])/(TYPE_TWO);
        TYPE rad = (maxs[d] - mins[d])/(TYPE_TWO);
        max_rad = MAX(rad, max_rad);
    }
    // need to add EPS for points that lie on border of node
    max_rad += TYPE_EPS;

    params->root = construct_tree(center[0], center[1], center[2], max_rad, 0, params->num_points, params->ncrit, params->points, params->weights,
            params->points_ordered, params->weights_ordered, params->reorder, params->reorder_ordered, params->acc, params->pot, 0, params->num_points, 0, params);

    //printf("Tree has %zu nodes\n", params->num_nodes);
    //count_desc_nodes(params->root);
    if(params->dodebug)
    {
        printf("num nodes = %zu\n", params->num_nodes);
        fflush(stdout);
    }

    size_t np = params->num_points;
    TYPE* ptr = params->points_ordered;
    params->x.assign(ptr, ptr+np); ptr += np;
    params->y.assign(ptr, ptr+np); ptr += np;
    params->z.assign(ptr, ptr+np); 
    params->w.assign(params->weights_ordered, params->weights_ordered+np); 
    params->ax.resize(np);
    params->ay.resize(np);
    params->az.resize(np);
    params->p.resize(np);
    for (size_t i = 0; i < np; ++i)
    {
        params->ax[i] = params->ay[i] = params->az[i] = params->p[i] = 0.0;
    }
}

//void free_tree_core(t_node* node)
//{
//    for (size_t i = 0; i < node->num_children; ++i) free_tree_core(node->child[i]);
//    free_node(node);
//}

//void free_tree(t_fmm_params* params)
//{
//    free_tree_core(params->root);
//    free(__node_array);
//}

