#include <stdio.h>
#include <stdlib.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "params.h"
#include "rng.h"
#include "parse_args.h"
#include "traversal.h"
#include "type.h"
#include "kernels.h"
#include "timer.h"
#include "verify.h"
#include "initialise.h"
#include "finalise.h"
#include "util.h"
namespace py = pybind11;


py::array_t<TYPE, py::array::c_style | py::array::forcecast> minifmm_pybind(py::size_t num_points,
 py::array_t<TYPE, py::array::c_style | py::array::forcecast> locations,
 py::array_t<TYPE, py::array::c_style | py::array::forcecast> weights, int num_terms, int bin_size, TYPE theta, int num_samples, bool dodebug)
{
    py::buffer_info loc = locations.request();
    TYPE *points = (TYPE *) loc.ptr;
    py::buffer_info wt_buffer = weights.request();
    TYPE *wt = (TYPE *) wt_buffer.ptr;
    // Output Array
    py::array_t<TYPE, py::array::c_style | py::array::forcecast> output({num_points*5});
    py::buffer_info out = output.request();
    TYPE *ptrOut = (TYPE *) out.ptr;
    struct t_fmm_params params = t_fmm_params();
    params.dodebug=dodebug;
    init_with_known_data(&params, wt, points, num_points, num_terms, bin_size, theta, num_samples);
    if(dodebug)
    {
        print_args(&params);
    }
    perform_fmm(&params);
    if(dodebug)
    {
        printf("Some random values : Pot[0] = %f Pot[1] = %f Pot[num_points] = %f Pot[30] = %f\n",params.pot[0],params.pot[1],params.pot[params.num_points],params.pot[30]);
        TYPE a_err=0, p_err=0;
        verify(&params, &a_err, &p_err);
        printf("force err.     = %e\n", a_err);
        printf("potential err. = %e\n", p_err);
    }
    for(int l = 0 ; l < num_points*5 ; l++)
    {
        if(l<num_points)
            ptrOut[l] = params.pot[l];
        else if(l<2*num_points)
            ptrOut[l] = params.reorder[l-num_points];
        else if(l<4*num_points)
            ptrOut[l] = params.points_ordered[l-2*num_points];
    }
    finalise(&params);
    return output;
}

PYBIND11_MODULE(minifmm, m)
{
    m.doc() = "Calculate the Repulsive Force gradient with MiniFMM on GPU"; // optional
    m.def("minifmm_pybind", &minifmm_pybind, py::return_value_policy::take_ownership);
}